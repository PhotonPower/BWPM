// CPU core, see cpp/include/bwpm/solver.hpp.
#include "bwpm/solver.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "pocketfft_hdronly.h"

namespace bwpm {

namespace {
double fftfreq(std::size_t i, std::size_t n, double d) {
  // 2π·numpy.fft.fftfreq(n, d) with the same rounding: 2π·(m·(1/(n·d))), m = 0 … ⌈n/2⌉−1,
  // −⌊n/2⌋ … −1. Grid points within rounding of |k⊥| = n·k₀ make k_z cancellation-sensitive,
  // so the rounding of k must match the reference exactly.
  const auto m = static_cast<long long>(i) - (i < (n + 1) / 2 ? 0LL : static_cast<long long>(n));
  const double val = 1.0 / (static_cast<double>(n) * d);
  return 2.0 * std::numbers::pi * (static_cast<double>(m) * val);
}
}  // namespace

Solver::Solver(Grid grid, std::vector<cplx> indices, std::vector<std::int32_t> labels,
               std::size_t nz, cplx index_before, std::vector<double> mask, int threads,
               std::vector<cplx> fluctuation)
    : g_(grid), indices_(std::move(indices)), labels_(std::move(labels)), nz_(nz),
      index_before_(index_before), mask_(std::move(mask)), threads_(std::max(threads, 1)),
      fluct_(std::move(fluctuation)) {
  const std::size_t n = g_.ny * g_.nx;
  if (labels_.size() != nz_ * n) throw std::invalid_argument("labels must have shape (nz, ny, nx)");
  if (!mask_.empty() && mask_.size() != n) throw std::invalid_argument("mask must have shape (ny, nx)");
  for (auto l : labels_)
    if (l < 0 || static_cast<std::size_t>(l) >= indices_.size())
      throw std::invalid_argument("labels must index into `indices`");
  before_.assign(n, -1);
  kx_.resize(n);
  ky_.resize(n);
  for (std::size_t j = 0; j < g_.ny; ++j)
    for (std::size_t i = 0; i < g_.nx; ++i) {
      kx_[j * g_.nx + i] = fftfreq(i, g_.nx, g_.dx);
      ky_[j * g_.nx + i] = fftfreq(j, g_.ny, g_.dy);
    }
  has_host_.assign(nz_, 0);
  host_.assign(nz_, cplx(0.0, 0.0));
  if (!fluct_.empty()) {
    if (fluct_.size() != nz_ * n)
      throw std::invalid_argument("fluctuation must have the shape of labels (nz, ny, nx)");
    for (std::size_t j = 0; j < nz_; ++j)
      for (std::size_t q = 0; q < n; ++q) {
        if (fluct_[j * n + q] == cplx(0.0, 0.0)) continue;
        const cplx idx = indices_[static_cast<std::size_t>(labels_[j * n + q])];
        if (has_host_[j] && host_[j] != idx)
          throw std::invalid_argument("slab " + std::to_string(j) +
                                      ": fluctuation must lie in one host index");
        has_host_[j] = 1;
        host_[j] = idx;
      }
  }
}

const std::int32_t* Solver::plane_labels(std::size_t p) const {
  return p == 0 ? before_.data() : labels_.data() + (p - 1) * g_.ny * g_.nx;
}

cplx Solver::index_of(std::int32_t label) const {
  return label < 0 ? index_before_ : indices_[static_cast<std::size_t>(label)];
}

const std::vector<cplx>& Solver::op(Kind kind, std::int32_t r, std::int32_t s) {
  const Key key{kind, r, s};
  auto it = cache_.find(key);
  if (it != cache_.end()) return it->second;
  const std::size_t n = g_.ny * g_.nx;
  std::vector<cplx> m(4 * n);
  const cplx nr = index_of(r), ns = index_of(s);
  const cplx im(0.0, 1.0);
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    const double kx = kx_[q], ky = ky_[q];
    const double kp2 = kx * kx + ky * ky;
    Matrix2 a;
    cplx phase(1.0, 0.0);
    switch (kind) {  // table 1 of wende2023
      case FwdT:
        a = transmission_matrix(kx, ky, nr, ns, g_.k0);
        phase = std::exp(im * longitudinal_wavenumber(ns, g_.k0, kp2) * g_.dz);
        break;
      case FwdR:
        a = reflection_matrix(kx, ky, nr, ns, g_.k0);
        break;
      case BwdT:
        a = transmission_matrix(kx, ky, nr, ns, g_.k0);
        phase = std::exp(im * longitudinal_wavenumber(nr, g_.k0, kp2) * g_.dz);
        break;
      default:  // BwdR
        a = reflection_matrix(kx, ky, nr, ns, g_.k0);
        phase = std::exp(2.0 * im * longitudinal_wavenumber(nr, g_.k0, kp2) * g_.dz);
    }
    for (int j = 0; j < 4; ++j) m[4 * q + j] = a[j] * phase;
  }
  return cache_.emplace(key, std::move(m)).first->second;
}

void Solver::fft(std::vector<cplx>& data, bool forward) const {
  const pocketfft::shape_t shape{g_.ny, g_.nx, 2};
  const auto c = static_cast<std::ptrdiff_t>(sizeof(cplx));
  const pocketfft::stride_t stride{static_cast<std::ptrdiff_t>(2 * g_.nx) * c, 2 * c, c};
  const double fct = forward ? 1.0 : 1.0 / static_cast<double>(g_.ny * g_.nx);
  pocketfft::c2c<double>(shape, stride, stride, {0, 1}, forward, data.data(), data.data(), fct,
                         static_cast<std::size_t>(threads_));
}

void Solver::fft_components(std::vector<cplx>& data, std::size_t ncomp, bool forward) const {
  const pocketfft::shape_t shape{g_.ny, g_.nx, ncomp};
  const auto c = static_cast<std::ptrdiff_t>(sizeof(cplx));
  const auto nc = static_cast<std::ptrdiff_t>(ncomp);
  const pocketfft::stride_t stride{static_cast<std::ptrdiff_t>(g_.nx) * nc * c, nc * c, c};
  const double fct = forward ? 1.0 : 1.0 / static_cast<double>(g_.ny * g_.nx);
  pocketfft::c2c<double>(shape, stride, stride, {0, 1}, forward, data.data(), data.data(), fct,
                         static_cast<std::size_t>(threads_));
}

namespace {
// GMRES with restarts for A x = b (complex, matrix-free), modified Gram–Schmidt and Givens
// rotations; x holds the start value and the solution. Returns false if not converged.
template <class Op>
bool gmres(const Op& apply, const std::vector<cplx>& b, std::vector<cplx>& x, double tol,
           int restart, int max_restarts) {
  const std::size_t n = b.size();
  auto dot = [n](const std::vector<cplx>& u, const std::vector<cplx>& v) {
    cplx s(0.0, 0.0);
    for (std::size_t i = 0; i < n; ++i) s += std::conj(u[i]) * v[i];
    return s;
  };
  auto norm = [&](const std::vector<cplx>& u) { return std::sqrt(std::real(dot(u, u))); };
  const double bnorm = norm(b);
  if (bnorm == 0.0) {
    std::fill(x.begin(), x.end(), cplx(0.0, 0.0));
    return true;
  }
  const auto m = static_cast<std::size_t>(restart);
  std::vector<std::vector<cplx>> basis(m + 1, std::vector<cplx>(n));
  std::vector<cplx> h((m + 1) * m), g(m + 1), cs(m), sn(m), w(n);
  for (int outer = 0; outer < max_restarts; ++outer) {
    apply(x, w);
    for (std::size_t i = 0; i < n; ++i) basis[0][i] = b[i] - w[i];
    const double beta = norm(basis[0]);
    if (beta <= tol * bnorm) return true;
    for (std::size_t i = 0; i < n; ++i) basis[0][i] /= beta;
    std::fill(g.begin(), g.end(), cplx(0.0, 0.0));
    g[0] = beta;
    std::size_t used = 0;
    bool done = false;
    for (std::size_t j = 0; j < m; ++j) {
      apply(basis[j], w);
      for (std::size_t i = 0; i <= j; ++i) {
        const cplx hij = dot(basis[i], w);
        h[i * m + j] = hij;
        for (std::size_t q = 0; q < n; ++q) w[q] -= hij * basis[i][q];
      }
      const double hnext = norm(w);
      if (hnext > 0.0)
        for (std::size_t q = 0; q < n; ++q) basis[j + 1][q] = w[q] / hnext;
      for (std::size_t i = 0; i < j; ++i) {  // previous rotations on column j
        const cplx a = h[i * m + j], c = h[(i + 1) * m + j];
        h[i * m + j] = cs[i] * a + sn[i] * c;
        h[(i + 1) * m + j] = -std::conj(sn[i]) * a + cs[i] * c;
      }
      const cplx a = h[j * m + j];
      const double r = std::sqrt(std::norm(a) + hnext * hnext);
      if (std::abs(a) == 0.0) {
        cs[j] = 0.0;
        sn[j] = 1.0;
      } else {
        cs[j] = std::abs(a) / r;
        sn[j] = (a / std::abs(a)) * hnext / r;
      }
      h[j * m + j] = cs[j] * a + sn[j] * hnext;
      g[j + 1] = -std::conj(sn[j]) * g[j];
      g[j] = cs[j] * g[j];
      used = j + 1;
      if (std::abs(g[j + 1]) <= tol * bnorm || hnext == 0.0) {
        done = true;
        break;
      }
    }
    std::vector<cplx> y(used);  // back substitution
    for (std::size_t i = used; i-- > 0;) {
      cplx s = g[i];
      for (std::size_t k = i + 1; k < used; ++k) s -= h[i * m + k] * y[k];
      y[i] = s / h[i * m + i];
    }
    for (std::size_t i = 0; i < used; ++i)
      for (std::size_t q = 0; q < n; ++q) x[q] += y[i] * basis[i][q];
    if (done) {
      apply(x, w);
      double res = 0.0;
      for (std::size_t q = 0; q < n; ++q) res += std::norm(b[q] - w[q]);
      if (std::sqrt(res) <= tol * bnorm * 10.0) return true;
    }
  }
  return false;
}

constexpr double kGmresTolerance = 1e-12;  // §12, G2 (as GMRES_TOLERANCE in bwpm/fluctuation.py)

// Gauss–Legendre points ζ/dz of the slab integral, weights 1/2 (§11, F1, #53)
const double kGauss[2] = {(1.0 - 1.0 / std::sqrt(3.0)) / 2.0, (1.0 + 1.0 / std::sqrt(3.0)) / 2.0};
}  // namespace

std::vector<cplx> Solver::born_sources(const std::vector<cplx>& in, int direction,
                                       std::size_t slab) {
  // Ŝ_m = F{V·E_m}, E_m = F⁻¹{e^(i k_z s_m)(Ê⊥, Ê_z)}, s_m = ζ_m (d = +1) or dz − ζ_m (d = −1),
  // Ê_z = −d·k⊥·Ê⊥/k_z (wende2023 eq. (5)), V = 2ñ₀k₀²δn (E11); propagating, non-grazing only.
  const std::size_t n = g_.ny * g_.nx;
  const cplx host = host_[slab];
  const cplx host2k2 = host * host * g_.k0 * g_.k0;  // as ñ²k₀² in longitudinal_wavenumber
  const double kabs2 = std::abs(host2k2);
  const cplx im(0.0, 1.0);
  std::vector<cplx> spectrum = in;
  fft(spectrum, true);
  std::vector<cplx> full(3 * n), kz(n);
  std::vector<char> prop(n);
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    const double kx = kx_[q], ky = ky_[q];
    const double kp2 = kx * kx + ky * ky;
    kz[q] = longitudinal_wavenumber(host, g_.k0, kp2);
    prop[q] = kp2 < kabs2 && std::abs(host2k2 - kp2) > kGrazingTolerance * kabs2;  // §3 (#25)
    const cplx ex = spectrum[2 * q], ey = spectrum[2 * q + 1];
    full[3 * q] = ex;
    full[3 * q + 1] = ey;
    full[3 * q + 2] = prop[q] ? -static_cast<double>(direction) * (kx * ex + ky * ey) / kz[q]
                              : cplx(0.0, 0.0);
  }
  const cplx* delta = fluct_.data() + slab * n;
  const cplx scale = 2.0 * host * g_.k0 * g_.k0;
  std::vector<cplx> sources(6 * n), work(3 * n);
  for (int m = 0; m < 2; ++m) {
    const double distance = (direction > 0 ? kGauss[m] : 1.0 - kGauss[m]) * g_.dz;
#pragma omp parallel for num_threads(threads_) schedule(static)
    for (long long q = 0; q < static_cast<long long>(n); ++q) {
      const cplx phase = prop[q] ? std::exp(im * kz[q] * distance) : cplx(0.0, 0.0);
      for (int c = 0; c < 3; ++c) work[3 * q + c] = full[3 * q + c] * phase;
    }
    fft_components(work, 3, false);
#pragma omp parallel for num_threads(threads_) schedule(static)
    for (long long q = 0; q < static_cast<long long>(n); ++q) {
      const cplx v = scale * delta[q];
      for (int c = 0; c < 3; ++c) work[3 * q + c] = v * work[3 * q + c];
    }
    fft_components(work, 3, true);
    std::copy(work.begin(), work.end(), sources.begin() + static_cast<std::ptrdiff_t>(3 * n * m));
  }
  return sources;
}

void Solver::born_operator(const std::vector<cplx>& spectrum, int direction, std::size_t slab,
                           std::vector<cplx>& out) {
  // G_d·S = i dz/(2k_z)·[P(k_d)·F{V·e}]⊥, e = F⁻¹(Ŝ⊥, Ŝ_z), Ŝ_z = −d·k⊥·Ŝ⊥/k_z (§12)
  const std::size_t n = g_.ny * g_.nx;
  const cplx host = host_[slab];
  const cplx k = host * g_.k0;
  const cplx k2 = k * k;
  const cplx host2k2 = host * host * g_.k0 * g_.k0;
  const double kabs2 = std::abs(host2k2);
  const cplx im(0.0, 1.0);
  const double d = static_cast<double>(direction);
  std::vector<cplx> work(3 * n), kz(n);
  std::vector<char> prop(n);
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    const double kx = kx_[q], ky = ky_[q];
    const double kp2 = kx * kx + ky * ky;
    kz[q] = longitudinal_wavenumber(host, g_.k0, kp2);
    prop[q] = kp2 < kabs2 && std::abs(host2k2 - kp2) > kGrazingTolerance * kabs2;  // §3 (#25)
    const cplx ex = prop[q] ? spectrum[2 * q] : cplx(0.0, 0.0);
    const cplx ey = prop[q] ? spectrum[2 * q + 1] : cplx(0.0, 0.0);
    work[3 * q] = ex;
    work[3 * q + 1] = ey;
    work[3 * q + 2] = prop[q] ? -d * (kx * ex + ky * ey) / kz[q] : cplx(0.0, 0.0);
  }
  fft_components(work, 3, false);
  const cplx* delta = fluct_.data() + slab * n;
  const cplx scale = 2.0 * host * g_.k0 * g_.k0;
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    const cplx v = scale * delta[q];
    for (int c = 0; c < 3; ++c) work[3 * q + c] = v * work[3 * q + c];
  }
  fft_components(work, 3, true);
  out.resize(2 * n);
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    if (!prop[q]) {
      out[2 * q] = 0.0;
      out[2 * q + 1] = 0.0;
      continue;
    }
    const double kx = kx_[q], ky = ky_[q];
    const cplx green = im * g_.dz / (2.0 * kz[q]);
    const cplx s0 = work[3 * q], s1 = work[3 * q + 1], s2 = work[3 * q + 2];
    const cplx longitudinal = (kx * s0 + ky * s1 + d * kz[q] * s2) / k2;
    out[2 * q] = green * (s0 - kx * longitudinal);
    out[2 * q + 1] = green * (s1 - ky * longitudinal);
  }
}

std::vector<cplx> Solver::born_step(const std::vector<cplx>& in, int direction,
                                    std::size_t slab) {
  // x = P_½·Ê, (1 − G_d/2)·y′ = G_d·x (GMRES, start G_d·x), ΔÊ = P_½^(±)·y′ (§12)
  const std::size_t n = g_.ny * g_.nx;
  const cplx host = host_[slab];
  const cplx host2k2 = host * host * g_.k0 * g_.k0;
  const double kabs2 = std::abs(host2k2);
  const cplx im(0.0, 1.0);
  std::vector<cplx> x = in;
  fft(x, true);
  std::vector<cplx> half(n);
  for (std::size_t q = 0; q < n; ++q) {
    const double kp2 = kx_[q] * kx_[q] + ky_[q] * ky_[q];
    const bool prop = kp2 < kabs2 && std::abs(host2k2 - kp2) > kGrazingTolerance * kabs2;
    half[q] = prop ? std::exp(im * longitudinal_wavenumber(host, g_.k0, kp2) * (g_.dz / 2.0))
                   : cplx(0.0, 0.0);
    x[2 * q] *= half[q];
    x[2 * q + 1] *= half[q];
  }
  std::vector<cplx> rhs;
  born_operator(x, direction, slab, rhs);
  std::vector<cplx> y = rhs;
  std::vector<cplx> tmp;
  auto apply = [&](const std::vector<cplx>& v, std::vector<cplx>& av) {
    born_operator(v, direction, slab, tmp);
    av.resize(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) av[i] = v[i] - 0.5 * tmp[i];
  };
  if (!gmres(apply, rhs, y, kGmresTolerance, 30, 100))
    throw std::runtime_error("GMRES did not converge in the Born step");
  for (std::size_t q = 0; q < n; ++q) {
    const cplx f = direction > 0 ? half[q] : (half[q] == cplx(0.0, 0.0) ? half[q] : 1.0 / half[q]);
    y[2 * q] *= f;
    y[2 * q + 1] *= f;
  }
  fft(y, false);
  return y;
}

std::vector<cplx> Solver::born_emission(const std::vector<cplx>& sources, int direction,
                                        std::size_t slab) {
  // ΔE⊥ = F⁻¹{Σ_m ½·e^(i d′k_z(dz − ζ_m))·i dz/(2k_z)·[P(k_d′)·Ŝ_m]⊥}, k_d′ = (k_x, k_y, d′k_z).
  const std::size_t n = g_.ny * g_.nx;
  const cplx host = host_[slab];
  const cplx k = host * g_.k0;
  const cplx k2 = k * k;
  const cplx host2k2 = host * host * g_.k0 * g_.k0;
  const double kabs2 = std::abs(host2k2);
  const cplx im(0.0, 1.0);
  const double d = static_cast<double>(direction);
  std::vector<cplx> total(2 * n, cplx(0.0, 0.0));
#pragma omp parallel for num_threads(threads_) schedule(static)
  for (long long q = 0; q < static_cast<long long>(n); ++q) {
    const double kx = kx_[q], ky = ky_[q];
    const double kp2 = kx * kx + ky * ky;
    if (!(kp2 < kabs2 && std::abs(host2k2 - kp2) > kGrazingTolerance * kabs2)) continue;
    const cplx kz = longitudinal_wavenumber(host, g_.k0, kp2);
    const cplx signed_kz = d * kz;
    for (int m = 0; m < 2; ++m) {
      const cplx green = 0.5 * std::exp(im * signed_kz * ((1.0 - kGauss[m]) * g_.dz)) *
                         (im * g_.dz / (2.0 * kz));
      const cplx* s = sources.data() + 3 * n * static_cast<std::size_t>(m) + 3 * q;
      const cplx longitudinal = (kx * s[0] + ky * s[1] + signed_kz * s[2]) / k2;
      total[2 * q] += green * (s[0] - kx * longitudinal);
      total[2 * q + 1] += green * (s[1] - ky * longitudinal);
    }
  }
  fft(total, false);
  return total;
}

void Solver::apply(Kind kind, const std::vector<cplx>& in, std::size_t p_r, std::size_t p_s,
                   std::vector<cplx>& out) {
  const std::size_t n = g_.ny * g_.nx;
  const std::int32_t* lr = plane_labels(p_r);
  const std::int32_t* ls = plane_labels(p_s);
  std::set<std::pair<std::int32_t, std::int32_t>> pairs;
  for (std::size_t q = 0; q < n; ++q) pairs.emplace(lr[q], ls[q]);
  std::vector<cplx> spectrum = in;
  fft(spectrum, true);
  out.assign(2 * n, cplx(0.0, 0.0));
  std::vector<cplx> work(2 * n);
  const bool reflection = (kind == FwdR || kind == BwdR);
  for (const auto& [r, s] : pairs) {
    if (reflection && index_of(r) == index_of(s)) continue;
    const std::vector<cplx>& m = op(kind, r, s);
#pragma omp parallel for num_threads(threads_) schedule(static)
    for (long long q = 0; q < static_cast<long long>(n); ++q) {
      const cplx ex = spectrum[2 * q], ey = spectrum[2 * q + 1];
      work[2 * q] = m[4 * q] * ex + m[4 * q + 1] * ey;
      work[2 * q + 1] = m[4 * q + 2] * ex + m[4 * q + 3] * ey;
    }
    fft(work, false);
    for (std::size_t q = 0; q < n; ++q)
      if (lr[q] == r && ls[q] == s) {
        out[2 * q] = work[2 * q];
        out[2 * q + 1] = work[2 * q + 1];
      }
  }
  if (!reflection && !mask_.empty())
    for (std::size_t q = 0; q < n; ++q) {
      out[2 * q] *= mask_[q];
      out[2 * q + 1] *= mask_[q];
    }
}

namespace {
struct Sparse {
  std::vector<std::size_t> index;
  std::vector<cplx> value;  // 2 per index
};

Sparse store(const std::vector<cplx>& v, const std::vector<char>& where) {
  Sparse s;
  for (std::size_t q = 0; q < where.size(); ++q)
    if (where[q]) {
      s.index.push_back(q);
      s.value.push_back(v[2 * q]);
      s.value.push_back(v[2 * q + 1]);
    }
  return s;
}

void add(std::vector<cplx>& v, const Sparse& s) {
  for (std::size_t j = 0; j < s.index.size(); ++j) {
    v[2 * s.index[j]] += s.value[2 * j];
    v[2 * s.index[j] + 1] += s.value[2 * j + 1];
  }
}
}  // namespace

RunResult Solver::run(const std::vector<cplx>& field_in, int iterations,
                      const std::vector<std::size_t>& keep) {
  const std::size_t n = g_.ny * g_.nx;
  if (field_in.size() != 2 * n) throw std::invalid_argument("field must have shape (ny, nx, 2)");
  const std::set<std::size_t> keep_set(keep.begin(), keep.end());
  // interface pixels between plane p (left medium) and plane p + 1 (slab p)
  std::vector<std::vector<char>> interface(nz_);
  for (std::size_t p = 0; p < nz_; ++p) {
    interface[p].assign(n, 0);
    const std::int32_t* a = plane_labels(p);
    const std::int32_t* b = plane_labels(p + 1);
    for (std::size_t q = 0; q < n; ++q) interface[p][q] = index_of(a[q]) != index_of(b[q]);
  }
  std::vector<bool> any(nz_);
  for (std::size_t p = 0; p < nz_; ++p)
    any[p] = std::find(interface[p].begin(), interface[p].end(), 1) != interface[p].end();

  RunResult res;
  std::map<std::size_t, Sparse> src_fwd, src_bwd;
  // Born sources of δn as full planes (F3): R1 of the forward pass, R3 of the backward pass
  std::map<std::size_t, std::vector<cplx>> born_fwd, born_bwd;
  auto add_plane = [](std::vector<cplx>& v, const std::vector<cplx>& w) {
    for (std::size_t i = 0; i < v.size(); ++i) v[i] += w[i];
  };
  std::vector<cplx> e, next;
  for (int it = 0; it < std::max(iterations, 1); ++it) {
    src_bwd.clear();
    born_bwd.clear();
    e = field_in;
    if (auto f = src_fwd.find(0); f != src_fwd.end()) add(e, f->second);
    if (keep_set.count(0)) res.forward[0] = e;
    for (std::size_t p = 0; p < nz_; ++p) {
      if (iterations > 0 && any[p]) {
        apply(FwdR, e, p, p + 1, next);
        src_bwd[p] = store(next, interface[p]);
      }
      apply(FwdT, e, p, p + 1, next);
      if (has_host_[p]) {
        const std::vector<cplx> increment = born_step(e, +1, p);
        for (std::size_t q = 0; q < n; ++q) {
          const double m = mask_.empty() ? 1.0 : mask_[q];
          next[2 * q] += increment[2 * q] * m;
          next[2 * q + 1] += increment[2 * q + 1] * m;
        }
        if (iterations > 0) born_bwd[p + 1] = born_emission(born_sources(e, +1, p), -1, p);  // R1
      }
      e.swap(next);
      if (auto f = src_fwd.find(p + 1); f != src_fwd.end()) add(e, f->second);
      if (auto f = born_fwd.find(p + 1); f != born_fwd.end()) add_plane(e, f->second);
      if (keep_set.count(p + 1)) res.forward[p + 1] = e;
    }
    res.transmitted = e;
    if (iterations == 0) break;
    src_fwd.clear();
    born_fwd.clear();
    e.assign(2 * n, cplx(0.0, 0.0));
    if (auto f = src_bwd.find(nz_); f != src_bwd.end()) add(e, f->second);
    if (auto f = born_bwd.find(nz_); f != born_bwd.end()) add_plane(e, f->second);
    if (keep_set.count(nz_)) res.backward[nz_] = e;
    for (std::size_t p = nz_; p > 0; --p) {
      if (has_host_[p - 1]) {
        born_fwd[p] = born_emission(born_sources(e, -1, p - 1), +1, p - 1);  // R3
        add_plane(e, born_step(e, -1, p - 1));  // R2 (Cayley, §12), referred to plane p
      }
      if (any[p - 1]) {
        apply(BwdR, e, p, p - 1, next);
        src_fwd[p] = store(next, interface[p - 1]);
      }
      apply(BwdT, e, p, p - 1, next);
      e.swap(next);
      if (auto f = src_bwd.find(p - 1); f != src_bwd.end()) add(e, f->second);
      if (auto f = born_bwd.find(p - 1); f != born_bwd.end()) add_plane(e, f->second);
      if (keep_set.count(p - 1)) res.backward[p - 1] = e;
    }
    res.reflected = e;
    res.transmitted_iter.push_back(res.transmitted);
    res.reflected_iter.push_back(res.reflected);
    std::size_t stored = 0;
    for (const auto& [p, s] : src_bwd) stored += s.index.size();
    for (const auto& [p, s] : src_fwd) stored += s.index.size();
    res.stored_values = std::max(res.stored_values, stored);
  }
  if (res.reflected.empty()) res.reflected.assign(2 * n, cplx(0.0, 0.0));
  return res;
}

}  // namespace bwpm
