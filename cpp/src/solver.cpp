// CPU core, see cpp/include/bwpm/solver.hpp.
#include "bwpm/solver.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <stdexcept>
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
               std::size_t nz, cplx index_before, std::vector<double> mask, int threads)
    : g_(grid), indices_(std::move(indices)), labels_(std::move(labels)), nz_(nz),
      index_before_(index_before), mask_(std::move(mask)), threads_(std::max(threads, 1)) {
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
  std::vector<cplx> e, next;
  for (int it = 0; it < std::max(iterations, 1); ++it) {
    src_bwd.clear();
    e = field_in;
    if (auto f = src_fwd.find(0); f != src_fwd.end()) add(e, f->second);
    if (keep_set.count(0)) res.forward[0] = e;
    for (std::size_t p = 0; p < nz_; ++p) {
      if (iterations > 0 && any[p]) {
        apply(FwdR, e, p, p + 1, next);
        src_bwd[p] = store(next, interface[p]);
      }
      apply(FwdT, e, p, p + 1, next);
      e.swap(next);
      if (auto f = src_fwd.find(p + 1); f != src_fwd.end()) add(e, f->second);
      if (keep_set.count(p + 1)) res.forward[p + 1] = e;
    }
    res.transmitted = e;
    if (iterations == 0) break;
    src_fwd.clear();
    e.assign(2 * n, cplx(0.0, 0.0));
    if (auto f = src_bwd.find(nz_); f != src_bwd.end()) add(e, f->second);
    if (keep_set.count(nz_)) res.backward[nz_] = e;
    for (std::size_t p = nz_; p > 0; --p) {
      if (any[p - 1]) {
        apply(BwdR, e, p, p - 1, next);
        src_fwd[p] = store(next, interface[p - 1]);
      }
      apply(BwdT, e, p, p - 1, next);
      e.swap(next);
      if (auto f = src_bwd.find(p - 1); f != src_bwd.end()) add(e, f->second);
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
