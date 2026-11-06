// CUDA core, see cpp/include/bwpm/solver_cuda.hpp.
//
// Device complex type: thrust::complex<double> (CCCL, part of the CUDA toolkit). The interface
// formulas are a line-by-line copy of cpp/src/interface.cpp; keep both in step.
#include "bwpm/solver_cuda.hpp"

#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace bwpm {

namespace {

using dcplx = thrust::complex<double>;

void check(cudaError_t e, const char* what) {
  if (e != cudaSuccess)
    throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e));
}

void check(cufftResult e, const char* what) {
  if (e != CUFFT_SUCCESS)
    throw std::runtime_error(std::string("cuFFT error ") + std::to_string(static_cast<int>(e)) +
                             " in " + what);
}

template <class T>
struct DeviceBuffer {
  T* ptr = nullptr;
  std::size_t size = 0;
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t n) : size(n) {
    if (n) check(cudaMalloc(&ptr, n * sizeof(T)), "cudaMalloc");
  }
  ~DeviceBuffer() {
    if (ptr) cudaFree(ptr);
  }
  DeviceBuffer(DeviceBuffer&& o) noexcept : ptr(o.ptr), size(o.size) {
    o.ptr = nullptr;
    o.size = 0;
  }
  DeviceBuffer& operator=(DeviceBuffer&& o) noexcept {
    std::swap(ptr, o.ptr);
    std::swap(size, o.size);
    return *this;
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
  void upload(const T* host, std::size_t n) {
    check(cudaMemcpy(ptr, host, n * sizeof(T), cudaMemcpyHostToDevice), "upload");
  }
  void download(T* host, std::size_t n) const {
    check(cudaMemcpy(host, ptr, n * sizeof(T), cudaMemcpyDeviceToHost), "download");
  }
};

// ---- interface formulas (copy of cpp/src/interface.cpp, device version) --------------------

__device__ dcplx kz_of(dcplx index, double k0, double kp2) {
  dcplx kz = thrust::sqrt(index * index * (k0 * k0) - kp2);
  if (kz.imag() < 0.0 || (kz.imag() == 0.0 && kz.real() < 0.0)) kz = -kz;
  const double tiny = 1e-12 * k0;
  if (thrust::abs(kz) < tiny) kz = dcplx(0.0, tiny);
  return kz;
}

struct DFresnel {
  dcplx r_s, r_p, t_s, t_p, kz_in, kz_out;
};

__device__ DFresnel fresnel(double kx, double ky, dcplx n1, dcplx n2, double k0) {
  const double kp2 = kx * kx + ky * ky;
  const dcplx kz1 = kz_of(n1, k0, kp2);
  const dcplx kz2 = kz_of(n2, k0, kp2);
  const dcplx q1 = kz1 / k0, q2 = kz2 / k0;
  const dcplx den_p = n2 * n2 * q1 + n1 * n1 * q2;
  return {(q1 - q2) / (q1 + q2), (n2 * n2 * q1 - n1 * n1 * q2) / den_p, 2.0 * q1 / (q1 + q2),
          2.0 * n1 * n2 * q1 / den_p, kz1, kz2};
}

// kind: 0 fwd_t, 1 fwd_r, 2 bwd_t, 3 bwd_r (table 1 of wende2023), row-major 2×2 into m[4]
__device__ void operator_matrix(int kind, double kx, double ky, dcplx n1, dcplx n2, double k0,
                                double dz, dcplx m[4]) {
  const DFresnel c = fresnel(kx, ky, n1, n2, k0);
  const double kp2 = kx * kx + ky * ky;
  const dcplx im(0.0, 1.0);
  const bool reflection = (kind == 1 || kind == 3);
  if (kp2 == 0.0) {
    const dcplx d = reflection ? c.r_s : c.t_s;
    m[0] = d;
    m[1] = 0.0;
    m[2] = 0.0;
    m[3] = d;
  } else if (!reflection) {
    const dcplx kr = n1 * k0, ks = n2 * k0;
    const dcplx gamma = c.kz_in * c.kz_out + c.kz_out * kp2 / c.kz_in;
    const dcplx scale = 1.0 / (kp2 * kr * ks);
    const dcplx off = (c.t_p * gamma - c.t_s * kr * ks) * kx * ky * scale;
    m[0] = (c.t_p * kx * kx * gamma + c.t_s * ky * ky * kr * ks) * scale;
    m[1] = off;
    m[2] = off;
    m[3] = (c.t_p * ky * ky * gamma + c.t_s * kx * kx * kr * ks) * scale;
  } else {
    const dcplx kr = n1 * k0;
    const dcplx gamma = c.kz_in * c.kz_in + kp2;
    const dcplx scale = 1.0 / (kp2 * kr * kr);
    const dcplx off = -(c.r_p * gamma + c.r_s * kr * kr) * kx * ky * scale;
    m[0] = (-c.r_p * kx * kx * gamma + c.r_s * ky * ky * kr * kr) * scale;
    m[1] = off;
    m[2] = off;
    m[3] = (-c.r_p * ky * ky * gamma + c.r_s * kx * kx * kr * kr) * scale;
  }
  dcplx phase(1.0, 0.0);
  if (kind == 0) phase = thrust::exp(im * c.kz_out * dz);        // i = 1: k_z,s
  if (kind == 2) phase = thrust::exp(im * c.kz_in * dz);         // i = 3: k_z,r
  if (kind == 3) phase = thrust::exp(2.0 * im * c.kz_in * dz);   // i = 4: 2 k_z,r
  for (int j = 0; j < 4; ++j) m[j] *= phase;
}

// 2π·numpy.fft.fftfreq rounding: 2π·(m·(1/(n·d))), as cpp/src/solver.cpp
__device__ double fftfreq(long long i, long long n, double val) {
  const long long m = i - (i < (n + 1) / 2 ? 0LL : n);
  return 2.0 * 3.141592653589793 * (static_cast<double>(m) * val);
}

// ---- kernels ---------------------------------------------------------------------------------

__global__ void k_apply_matrix(const dcplx* spectrum, dcplx* work, std::size_t ny, std::size_t nx,
                               double val_x, double val_y, int kind, dcplx n1, dcplx n2,
                               double k0, double dz) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  dcplx m[4];
  operator_matrix(kind, kx, ky, n1, n2, k0, dz, m);
  const dcplx ex = spectrum[2 * q], ey = spectrum[2 * q + 1];
  work[2 * q] = m[0] * ex + m[1] * ey;
  work[2 * q + 1] = m[2] * ex + m[3] * ey;
}

__global__ void k_scatter(const dcplx* work, dcplx* out, const std::int32_t* lr,
                          const std::int32_t* ls, std::int32_t r, std::int32_t s, double norm,
                          std::size_t n) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= n) return;
  if (lr[q] == r && ls[q] == s) {
    out[2 * q] = work[2 * q] * norm;
    out[2 * q + 1] = work[2 * q + 1] * norm;
  }
}

__global__ void k_mask(dcplx* out, const double* mask, std::size_t n) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= n) return;
  out[2 * q] *= mask[q];
  out[2 * q + 1] *= mask[q];
}

__global__ void k_gather(const dcplx* field, const int* index, dcplx* values, std::size_t m) {
  const std::size_t j = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (j >= m) return;
  values[2 * j] = field[2 * static_cast<std::size_t>(index[j])];
  values[2 * j + 1] = field[2 * static_cast<std::size_t>(index[j]) + 1];
}

__global__ void k_add(dcplx* field, const int* index, const dcplx* values, std::size_t m) {
  const std::size_t j = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (j >= m) return;
  field[2 * static_cast<std::size_t>(index[j])] += values[2 * j];
  field[2 * static_cast<std::size_t>(index[j]) + 1] += values[2 * j + 1];
}

constexpr int kThreads = 256;
unsigned blocks(std::size_t n) { return static_cast<unsigned>((n + kThreads - 1) / kThreads); }

}  // namespace

int cuda_device_count() {
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess) {
    cudaGetLastError();  // clear the error state
    return 0;
  }
  return n;
}

struct CudaSolver::Impl {
  Grid g;
  std::vector<cplx> indices;
  std::vector<std::int32_t> labels;
  std::vector<std::int32_t> before;
  std::size_t nz;
  cplx index_before;
  bool has_mask;
  DeviceBuffer<double> mask;
  DeviceBuffer<std::int32_t> lab_r, lab_s;
  DeviceBuffer<dcplx> spectrum, work;
  cufftHandle plan{};
  bool plan_ok = false;
  std::vector<std::set<std::pair<std::int32_t, std::int32_t>>> pairs_fwd;  // pairs of step p

  std::size_t n() const { return g.ny * g.nx; }
  const std::int32_t* plane(std::size_t p) const {
    return p == 0 ? before.data() : labels.data() + (p - 1) * n();
  }
  cplx index_of(std::int32_t l) const {
    return l < 0 ? index_before : indices[static_cast<std::size_t>(l)];
  }

  // out = Σ_pairs I_r,s F⁻¹{M F{in}} (reflections skip equal indices), mask on transmissions
  void apply(int kind, const DeviceBuffer<dcplx>& in, std::size_t p_r, std::size_t p_s,
             DeviceBuffer<dcplx>& out) {
    const std::size_t nn = n();
    lab_r.upload(plane(p_r), nn);
    lab_s.upload(plane(p_s), nn);
    check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(in.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(spectrum.ptr), CUFFT_FORWARD),
          "forward FFT");
    check(cudaMemset(out.ptr, 0, 2 * nn * sizeof(dcplx)), "memset");
    const bool reflection = (kind == 1 || kind == 3);
    const std::size_t step = std::min(p_r, p_s);  // pairs are symmetric in the two planes
    const auto& pairs = pairs_fwd[step];
    const double val_x = 1.0 / (static_cast<double>(g.nx) * g.dx);
    const double val_y = 1.0 / (static_cast<double>(g.ny) * g.dy);
    const double norm = 1.0 / static_cast<double>(nn);
    for (const auto& [a, b] : pairs) {
      // pairs_fwd holds (plane step, plane step + 1); backward steps swap the roles
      const std::int32_t r = (p_r <= p_s) ? a : b;
      const std::int32_t s = (p_r <= p_s) ? b : a;
      const cplx nr = index_of(r), ns = index_of(s);
      if (reflection && nr == ns) continue;
      k_apply_matrix<<<blocks(nn), kThreads>>>(spectrum.ptr, work.ptr, g.ny, g.nx, val_x, val_y,
                                               kind, dcplx(nr.real(), nr.imag()),
                                               dcplx(ns.real(), ns.imag()), g.k0, g.dz);
      check(cudaGetLastError(), "k_apply_matrix");
      check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(work.ptr),
                         reinterpret_cast<cufftDoubleComplex*>(work.ptr), CUFFT_INVERSE),
            "inverse FFT");
      k_scatter<<<blocks(nn), kThreads>>>(work.ptr, out.ptr, lab_r.ptr, lab_s.ptr, r, s, norm, nn);
      check(cudaGetLastError(), "k_scatter");
    }
    if (!reflection && has_mask) {
      k_mask<<<blocks(nn), kThreads>>>(out.ptr, mask.ptr, nn);
      check(cudaGetLastError(), "k_mask");
    }
  }
};

CudaSolver::CudaSolver(Grid grid, std::vector<cplx> indices, std::vector<std::int32_t> labels,
                       std::size_t nz, cplx index_before, std::vector<double> mask, int device)
    : impl_(new Impl) {
  Impl& s = *impl_;
  try {
    check(cudaSetDevice(device), "cudaSetDevice");
    s.g = grid;
    s.indices = std::move(indices);
    s.labels = std::move(labels);
    s.nz = nz;
    s.index_before = index_before;
    const std::size_t nn = s.n();
    if (s.labels.size() != nz * nn) throw std::invalid_argument("labels must have shape (nz, ny, nx)");
    if (!mask.empty() && mask.size() != nn) throw std::invalid_argument("mask must have shape (ny, nx)");
    for (auto l : s.labels)
      if (l < 0 || static_cast<std::size_t>(l) >= s.indices.size())
        throw std::invalid_argument("labels must index into `indices`");
    s.before.assign(nn, -1);
    s.has_mask = !mask.empty();
    if (s.has_mask) {
      s.mask = DeviceBuffer<double>(nn);
      s.mask.upload(mask.data(), nn);
    }
    s.lab_r = DeviceBuffer<std::int32_t>(nn);
    s.lab_s = DeviceBuffer<std::int32_t>(nn);
    s.spectrum = DeviceBuffer<dcplx>(2 * nn);
    s.work = DeviceBuffer<dcplx>(2 * nn);
    // 2D FFT over (ny, nx) of both interleaved components: stride 2, distance 1, batch 2
    int dims[2] = {static_cast<int>(grid.ny), static_cast<int>(grid.nx)};
    check(cufftPlanMany(&s.plan, 2, dims, dims, 2, 1, dims, 2, 1, CUFFT_Z2Z, 2), "cufftPlanMany");
    s.plan_ok = true;
    s.pairs_fwd.resize(nz);
    for (std::size_t p = 0; p < nz; ++p) {
      const std::int32_t* a = s.plane(p);
      const std::int32_t* b = s.plane(p + 1);
      for (std::size_t q = 0; q < nn; ++q) s.pairs_fwd[p].emplace(a[q], b[q]);
    }
  } catch (...) {
    if (impl_->plan_ok) cufftDestroy(impl_->plan);
    delete impl_;
    impl_ = nullptr;
    throw;
  }
}

CudaSolver::~CudaSolver() {
  if (impl_) {
    if (impl_->plan_ok) cufftDestroy(impl_->plan);
    delete impl_;
  }
}

RunResult CudaSolver::run(const std::vector<cplx>& field_in, int iterations,
                          const std::vector<std::size_t>& keep) {
  Impl& s = *impl_;
  const std::size_t nn = s.n();
  if (field_in.size() != 2 * nn) throw std::invalid_argument("field must have shape (ny, nx, 2)");
  const std::set<std::size_t> keep_set(keep.begin(), keep.end());

  // interface pixels per step p (index changes between plane p and plane p + 1), on the device
  struct Sources {
    DeviceBuffer<int> index;
    DeviceBuffer<dcplx> value;
    std::size_t m = 0;
  };
  std::vector<std::vector<int>> interface_host(s.nz);
  std::vector<DeviceBuffer<int>> interface_dev(s.nz);
  for (std::size_t p = 0; p < s.nz; ++p) {
    const std::int32_t* a = s.plane(p);
    const std::int32_t* b = s.plane(p + 1);
    for (std::size_t q = 0; q < nn; ++q)
      if (s.index_of(a[q]) != s.index_of(b[q])) interface_host[p].push_back(static_cast<int>(q));
    if (!interface_host[p].empty()) {
      interface_dev[p] = DeviceBuffer<int>(interface_host[p].size());
      interface_dev[p].upload(interface_host[p].data(), interface_host[p].size());
    }
  }

  auto gather = [&](const DeviceBuffer<dcplx>& field, std::size_t p) {
    Sources src;
    src.m = interface_host[p].size();
    src.value = DeviceBuffer<dcplx>(2 * src.m);
    k_gather<<<blocks(src.m), kThreads>>>(field.ptr, interface_dev[p].ptr, src.value.ptr, src.m);
    check(cudaGetLastError(), "k_gather");
    return src;
  };
  auto add = [&](DeviceBuffer<dcplx>& field, const Sources& src, std::size_t p) {
    k_add<<<blocks(src.m), kThreads>>>(field.ptr, interface_dev[p].ptr, src.value.ptr, src.m);
    check(cudaGetLastError(), "k_add");
  };
  auto to_host = [&](const DeviceBuffer<dcplx>& field) {
    std::vector<cplx> host(2 * nn);
    field.download(reinterpret_cast<dcplx*>(host.data()), 2 * nn);
    return host;
  };

  DeviceBuffer<dcplx> input(2 * nn), e(2 * nn), next(2 * nn), reflected(2 * nn);
  input.upload(reinterpret_cast<const dcplx*>(field_in.data()), 2 * nn);
  RunResult res;
  // sources: from backward reflections at plane p (interface p − 1) and forward reflections at
  // plane p (interface p); both keyed by the interface index used to gather them
  std::map<std::size_t, Sources> src_fwd, src_bwd;
  for (int it = 0; it < std::max(iterations, 1); ++it) {
    src_bwd.clear();
    check(cudaMemcpy(e.ptr, input.ptr, 2 * nn * sizeof(dcplx), cudaMemcpyDeviceToDevice), "copy");
    if (auto f = src_fwd.find(0); f != src_fwd.end()) add(e, f->second, 0);
    if (keep_set.count(0)) res.forward[0] = to_host(e);
    for (std::size_t p = 0; p < s.nz; ++p) {
      if (iterations > 0 && !interface_host[p].empty()) {
        s.apply(1, e, p, p + 1, next);
        src_bwd[p] = gather(next, p);
      }
      s.apply(0, e, p, p + 1, next);
      std::swap(e, next);
      // backward reflection at interface p was stored for plane p + 1
      if (auto f = src_fwd.find(p + 1); f != src_fwd.end()) add(e, f->second, p);
      if (keep_set.count(p + 1)) res.forward[p + 1] = to_host(e);
    }
    res.transmitted = to_host(e);
    if (iterations == 0) break;
    src_fwd.clear();
    check(cudaMemset(e.ptr, 0, 2 * nn * sizeof(dcplx)), "memset");
    if (keep_set.count(s.nz)) res.backward[s.nz] = to_host(e);
    for (std::size_t p = s.nz; p > 0; --p) {
      if (!interface_host[p - 1].empty()) {
        s.apply(3, e, p, p - 1, next);
        src_fwd[p] = gather(next, p - 1);
      }
      s.apply(2, e, p, p - 1, next);
      std::swap(e, next);
      if (auto f = src_bwd.find(p - 1); f != src_bwd.end()) add(e, f->second, p - 1);
      if (keep_set.count(p - 1)) res.backward[p - 1] = to_host(e);
    }
    res.reflected = to_host(e);
    res.transmitted_iter.push_back(res.transmitted);
    res.reflected_iter.push_back(res.reflected);
    std::size_t stored = 0;
    for (const auto& [p, src] : src_bwd) stored += src.m;
    for (const auto& [p, src] : src_fwd) stored += src.m;
    res.stored_values = std::max(res.stored_values, stored);
  }
  if (res.reflected.empty()) res.reflected.assign(2 * nn, cplx(0.0, 0.0));
  check(cudaDeviceSynchronize(), "synchronize");
  return res;
}

}  // namespace bwpm
