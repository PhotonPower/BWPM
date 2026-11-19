// CUDA core, see cpp/include/bwpm/solver_cuda.hpp.
//
// Device complex type: thrust::complex<double> (CCCL, part of the CUDA toolkit). The interface
// formulas are a line-by-line copy of cpp/src/interface.cpp; keep both in step.
#include "bwpm/solver_cuda.hpp"

#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>

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
  const dcplx n2k2 = index * index * (k0 * k0);
  dcplx q = n2k2 - kp2;
  if (thrust::abs(q) <= kGrazingTolerance * thrust::abs(n2k2)) q = 0.0;  // grazing rule (#25)
  dcplx kz = thrust::sqrt(q);
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

// ---- Born increment for weak fluctuations (docs/spezifikation.md §10, as cpp/src/solver.cpp) --

__device__ bool born_propagating(dcplx host2k2, double kp2) {
  const double kabs2 = thrust::abs(host2k2);
  return kp2 < kabs2 && thrust::abs(host2k2 - kp2) > kGrazingTolerance * kabs2;  // §3 (#25)
}

// Gauss–Legendre points ζ/dz of the slab integral, weights 1/2 (§11, F1, #53)
__device__ double gauss_point(int m) {
  const double r = 0.5773502691896257;  // 1/√3, as (1 ∓ 1/√3)/2 in cpp/src/solver.cpp
  return m == 0 ? (1.0 - r) / 2.0 : (1.0 + r) / 2.0;
}

// spectrum: F{E⊥} -> full3 = (Ê_x, Ê_y, Ê_z), Ê_z = −d·k⊥·Ê⊥/k_z, propagating only
__global__ void k_born_full(const dcplx* spectrum, dcplx* full3, std::size_t ny, std::size_t nx,
                            double val_x, double val_y, dcplx host, double k0, int direction) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  const dcplx kz = kz_of(host, k0, kp2);
  const bool propagating = born_propagating(host * host * (k0 * k0), kp2);
  const dcplx ex = spectrum[2 * q], ey = spectrum[2 * q + 1];
  full3[3 * q] = ex;
  full3[3 * q + 1] = ey;
  full3[3 * q + 2] =
      propagating ? -static_cast<double>(direction) * (kx * ex + ky * ey) / kz : dcplx(0.0, 0.0);
}

// work3 = e^(i k_z s)·full3 (s = distance to the Gauss point), propagating only
__global__ void k_born_phase(const dcplx* full3, dcplx* work3, std::size_t ny, std::size_t nx,
                             double val_x, double val_y, dcplx host, double k0, double distance) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  const dcplx kz = kz_of(host, k0, kp2);
  const dcplx im(0.0, 1.0);
  const dcplx phase = born_propagating(host * host * (k0 * k0), kp2)
                          ? thrust::exp(im * kz * distance)
                          : dcplx(0.0, 0.0);
  for (int c = 0; c < 3; ++c) work3[3 * q + c] = full3[3 * q + c] * phase;
}

// work3 = V·E_m with V = 2ñ₀k₀²δn (E11); `norm` completes the unnormalised inverse FFT
__global__ void k_born_potential(dcplx* work3, const dcplx* delta, dcplx scale, double norm,
                                 std::size_t n) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= n) return;
  const dcplx v = scale * delta[q] * norm;
  for (int c = 0; c < 3; ++c) work3[3 * q + c] = v * work3[3 * q + c];
}

// sources (2 Gauss points × 3) -> total = Σ_m ½·e^(i d′k_z(dz − ζ_m))·i dz/(2k_z)·[P(k_d′)·Ŝ_m]⊥
__global__ void k_born_emit(const dcplx* sources, dcplx* total, std::size_t ny, std::size_t nx,
                            double val_x, double val_y, dcplx host, double k0, double dz,
                            int direction) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  const std::size_t n = ny * nx;
  if (q >= n) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  total[2 * q] = 0.0;
  total[2 * q + 1] = 0.0;
  if (!born_propagating(host * host * (k0 * k0), kp2)) return;
  const dcplx kz = kz_of(host, k0, kp2);
  const dcplx signed_kz = static_cast<double>(direction) * kz;
  const dcplx im(0.0, 1.0);
  const dcplx k = host * k0;
  const dcplx k2 = k * k;
  for (int m = 0; m < 2; ++m) {
    const dcplx green =
        0.5 * thrust::exp(im * signed_kz * ((1.0 - gauss_point(m)) * dz)) * (im * dz / (2.0 * kz));
    const dcplx* src = sources + 3 * n * static_cast<std::size_t>(m) + 3 * q;
    const dcplx longitudinal = (kx * src[0] + ky * src[1] + signed_kz * src[2]) / k2;
    total[2 * q] += green * (src[0] - kx * longitudinal);
    total[2 * q + 1] += green * (src[1] - ky * longitudinal);
  }
}

// ---- flux-conserving Born step (§12, B9b): G_d, phases and vector operations for GMRES ------

// G_d, first half: work3 = (Ŝ_x, Ŝ_y, Ŝ_z) on propagating waves, Ŝ_z = −d·k⊥·Ŝ⊥/k_z
__global__ void k_g_prepare(const dcplx* spectrum, dcplx* work3, std::size_t ny, std::size_t nx,
                            double val_x, double val_y, dcplx host, double k0, int direction) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  if (!born_propagating(host * host * (k0 * k0), kp2)) {
    work3[3 * q] = work3[3 * q + 1] = work3[3 * q + 2] = 0.0;
    return;
  }
  const dcplx kz = kz_of(host, k0, kp2);
  const dcplx ex = spectrum[2 * q], ey = spectrum[2 * q + 1];
  work3[3 * q] = ex;
  work3[3 * q + 1] = ey;
  work3[3 * q + 2] = -static_cast<double>(direction) * (kx * ex + ky * ey) / kz;
}

// G_d, second half: out = i dz/(2k_z)·[P(k_d)·Ŝ]⊥ on propagating waves
__global__ void k_g_finish(const dcplx* work3, dcplx* out, std::size_t ny, std::size_t nx,
                           double val_x, double val_y, dcplx host, double k0, double dz,
                           int direction) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  if (!born_propagating(host * host * (k0 * k0), kp2)) {
    out[2 * q] = out[2 * q + 1] = 0.0;
    return;
  }
  const dcplx kz = kz_of(host, k0, kp2);
  const dcplx im(0.0, 1.0);
  const dcplx k = host * k0;
  const dcplx green = im * dz / (2.0 * kz);
  const dcplx s0 = work3[3 * q], s1 = work3[3 * q + 1], s2 = work3[3 * q + 2];
  const dcplx longitudinal = (kx * s0 + ky * s1 + static_cast<double>(direction) * kz * s2) / (k * k);
  out[2 * q] = green * (s0 - kx * longitudinal);
  out[2 * q + 1] = green * (s1 - ky * longitudinal);
}

// data *= e^(i·sign·k_z dz/2) on propagating waves, 0 elsewhere (P_½ and its inverse)
__global__ void k_half(dcplx* data, std::size_t ny, std::size_t nx, double val_x, double val_y,
                       dcplx host, double k0, double dz, double sign) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= ny * nx) return;
  const double kx = fftfreq(static_cast<long long>(q % nx), static_cast<long long>(nx), val_x);
  const double ky = fftfreq(static_cast<long long>(q / nx), static_cast<long long>(ny), val_y);
  const double kp2 = kx * kx + ky * ky;
  dcplx f(0.0, 0.0);
  if (born_propagating(host * host * (k0 * k0), kp2))
    f = thrust::exp(dcplx(0.0, sign) * kz_of(host, k0, kp2) * (dz / 2.0));
  data[2 * q] *= f;
  data[2 * q + 1] *= f;
}

__global__ void k_axpby(dcplx* y, const dcplx* x, dcplx a, dcplx b, std::size_t m) {
  // y = a·x + b·y
  const std::size_t i = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (i < m) y[i] = a * x[i] + b * y[i];
}

struct ConjMul {
  __host__ __device__ dcplx operator()(const dcplx& u, const dcplx& v) const {
    return thrust::conj(u) * v;
  }
};

// out += inc·norm·mask (mask may be null)
__global__ void k_born_add(dcplx* out, const dcplx* inc, const double* mask, double norm,
                           std::size_t n) {
  const std::size_t q = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
  if (q >= n) return;
  const double m = mask ? mask[q] * norm : norm;
  out[2 * q] += inc[2 * q] * m;
  out[2 * q + 1] += inc[2 * q + 1] * m;
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
  // weak fluctuations (B7b): δn on the host, uploaded per slab
  std::vector<cplx> fluct;
  std::vector<char> has_host;
  std::vector<cplx> host;
  DeviceBuffer<dcplx> delta, full3, work3, sources, inc;
  cufftHandle plan3{};
  bool plan3_ok = false;
  // GMRES (B9b): Krylov basis and work vectors of 2·n entries
  std::vector<DeviceBuffer<dcplx>> basis;
  DeviceBuffer<dcplx> g_x, g_rhs, g_y, g_w, g_tmp;
  std::vector<std::set<std::pair<std::int32_t, std::int32_t>>> pairs_fwd;  // pairs of step p

  std::size_t n() const { return g.ny * g.nx; }
  const std::int32_t* plane(std::size_t p) const {
    return p == 0 ? before.data() : labels.data() + (p - 1) * n();
  }
  cplx index_of(std::int32_t l) const {
    return l < 0 ? index_before : indices[static_cast<std::size_t>(l)];
  }

  // Born sources of slab `slab` at the Gauss points (§11) into `sources`; `spectrum` must hold
  // F{E⊥} of the incident field of `direction` (+1 left slab edge, −1 right slab edge)
  void born_sources(std::size_t slab, int direction) {
    const std::size_t nn = n();
    const double val_x = 1.0 / (static_cast<double>(g.nx) * g.dx);
    const double val_y = 1.0 / (static_cast<double>(g.ny) * g.dy);
    const double norm = 1.0 / static_cast<double>(nn);
    const dcplx h(host[slab].real(), host[slab].imag());
    delta.upload(reinterpret_cast<const dcplx*>(fluct.data()) + slab * nn, nn);
    k_born_full<<<blocks(nn), kThreads>>>(spectrum.ptr, full3.ptr, g.ny, g.nx, val_x, val_y, h,
                                          g.k0, direction);
    check(cudaGetLastError(), "k_born_full");
    const double r = 0.5773502691896257;
    const double gauss[2] = {(1.0 - r) / 2.0, (1.0 + r) / 2.0};
    for (int m = 0; m < 2; ++m) {
      const double distance = (direction > 0 ? gauss[m] : 1.0 - gauss[m]) * g.dz;
      k_born_phase<<<blocks(nn), kThreads>>>(full3.ptr, work3.ptr, g.ny, g.nx, val_x, val_y, h,
                                             g.k0, distance);
      check(cudaGetLastError(), "k_born_phase");
      check(cufftExecZ2Z(plan3, reinterpret_cast<cufftDoubleComplex*>(work3.ptr),
                         reinterpret_cast<cufftDoubleComplex*>(work3.ptr), CUFFT_INVERSE),
            "inverse FFT (Born)");
      k_born_potential<<<blocks(nn), kThreads>>>(work3.ptr, delta.ptr, 2.0 * h * (g.k0 * g.k0),
                                                 norm, nn);
      check(cudaGetLastError(), "k_born_potential");
      check(cufftExecZ2Z(plan3, reinterpret_cast<cufftDoubleComplex*>(work3.ptr),
                         reinterpret_cast<cufftDoubleComplex*>(sources.ptr + 3 * nn * m),
                         CUFFT_FORWARD),
            "forward FFT (Born)");
    }
  }

  // out += (field radiated by `sources` in `direction` at the right slab edge)·mask
  void add_born_emission(std::size_t slab, int direction, DeviceBuffer<dcplx>& out,
                         bool masked) {
    const std::size_t nn = n();
    const double val_x = 1.0 / (static_cast<double>(g.nx) * g.dx);
    const double val_y = 1.0 / (static_cast<double>(g.ny) * g.dy);
    const double norm = 1.0 / static_cast<double>(nn);
    const dcplx h(host[slab].real(), host[slab].imag());
    k_born_emit<<<blocks(nn), kThreads>>>(sources.ptr, inc.ptr, g.ny, g.nx, val_x, val_y, h,
                                          g.k0, g.dz, direction);
    check(cudaGetLastError(), "k_born_emit");
    check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(inc.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(inc.ptr), CUFFT_INVERSE),
          "inverse FFT (Born emission)");
    k_born_add<<<blocks(nn), kThreads>>>(out.ptr, inc.ptr,
                                         masked && has_mask ? mask.ptr : nullptr, norm, nn);
    check(cudaGetLastError(), "k_born_add");
  }

  // out = G_d·in on spectra (2·n), §12; δn of the slab must be uploaded
  void born_operator(const dcplx* in, dcplx* out, std::size_t slab, int direction) {
    const std::size_t nn = n();
    const double val_x = 1.0 / (static_cast<double>(g.nx) * g.dx);
    const double val_y = 1.0 / (static_cast<double>(g.ny) * g.dy);
    const double norm = 1.0 / static_cast<double>(nn);
    const dcplx h(host[slab].real(), host[slab].imag());
    k_g_prepare<<<blocks(nn), kThreads>>>(in, work3.ptr, g.ny, g.nx, val_x, val_y, h, g.k0,
                                          direction);
    check(cudaGetLastError(), "k_g_prepare");
    check(cufftExecZ2Z(plan3, reinterpret_cast<cufftDoubleComplex*>(work3.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(work3.ptr), CUFFT_INVERSE),
          "inverse FFT (G)");
    k_born_potential<<<blocks(nn), kThreads>>>(work3.ptr, delta.ptr, 2.0 * h * (g.k0 * g.k0),
                                               norm, nn);
    check(cudaGetLastError(), "k_born_potential (G)");
    check(cufftExecZ2Z(plan3, reinterpret_cast<cufftDoubleComplex*>(work3.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(work3.ptr), CUFFT_FORWARD),
          "forward FFT (G)");
    k_g_finish<<<blocks(nn), kThreads>>>(work3.ptr, out, g.ny, g.nx, val_x, val_y, h, g.k0, g.dz,
                                         direction);
    check(cudaGetLastError(), "k_g_finish");
  }

  dcplx dot(const dcplx* u, const dcplx* v, std::size_t m) {
    return thrust::inner_product(thrust::device, u, u + m, v, dcplx(0.0, 0.0),
                                 thrust::plus<dcplx>(), ConjMul());
  }

  // GMRES(30) on the device for (1 − G/2)·y = rhs, start y = rhs, as cpp/src/solver.cpp
  void gmres_cayley(std::size_t slab, int direction) {
    const std::size_t m2 = 2 * n();
    const int restart = 30;
    const double tol = 1e-12;  // §12, G2
    if (basis.empty())
      for (int i = 0; i <= restart; ++i) basis.emplace_back(m2);
    auto apply = [&](const dcplx* v, dcplx* av) {  // av = v − G v/2
      born_operator(v, g_tmp.ptr, slab, direction);
      check(cudaMemcpy(av, v, m2 * sizeof(dcplx), cudaMemcpyDeviceToDevice), "copy");
      k_axpby<<<blocks(m2), kThreads>>>(av, g_tmp.ptr, dcplx(-0.5, 0.0), dcplx(1.0, 0.0), m2);
      check(cudaGetLastError(), "k_axpby");
    };
    auto norm2 = [&](const dcplx* u) { return std::sqrt(dot(u, u, m2).real()); };
    check(cudaMemcpy(g_y.ptr, g_rhs.ptr, m2 * sizeof(dcplx), cudaMemcpyDeviceToDevice), "copy");
    const double bnorm = norm2(g_rhs.ptr);
    if (bnorm == 0.0) return;
    const auto m = static_cast<std::size_t>(restart);
    std::vector<dcplx> h((m + 1) * m), gv(m + 1), cs(m), sn(m);
    for (int outer = 0; outer < 100; ++outer) {
      apply(g_y.ptr, g_w.ptr);
      check(cudaMemcpy(basis[0].ptr, g_rhs.ptr, m2 * sizeof(dcplx), cudaMemcpyDeviceToDevice),
            "copy");
      k_axpby<<<blocks(m2), kThreads>>>(basis[0].ptr, g_w.ptr, dcplx(-1.0, 0.0),
                                        dcplx(1.0, 0.0), m2);
      const double beta = norm2(basis[0].ptr);
      if (beta <= tol * bnorm) return;
      k_axpby<<<blocks(m2), kThreads>>>(basis[0].ptr, basis[0].ptr, dcplx(0.0, 0.0),
                                        dcplx(1.0 / beta, 0.0), m2);
      std::fill(gv.begin(), gv.end(), dcplx(0.0, 0.0));
      gv[0] = beta;
      std::size_t used = 0;
      bool done = false;
      for (std::size_t j = 0; j < m; ++j) {
        apply(basis[j].ptr, g_w.ptr);
        for (std::size_t i = 0; i <= j; ++i) {
          const dcplx hij = dot(basis[i].ptr, g_w.ptr, m2);
          h[i * m + j] = hij;
          k_axpby<<<blocks(m2), kThreads>>>(g_w.ptr, basis[i].ptr, -hij, dcplx(1.0, 0.0), m2);
        }
        const double hnext = norm2(g_w.ptr);
        if (hnext > 0.0) {
          check(cudaMemcpy(basis[j + 1].ptr, g_w.ptr, m2 * sizeof(dcplx), cudaMemcpyDeviceToDevice),
                "copy");
          k_axpby<<<blocks(m2), kThreads>>>(basis[j + 1].ptr, basis[j + 1].ptr, dcplx(0.0, 0.0),
                                            dcplx(1.0 / hnext, 0.0), m2);
        }
        for (std::size_t i = 0; i < j; ++i) {
          const dcplx a = h[i * m + j], c = h[(i + 1) * m + j];
          h[i * m + j] = cs[i] * a + sn[i] * c;
          h[(i + 1) * m + j] = -thrust::conj(sn[i]) * a + cs[i] * c;
        }
        const dcplx a = h[j * m + j];
        const double r = std::sqrt(thrust::norm(a) + hnext * hnext);
        if (thrust::abs(a) == 0.0) {
          cs[j] = 0.0;
          sn[j] = 1.0;
        } else {
          cs[j] = thrust::abs(a) / r;
          sn[j] = (a / thrust::abs(a)) * hnext / r;
        }
        h[j * m + j] = cs[j] * a + sn[j] * hnext;
        gv[j + 1] = -thrust::conj(sn[j]) * gv[j];
        gv[j] = cs[j] * gv[j];
        used = j + 1;
        if (thrust::abs(gv[j + 1]) <= tol * bnorm || hnext == 0.0) {
          done = true;
          break;
        }
      }
      std::vector<dcplx> y(used);
      for (std::size_t i = used; i-- > 0;) {
        dcplx s = gv[i];
        for (std::size_t k = i + 1; k < used; ++k) s -= h[i * m + k] * y[k];
        y[i] = s / h[i * m + i];
      }
      for (std::size_t i = 0; i < used; ++i)
        k_axpby<<<blocks(m2), kThreads>>>(g_y.ptr, basis[i].ptr, y[i], dcplx(1.0, 0.0), m2);
      check(cudaGetLastError(), "GMRES update");
      if (done) {
        apply(g_y.ptr, g_w.ptr);
        check(cudaMemcpy(g_tmp.ptr, g_rhs.ptr, m2 * sizeof(dcplx), cudaMemcpyDeviceToDevice),
              "copy");
        k_axpby<<<blocks(m2), kThreads>>>(g_tmp.ptr, g_w.ptr, dcplx(-1.0, 0.0), dcplx(1.0, 0.0),
                                          m2);
        if (norm2(g_tmp.ptr) <= tol * bnorm * 10.0) return;
      }
    }
    throw std::runtime_error("GMRES did not converge in the Born step (CUDA)");
  }

  // target += (flux-conserving Born increment of slab `slab` for `field` of `direction`)·mask
  void add_born_step(const DeviceBuffer<dcplx>& field, std::size_t slab, int direction,
                     DeviceBuffer<dcplx>& target, bool masked) {
    const std::size_t nn = n();
    const double val_x = 1.0 / (static_cast<double>(g.nx) * g.dx);
    const double val_y = 1.0 / (static_cast<double>(g.ny) * g.dy);
    const double norm = 1.0 / static_cast<double>(nn);
    const dcplx h(host[slab].real(), host[slab].imag());
    delta.upload(reinterpret_cast<const dcplx*>(fluct.data()) + slab * nn, nn);
    check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(field.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(g_x.ptr), CUFFT_FORWARD),
          "forward FFT (Born step)");
    k_half<<<blocks(nn), kThreads>>>(g_x.ptr, g.ny, g.nx, val_x, val_y, h, g.k0, g.dz, 1.0);
    check(cudaGetLastError(), "k_half");
    born_operator(g_x.ptr, g_rhs.ptr, slab, direction);
    gmres_cayley(slab, direction);
    k_half<<<blocks(nn), kThreads>>>(g_y.ptr, g.ny, g.nx, val_x, val_y, h, g.k0, g.dz,
                                     direction > 0 ? 1.0 : -1.0);
    check(cudaGetLastError(), "k_half");
    check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(g_y.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(g_y.ptr), CUFFT_INVERSE),
          "inverse FFT (Born step)");
    k_born_add<<<blocks(nn), kThreads>>>(target.ptr, g_y.ptr, masked && has_mask ? mask.ptr : nullptr,
                                         norm, nn);
    check(cudaGetLastError(), "k_born_add (step)");
  }

  // F{field} into `spectrum`
  void transform(const DeviceBuffer<dcplx>& field) {
    check(cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(field.ptr),
                       reinterpret_cast<cufftDoubleComplex*>(spectrum.ptr), CUFFT_FORWARD),
          "forward FFT (Born incident)");
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
                       std::size_t nz, cplx index_before, std::vector<double> mask, int device,
                       std::vector<cplx> fluctuation)
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
    s.fluct = std::move(fluctuation);
    s.has_host.assign(nz, 0);
    s.host.assign(nz, cplx(0.0, 0.0));
    if (!s.fluct.empty()) {
      if (s.fluct.size() != nz * nn)
        throw std::invalid_argument("fluctuation must have the shape of labels (nz, ny, nx)");
      for (std::size_t j = 0; j < nz; ++j)
        for (std::size_t q = 0; q < nn; ++q) {
          if (s.fluct[j * nn + q] == cplx(0.0, 0.0)) continue;
          const cplx idx = s.indices[static_cast<std::size_t>(s.labels[j * nn + q])];
          if (s.has_host[j] && s.host[j] != idx)
            throw std::invalid_argument("slab " + std::to_string(j) +
                                        ": fluctuation must lie in one host index");
          s.has_host[j] = 1;
          s.host[j] = idx;
        }
      s.delta = DeviceBuffer<dcplx>(nn);
      s.full3 = DeviceBuffer<dcplx>(3 * nn);
      s.work3 = DeviceBuffer<dcplx>(3 * nn);
      s.sources = DeviceBuffer<dcplx>(6 * nn);
      s.inc = DeviceBuffer<dcplx>(2 * nn);
      s.g_x = DeviceBuffer<dcplx>(2 * nn);
      s.g_rhs = DeviceBuffer<dcplx>(2 * nn);
      s.g_y = DeviceBuffer<dcplx>(2 * nn);
      s.g_w = DeviceBuffer<dcplx>(2 * nn);
      s.g_tmp = DeviceBuffer<dcplx>(2 * nn);
      // 2D FFT of three interleaved components (E_x, E_y, E_z): stride 3, distance 1, batch 3
      check(cufftPlanMany(&s.plan3, 2, dims, dims, 3, 1, dims, 3, 1, CUFFT_Z2Z, 3),
            "cufftPlanMany (Born)");
      s.plan3_ok = true;
    }
  } catch (...) {
    if (impl_->plan3_ok) cufftDestroy(impl_->plan3);
    if (impl_->plan_ok) cufftDestroy(impl_->plan);
    delete impl_;
    impl_ = nullptr;
    throw;
  }
}

CudaSolver::~CudaSolver() {
  if (impl_) {
    if (impl_->plan3_ok) cufftDestroy(impl_->plan3);
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
  // Born sources of δn as full planes on the device (F3): R1 of the forward pass, R3 backward
  std::map<std::size_t, DeviceBuffer<dcplx>> born_fwd, born_bwd;
  auto zero_plane = [&]() {
    DeviceBuffer<dcplx> b(2 * nn);
    check(cudaMemset(b.ptr, 0, 2 * nn * sizeof(dcplx)), "memset");
    return b;
  };
  auto add_plane = [&](DeviceBuffer<dcplx>& field, const DeviceBuffer<dcplx>& plane) {
    k_born_add<<<blocks(nn), kThreads>>>(field.ptr, plane.ptr, nullptr, 1.0, nn);
    check(cudaGetLastError(), "k_born_add (plane)");
  };
  for (int it = 0; it < std::max(iterations, 1); ++it) {
    src_bwd.clear();
    born_bwd.clear();
    check(cudaMemcpy(e.ptr, input.ptr, 2 * nn * sizeof(dcplx), cudaMemcpyDeviceToDevice), "copy");
    if (auto f = src_fwd.find(0); f != src_fwd.end()) add(e, f->second, 0);
    if (keep_set.count(0)) res.forward[0] = to_host(e);
    for (std::size_t p = 0; p < s.nz; ++p) {
      if (iterations > 0 && !interface_host[p].empty()) {
        s.apply(1, e, p, p + 1, next);
        src_bwd[p] = gather(next, p);
      }
      s.apply(0, e, p, p + 1, next);  // leaves F{e} in s.spectrum
      if (s.has_host[p]) {
        s.add_born_step(e, p, +1, next, true);  // Cayley step (§12)
        if (iterations > 0) {  // R1
          s.transform(e);
          s.born_sources(p, +1);
          born_bwd[p + 1] = zero_plane();
          s.add_born_emission(p, -1, born_bwd[p + 1], false);
        }
      }
      std::swap(e, next);
      // backward reflection at interface p was stored for plane p + 1
      if (auto f = src_fwd.find(p + 1); f != src_fwd.end()) add(e, f->second, p);
      if (auto f = born_fwd.find(p + 1); f != born_fwd.end()) add_plane(e, f->second);
      if (keep_set.count(p + 1)) res.forward[p + 1] = to_host(e);
    }
    res.transmitted = to_host(e);
    if (iterations == 0) break;
    src_fwd.clear();
    born_fwd.clear();
    check(cudaMemset(e.ptr, 0, 2 * nn * sizeof(dcplx)), "memset");
    if (auto f = born_bwd.find(s.nz); f != born_bwd.end()) add_plane(e, f->second);
    if (keep_set.count(s.nz)) res.backward[s.nz] = to_host(e);
    for (std::size_t p = s.nz; p > 0; --p) {
      if (s.has_host[p - 1]) {
        s.transform(e);
        s.born_sources(p - 1, -1);
        born_fwd[p] = zero_plane();
        s.add_born_emission(p - 1, +1, born_fwd[p], false);  // R3
        s.add_born_step(e, p - 1, -1, e, false);             // R2 (Cayley, §12), at plane p
      }
      if (!interface_host[p - 1].empty()) {
        s.apply(3, e, p, p - 1, next);
        src_fwd[p] = gather(next, p - 1);
      }
      s.apply(2, e, p, p - 1, next);
      std::swap(e, next);
      if (auto f = src_bwd.find(p - 1); f != src_bwd.end()) add(e, f->second, p - 1);
      if (auto f = born_bwd.find(p - 1); f != born_bwd.end()) add_plane(e, f->second);
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
