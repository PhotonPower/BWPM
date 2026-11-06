// Python module bwpm._core (pybind11): vectorised wrappers of the C++ core.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <memory>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include "bwpm/interface.hpp"
#include "bwpm/solver.hpp"
#ifdef BWPM_HAVE_CUDA
#include "bwpm/solver_cuda.hpp"
#endif

namespace py = pybind11;
using bwpm::cplx;

namespace {

using F64 = py::array_t<double, py::array::c_style | py::array::forcecast>;

// Apply a 2×2-matrix function to every (k_x, k_y); result shape (..., 2, 2), complex128.
template <class Fn>
py::array_t<cplx> matrices(const F64& kx, const F64& ky, cplx n_in, cplx n_out, double k0, Fn fn) {
  if (kx.ndim() != ky.ndim()) throw std::invalid_argument("k_x and k_y must have the same shape");
  std::vector<py::ssize_t> shape(kx.shape(), kx.shape() + kx.ndim());
  for (py::ssize_t d = 0; d < kx.ndim(); ++d)
    if (kx.shape(d) != ky.shape(d)) throw std::invalid_argument("k_x and k_y must have the same shape");
  shape.push_back(2);
  shape.push_back(2);
  py::array_t<cplx> out(shape);
  const double* px = kx.data();
  const double* py_ = ky.data();
  cplx* po = out.mutable_data();
  const py::ssize_t n = kx.size();
  {
    py::gil_scoped_release release;
    for (py::ssize_t i = 0; i < n; ++i) {
      const bwpm::Matrix2 m = fn(px[i], py_[i], n_in, n_out, k0);
      for (int j = 0; j < 4; ++j) po[4 * i + j] = m[j];
    }
  }
  return out;
}

using C128 = py::array_t<cplx, py::array::c_style | py::array::forcecast>;
using I32 = py::array_t<std::int32_t, py::array::c_style | py::array::forcecast>;

py::array_t<cplx> as_array(const std::vector<cplx>& v, std::vector<py::ssize_t> shape) {
  py::array_t<cplx> out(shape);
  std::copy(v.begin(), v.end(), out.mutable_data());
  return out;
}

template <class S>
py::dict run_solver(S& self, const C128& field, int iterations, std::vector<std::size_t> keep) {
  if (field.ndim() != 3 || field.shape(2) != 2)
    throw std::invalid_argument("field must have shape (ny, nx, 2)");
  const std::vector<py::ssize_t> shape{field.shape(0), field.shape(1), 2};
  std::vector<cplx> in(field.data(), field.data() + field.size());
  bwpm::RunResult r;
  {
    py::gil_scoped_release release;
    r = self.run(in, iterations, keep);
  }
  py::dict out;
  out["transmitted"] = as_array(r.transmitted, shape);
  out["reflected"] = as_array(r.reflected, shape);
  py::list ti, ri;
  for (const auto& v : r.transmitted_iter) ti.append(as_array(v, shape));
  for (const auto& v : r.reflected_iter) ri.append(as_array(v, shape));
  out["transmitted_iter"] = ti;
  out["reflected_iter"] = ri;
  py::dict fw, bw;
  for (const auto& [p, v] : r.forward) fw[py::int_(p)] = as_array(v, shape);
  for (const auto& [p, v] : r.backward) bw[py::int_(p)] = as_array(v, shape);
  out["forward"] = fw;
  out["backward"] = bw;
  out["stored_values"] = r.stored_values;
  return out;
}

struct SolverInput {
  bwpm::Grid grid;
  std::vector<cplx> indices;
  std::vector<std::int32_t> labels;
  std::size_t nz;
  std::vector<double> mask;
};

SolverInput solver_input(std::size_t ny, std::size_t nx, double dy, double dx, double dz,
                         double k0, const C128& indices, const I32& labels,
                         const std::optional<py::array_t<double, py::array::c_style | py::array::forcecast>>& mask) {
  if (labels.ndim() != 3) throw std::invalid_argument("labels must have shape (nz, ny, nx)");
  SolverInput in{{ny, nx, dy, dx, dz, k0},
                 std::vector<cplx>(indices.data(), indices.data() + indices.size()),
                 std::vector<std::int32_t>(labels.data(), labels.data() + labels.size()),
                 static_cast<std::size_t>(labels.shape(0)),
                 {}};
  if (mask) in.mask.assign(mask->data(), mask->data() + mask->size());
  return in;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "C++ core of the BWPM (stage B4). Same conventions as the Python reference.";
  m.def(
      "transmission_matrix",
      [](const F64& kx, const F64& ky, cplx n_in, cplx n_out, double k0) {
        return matrices(kx, ky, n_in, n_out, k0, bwpm::transmission_matrix);
      },
      py::arg("k_x"), py::arg("k_y"), py::arg("index_in"), py::arg("index_out"), py::arg("k0"),
      "T of wende2023 eqs. (2), (3) without propagation factor, shape (..., 2, 2).");
  m.def(
      "reflection_matrix",
      [](const F64& kx, const F64& ky, cplx n_in, cplx n_out, double k0) {
        return matrices(kx, ky, n_in, n_out, k0, bwpm::reflection_matrix);
      },
      py::arg("k_x"), py::arg("k_y"), py::arg("index_in"), py::arg("index_out"), py::arg("k0"),
      "R of wende2023 eqs. (8), (9) without propagation factor, shape (..., 2, 2).");
  using F64Opt = std::optional<py::array_t<double, py::array::c_style | py::array::forcecast>>;
  py::class_<bwpm::Solver>(m, "Solver",
                           "FPWPM/BWPM on the CPU in 2D (nx = 1) and 3D (stage B4c).")
      .def(py::init([](std::size_t ny, std::size_t nx, double dy, double dx, double dz, double k0,
                       const C128& indices, const I32& labels, cplx index_before,
                       const F64Opt& mask, int threads) {
             SolverInput in = solver_input(ny, nx, dy, dx, dz, k0, indices, labels, mask);
             return bwpm::Solver(in.grid, std::move(in.indices), std::move(in.labels), in.nz,
                                 index_before, std::move(in.mask), threads);
           }),
           py::arg("ny"), py::arg("nx"), py::arg("dy"), py::arg("dx"), py::arg("dz"), py::arg("k0"),
           py::arg("indices"), py::arg("labels"), py::arg("index_before"), py::arg("mask"),
           py::arg("threads") = 1)
      .def("run", &run_solver<bwpm::Solver>, py::arg("field"), py::arg("iterations"),
           py::arg("keep") = std::vector<std::size_t>{});

#ifdef BWPM_HAVE_CUDA
  m.attr("HAS_CUDA") = true;
  m.def("cuda_device_count", &bwpm::cuda_device_count);
  py::class_<bwpm::CudaSolver>(m, "CudaSolver", "FPWPM/BWPM on the GPU (stage B4d).")
      .def(py::init([](std::size_t ny, std::size_t nx, double dy, double dx, double dz, double k0,
                       const C128& indices, const I32& labels, cplx index_before,
                       const F64Opt& mask, int device) {
             SolverInput in = solver_input(ny, nx, dy, dx, dz, k0, indices, labels, mask);
             return std::make_unique<bwpm::CudaSolver>(in.grid, std::move(in.indices),
                                                       std::move(in.labels), in.nz, index_before,
                                                       std::move(in.mask), device);
           }),
           py::arg("ny"), py::arg("nx"), py::arg("dy"), py::arg("dx"), py::arg("dz"), py::arg("k0"),
           py::arg("indices"), py::arg("labels"), py::arg("index_before"), py::arg("mask"),
           py::arg("device") = 0)
      .def("run", &run_solver<bwpm::CudaSolver>, py::arg("field"), py::arg("iterations"),
           py::arg("keep") = std::vector<std::size_t>{});
#else
  m.attr("HAS_CUDA") = false;
  m.def("cuda_device_count", []() { return 0; });
#endif
  m.def("longitudinal_wavenumber", &bwpm::longitudinal_wavenumber, py::arg("index"),
        py::arg("k0"), py::arg("k_perp_sq"));
}
