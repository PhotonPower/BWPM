// Python module bwpm._core (pybind11): vectorised wrappers of the C++ core.
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <stdexcept>
#include <vector>

#include "bwpm/interface.hpp"

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
  m.def("longitudinal_wavenumber", &bwpm::longitudinal_wavenumber, py::arg("index"),
        py::arg("k0"), py::arg("k_perp_sq"));
}
