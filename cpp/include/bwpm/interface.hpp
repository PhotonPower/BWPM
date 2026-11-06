// Interface matrices of the BWPM for the transverse field (E_x, E_y).
//
// Specification: docs/spezifikation.md, section 3; source wende2023, eqs. (2), (3), (8), (9).
// Conventions (docs/konventionen.md): fields ∝ exp(i(k·r − ωt)), k_z with Im k_z ≥ 0, Fresnel
// Convention A (Byrnes eq. (6)), TE = s, TM = p. Units: k in 1/mm. Mirrors bwpm/interface.py.
#pragma once

#include <array>
#include <complex>

namespace bwpm {

using cplx = std::complex<double>;
using Matrix2 = std::array<cplx, 4>;  // row-major [[m00, m01], [m10, m11]]

// k_z = sqrt(ñ²k₀² − k⊥²), Im k_z ≥ 0 (Re k_z ≥ 0 if Im k_z = 0); |k_z| < 1e-12·k₀ -> i·1e-12·k₀.
cplx longitudinal_wavenumber(cplx index, double k0, double k_perp_sq);

struct Fresnel {
  cplx r_s, r_p, t_s, t_p, kz_in, kz_out;
};

// Convention A with q_j = k_z,j/k₀.
Fresnel fresnel_coefficients(double kx, double ky, cplx index_in, cplx index_out, double k0);

// T of eqs. (2), (3) without propagation factor; k⊥ = 0: t·I.
Matrix2 transmission_matrix(double kx, double ky, cplx index_in, cplx index_out, double k0);

// R of eqs. (8), (9) without propagation factor; k⊥ = 0: r_s·I.
Matrix2 reflection_matrix(double kx, double ky, cplx index_in, cplx index_out, double k0);

}  // namespace bwpm
