// Interface matrices, see cpp/include/bwpm/interface.hpp.
#include "bwpm/interface.hpp"

#include <cmath>

namespace bwpm {

cplx longitudinal_wavenumber(cplx index, double k0, double k_perp_sq) {
  cplx kz = std::sqrt(index * index * (k0 * k0) - k_perp_sq);
  if (kz.imag() < 0.0 || (kz.imag() == 0.0 && kz.real() < 0.0)) kz = -kz;
  const double tiny = 1e-12 * k0;  // exactly grazing: avoid 0/0 in Fresnel (as the Python code)
  if (std::abs(kz) < tiny) kz = cplx(0.0, tiny);
  return kz;
}

Fresnel fresnel_coefficients(double kx, double ky, cplx n1, cplx n2, double k0) {
  const double kp2 = kx * kx + ky * ky;
  const cplx kz1 = longitudinal_wavenumber(n1, k0, kp2);
  const cplx kz2 = longitudinal_wavenumber(n2, k0, kp2);
  const cplx q1 = kz1 / k0, q2 = kz2 / k0;
  const cplx den_p = n2 * n2 * q1 + n1 * n1 * q2;
  return {(q1 - q2) / (q1 + q2), (n2 * n2 * q1 - n1 * n1 * q2) / den_p, 2.0 * q1 / (q1 + q2),
          2.0 * n1 * n2 * q1 / den_p, kz1, kz2};
}

Matrix2 transmission_matrix(double kx, double ky, cplx n1, cplx n2, double k0) {
  const Fresnel c = fresnel_coefficients(kx, ky, n1, n2, k0);
  const double kp2 = kx * kx + ky * ky;
  if (kp2 == 0.0) return {c.t_s, 0.0, 0.0, c.t_s};
  const cplx kr = n1 * k0, ks = n2 * k0;
  const cplx gamma = c.kz_in * c.kz_out + c.kz_out * kp2 / c.kz_in;
  const cplx scale = 1.0 / (kp2 * kr * ks);
  const cplx off = (c.t_p * gamma - c.t_s * kr * ks) * kx * ky * scale;
  return {(c.t_p * kx * kx * gamma + c.t_s * ky * ky * kr * ks) * scale, off, off,
          (c.t_p * ky * ky * gamma + c.t_s * kx * kx * kr * ks) * scale};
}

Matrix2 reflection_matrix(double kx, double ky, cplx n1, cplx n2, double k0) {
  const Fresnel c = fresnel_coefficients(kx, ky, n1, n2, k0);
  const double kp2 = kx * kx + ky * ky;
  if (kp2 == 0.0) return {c.r_s, 0.0, 0.0, c.r_s};
  const cplx kr = n1 * k0;
  const cplx gamma = c.kz_in * c.kz_in + kp2;
  const cplx scale = 1.0 / (kp2 * kr * kr);
  const cplx off = -(c.r_p * gamma + c.r_s * kr * kr) * kx * ky * scale;
  return {(-c.r_p * kx * kx * gamma + c.r_s * ky * ky * kr * kr) * scale, off, off,
          (-c.r_p * ky * ky * gamma + c.r_s * kx * kx * kr * kr) * scale};
}

}  // namespace bwpm
