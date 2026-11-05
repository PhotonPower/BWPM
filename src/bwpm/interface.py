"""Transmission and reflection matrices of a plane interface for the transverse field (E_x, E_y).

Specification: docs/spezifikation.md, section 3; source wende2023, eqs. (2), (3), (8), (9).

A plane wave with transverse frequency k⊥ = (k_x, k_y) in medium r (index ñ_r) hits the plane
z = const with medium s (index ñ_s) behind it. The matrices map the transverse components
(E_x, E_y) of the incident wave to those of the transmitted (+z) and reflected (−z) wave, in the
global frame. The propagation factor 𝒫 of the source is not included here.

Conventions (docs/konventionen.md): fields ∝ exp(i(k·r − ωt)), k_z with Im k_z ≥ 0, Fresnel
Convention A (Byrnes, arXiv:1603.02720v5, eq. (6)); TE = s, TM = p. With these coefficients the
matrices of the source equal the vector Fresnel operator (tests/test_interface.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from bwpm._types import ComplexArray, FloatArray


def longitudinal_wavenumber(index: complex, k0: float, k_perp_sq: FloatArray) -> ComplexArray:
    """k_z = sqrt(ñ²k₀² − |k⊥|²) with Im k_z ≥ 0, and Re k_z ≥ 0 if Im k_z = 0 (1/mm).

    Args:
        index: complex refractive index ñ = n + iκ, κ ≥ 0.
        k0: vacuum wavenumber 2π/λ in 1/mm.
        k_perp_sq: |k⊥|² in 1/mm².

    |k_z| < 10⁻¹²·k₀ (exactly grazing) is replaced by i·10⁻¹²·k₀ to avoid 0/0 in Fresnel.
    """
    kz = np.sqrt(np.asarray(complex(index) ** 2 * k0 * k0 - k_perp_sq, dtype=np.complex128))
    flip = (kz.imag < 0) | ((kz.imag == 0) & (kz.real < 0))
    kz = np.where(flip, -kz, kz)
    # exactly grazing waves (k_z = 0, e.g. |k⊥| = k on the FFT grid) make Fresnel 0/0; they carry
    # no power. Shift them to a tiny evanescent value (as FourierGrating does for β).
    tiny = 1e-12 * k0
    return np.asarray(np.where(np.abs(kz) < tiny, 1j * tiny, kz))


@dataclass(frozen=True)
class FresnelCoefficients:
    """Amplitude coefficients in Convention A for every k⊥ (arrays of the shape of k_x).

    Attributes:
        r_s, r_p, t_s, t_p: reflection and transmission coefficients (TE = s, TM = p).
        kz_in, kz_out: k_z in the medium of incidence and of transmission (1/mm).
    """

    r_s: ComplexArray
    r_p: ComplexArray
    t_s: ComplexArray
    t_p: ComplexArray
    kz_in: ComplexArray
    kz_out: ComplexArray


def fresnel_coefficients(
    k_x: FloatArray, k_y: FloatArray, index_in: complex, index_out: complex, k0: float
) -> FresnelCoefficients:
    """Fresnel coefficients of Convention A at the transverse frequencies (k_x, k_y) (1/mm).

    With q_j = k_z,j/k₀ = ñ_j cos θ_j: r_s = (q₁ − q₂)/(q₁ + q₂), t_s = 2q₁/(q₁ + q₂),
    r_p = (ñ₂² q₁ − ñ₁² q₂)/(ñ₂² q₁ + ñ₁² q₂), t_p = 2ñ₁ñ₂ q₁/(ñ₂² q₁ + ñ₁² q₂)
    (Byrnes eq. (6) multiplied by ñ₁ñ₂). Normal incidence: r_p = −r_s.
    """
    k_perp_sq = np.asarray(k_x, dtype=np.float64) ** 2 + np.asarray(k_y, dtype=np.float64) ** 2
    n1, n2 = complex(index_in), complex(index_out)
    kz1 = longitudinal_wavenumber(n1, k0, k_perp_sq)
    kz2 = longitudinal_wavenumber(n2, k0, k_perp_sq)
    q1, q2 = kz1 / k0, kz2 / k0
    denominator_p = n2 * n2 * q1 + n1 * n1 * q2
    return FresnelCoefficients(
        r_s=np.asarray((q1 - q2) / (q1 + q2)),
        r_p=np.asarray((n2 * n2 * q1 - n1 * n1 * q2) / denominator_p),
        t_s=np.asarray(2.0 * q1 / (q1 + q2)),
        t_p=np.asarray(2.0 * n1 * n2 * q1 / denominator_p),
        kz_in=kz1,
        kz_out=kz2,
    )


def transmission_matrix(
    k_x: FloatArray, k_y: FloatArray, index_in: complex, index_out: complex, k0: float
) -> ComplexArray:
    """Matrix T of wende2023, eqs. (2) and (3), without propagation factor 𝒫, shape (..., 2, 2).

    T = 1/(k⊥² k_r k_s) [[t_TM k_x² Γ + t_TE k_y² k_r k_s, (t_TM Γ − t_TE k_r k_s) k_x k_y],
                         [(t_TM Γ − t_TE k_r k_s) k_x k_y, t_TM k_y² Γ + t_TE k_x² k_r k_s]],
    Γ = k_z,r k_z,s + k_z,s k⊥²/k_z,r, k_r = ñ_r k₀, k_s = ñ_s k₀. At k⊥ = 0 the limit t·I is
    used (t_TE = t_TM there).

    Args:
        k_x, k_y: transverse frequency in 1/mm (broadcastable arrays).
        index_in, index_out: ñ_r of the incident and ñ_s of the transmitting medium.
        k0: vacuum wavenumber in 1/mm.
    """
    c = fresnel_coefficients(k_x, k_y, index_in, index_out, k0)
    kx = np.asarray(k_x, dtype=np.float64)
    ky = np.asarray(k_y, dtype=np.float64)
    kp2 = kx * kx + ky * ky
    kr, ks = complex(index_in) * k0, complex(index_out) * k0
    zero = kp2 == 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = c.kz_in * c.kz_out + c.kz_out * kp2 / c.kz_in
        scale = 1.0 / (kp2 * kr * ks)
        t00 = (c.t_p * kx * kx * gamma + c.t_s * ky * ky * kr * ks) * scale
        t01 = (c.t_p * gamma - c.t_s * kr * ks) * kx * ky * scale
        t11 = (c.t_p * ky * ky * gamma + c.t_s * kx * kx * kr * ks) * scale
    t00 = np.where(zero, c.t_s, t00)
    t01 = np.where(zero, 0.0, t01)
    t11 = np.where(zero, c.t_s, t11)
    return np.asarray(np.stack([np.stack([t00, t01], -1), np.stack([t01, t11], -1)], -2))


def reflection_matrix(
    k_x: FloatArray, k_y: FloatArray, index_in: complex, index_out: complex, k0: float
) -> ComplexArray:
    """Matrix R of wende2023, eqs. (8) and (9), without propagation factor 𝒫, shape (..., 2, 2).

    R = 1/(k⊥² k_r²) [[−r_TM k_x² Γ_r + r_TE k_y² k_r², −(r_TM Γ_r + r_TE k_r²) k_x k_y],
                      [−(r_TM Γ_r + r_TE k_r²) k_x k_y, −r_TM k_y² Γ_r + r_TE k_x² k_r²]],
    Γ_r = k_z,r² + k⊥² = k_r². The reflected wave runs in −z; its transverse components are given
    in the global frame. At k⊥ = 0 the limit r_s·I is used (Convention A: r_p = −r_s there).
    """
    c = fresnel_coefficients(k_x, k_y, index_in, index_out, k0)
    kx = np.asarray(k_x, dtype=np.float64)
    ky = np.asarray(k_y, dtype=np.float64)
    kp2 = kx * kx + ky * ky
    kr = complex(index_in) * k0
    zero = kp2 == 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = c.kz_in * c.kz_in + kp2
        scale = 1.0 / (kp2 * kr * kr)
        r00 = (-c.r_p * kx * kx * gamma + c.r_s * ky * ky * kr * kr) * scale
        r01 = -(c.r_p * gamma + c.r_s * kr * kr) * kx * ky * scale
        r11 = (-c.r_p * ky * ky * gamma + c.r_s * kx * kx * kr * kr) * scale
    r00 = np.where(zero, c.r_s, r00)
    r01 = np.where(zero, 0.0, r01)
    r11 = np.where(zero, c.r_s, r11)
    return np.asarray(np.stack([np.stack([r00, r01], -1), np.stack([r01, r11], -1)], -2))
