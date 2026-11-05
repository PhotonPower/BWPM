"""K1 (B0): interface matrices of wende2023 against an independent vector Fresnel operator.

Reference: the incident plane wave is completed to a 3-vector by transversality, split into s and
p with the basis of Convention A (s = k_in × N/|k_in × N|, p = k̂ × s for every wave), multiplied
by Byrnes' coefficients written with angles (eq. (6) of arXiv:1603.02720v5) and recombined with
the basis vectors of the transmitted and reflected waves. Nothing here uses bwpm.interface.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm._types import ComplexArray
from bwpm.interface import (
    fresnel_coefficients,
    longitudinal_wavenumber,
    reflection_matrix,
    transmission_matrix,
)

K0 = 2 * np.pi / 0.5e-3  # λ = 0.5 µm, k0 in 1/mm
SILVER = 0.213 + 6.697j  # wende2023, section 2.1.4


def _kz(n: complex, kp2: float) -> complex:
    q = complex(np.sqrt(complex(n * n * K0 * K0 - kp2)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _reference(
    e_perp: ComplexArray, kx: float, ky: float, n1: complex, n2: complex
) -> tuple[ComplexArray, ComplexArray]:
    """Transverse (E_x, E_y) of the transmitted and reflected wave (vector Fresnel operator)."""
    kp2 = kx * kx + ky * ky
    k1z, k2z = _kz(n1, kp2), _kz(n2, kp2)
    k_in = np.array([kx, ky, k1z])
    k_t = np.array([kx, ky, k2z])
    k_r = np.array([kx, ky, -k1z])
    e_in = np.array([e_perp[0], e_perp[1], -(kx * e_perp[0] + ky * e_perp[1]) / k1z])
    s = np.cross(k_in, [0.0, 0.0, 1.0])
    s = s / np.sqrt(s @ s)
    p_in = np.cross(k_in / (n1 * K0), s)
    p_t = np.cross(k_t / (n2 * K0), s)
    p_r = np.cross(k_r / (n1 * K0), s)
    cos1, cos2 = k1z / (n1 * K0), k2z / (n2 * K0)
    r_s = (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    r_p = (n2 * cos1 - n1 * cos2) / (n2 * cos1 + n1 * cos2)
    t_s = 2 * n1 * cos1 / (n1 * cos1 + n2 * cos2)
    t_p = 2 * n1 * cos1 / (n2 * cos1 + n1 * cos2)
    a_s, a_p = e_in @ s, e_in @ p_in
    e_t = t_s * a_s * s + t_p * a_p * p_t
    e_r = r_s * a_s * s + r_p * a_p * p_r
    return np.asarray(e_t[:2]), np.asarray(e_r[:2])


def _directions(n1: complex, count: int, seed: int) -> list[tuple[float, float]]:
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0.0, np.deg2rad(85.0), count)
    phi = rng.uniform(0.0, 2 * np.pi, count)
    k = float(np.real(n1)) * K0 * np.sin(theta)
    return [(float(a), float(b)) for a, b in zip(k * np.cos(phi), k * np.sin(phi), strict=True)]


@pytest.mark.parametrize(
    ("n1", "n2"),
    [(1.0, 1.5), (1.5, 1.0), (1.0, SILVER), (1.505, SILVER), (1.5, 1.511), (1.0, 1.0)],
    ids=["air-glass", "glass-air(TIR)", "air-silver", "IP-S-silver", "index-matched", "same"],
)
def test_matrices_equal_vector_fresnel(n1: complex, n2: complex) -> None:
    rng = np.random.default_rng(7)
    worst_t, worst_r = 0.0, 0.0
    for kx, ky in _directions(n1, 300, seed=11):
        e = rng.normal(size=2) + 1j * rng.normal(size=2)
        ref_t, ref_r = _reference(e, kx, ky, n1, n2)
        t = transmission_matrix(np.array(kx), np.array(ky), n1, n2, K0)
        r = reflection_matrix(np.array(kx), np.array(ky), n1, n2, K0)
        # errors relative to the incident field (the reflection can be tiny, e.g. index-matched)
        scale = float(np.max(np.abs(e)))
        worst_t = max(worst_t, float(np.max(np.abs(t @ e - ref_t))) / scale)
        worst_r = max(worst_r, float(np.max(np.abs(r @ e - ref_r))) / scale)
    assert worst_t < 1e-13
    assert worst_r < 1e-13


@pytest.mark.parametrize(("n1", "n2"), [(1.0, 1.5), (1.5, 1.0), (1.0, SILVER)])
def test_normal_incidence_limit(n1: complex, n2: complex) -> None:
    # exact: r = (n1 - n2)/(n1 + n2), t = 2 n1/(n1 + n2) for both components; continuous at 0
    r_exact, t_exact = (n1 - n2) / (n1 + n2), 2 * n1 / (n1 + n2)
    zero = np.zeros(1)
    assert np.allclose(transmission_matrix(zero, zero, n1, n2, K0)[0], t_exact * np.eye(2))
    assert np.allclose(reflection_matrix(zero, zero, n1, n2, K0)[0], r_exact * np.eye(2))
    tiny = np.array([1e-6 * K0])
    t_near = transmission_matrix(tiny, 0.5 * tiny, n1, n2, K0)[0]
    r_near = reflection_matrix(tiny, 0.5 * tiny, n1, n2, K0)[0]
    assert np.max(np.abs(t_near - t_exact * np.eye(2))) < 1e-9
    assert np.max(np.abs(r_near - r_exact * np.eye(2))) < 1e-9


def test_vectorised_grid_and_power_balance() -> None:
    # Poynting flux S_z ∝ Re(Σ|E_j|² k_z*)/... for transverse waves: S_z = Re(k_z*)·|E|²/(2ωμ),
    # |E|² including E_z. Lossless media, propagating transmission: R + T = 1 for every sample.
    n1, n2 = 1.0, 1.5
    k = np.linspace(-0.9, 0.9, 41) * K0 / np.sqrt(2)
    kx, ky = np.meshgrid(k, k)
    t = transmission_matrix(kx, ky, n1, n2, K0)
    r = reflection_matrix(kx, ky, n1, n2, K0)
    rng = np.random.default_rng(3)
    e = rng.normal(size=(*kx.shape, 2)) + 1j * rng.normal(size=(*kx.shape, 2))
    kz1 = longitudinal_wavenumber(n1, K0, kx * kx + ky * ky)
    kz2 = longitudinal_wavenumber(n2, K0, kx * kx + ky * ky)

    def flux(field: ComplexArray, kz: ComplexArray, sign: float) -> ComplexArray:
        ez = -sign * (kx * field[..., 0] + ky * field[..., 1]) / kz
        return np.asarray((np.sum(np.abs(field) ** 2, -1) + np.abs(ez) ** 2) * kz.real)

    e_t = np.einsum("...ij,...j->...i", t, e)
    e_r = np.einsum("...ij,...j->...i", r, e)
    balance = (flux(e_r, kz1, -1.0) + flux(e_t, kz2, 1.0)) / flux(e, kz1, 1.0)
    assert np.max(np.abs(balance - 1.0)) < 1e-13


def test_coefficients_match_byrnes_at_brewster() -> None:
    # Brewster angle air -> glass: r_p = 0 (Convention A and every other convention)
    theta_b = np.arctan(1.5)
    c = fresnel_coefficients(np.array([K0 * np.sin(theta_b)]), np.zeros(1), 1.0, 1.5, K0)
    assert abs(c.r_p[0]) < 1e-15
    assert (
        abs(
            c.r_s[0]
            - (np.cos(theta_b) - 1.5 * np.cos(np.arcsin(np.sin(theta_b) / 1.5)))
            / (np.cos(theta_b) + 1.5 * np.cos(np.arcsin(np.sin(theta_b) / 1.5)))
        )
        < 1e-15
    )
