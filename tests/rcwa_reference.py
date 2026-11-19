"""Rigorous coupled-wave analysis (RCWA) of one lamellar layer: independent reference for K17.

Test-only code. A layer of thickness L, periodic in y with period Λ, index ñ₂ for
y ∈ [y₀, y₁) and ñ₁ elsewhere, lies between two half-spaces of index ñ₀. A plane wave of unit
amplitude falls in along +z (k_y = 0). Fields ∝ exp(i(k·r − ωt)).

TE (E_x): S″ = (K_y² − k₀²[[ε]]) S. TM (H_x): U″ = [[1/ε]]⁻¹(K_y [[ε]]⁻¹ K_y − k₀²) U, the
factorisation rules of li1996 for the lamellar layer (Laurent for 1/ε acting on ∂_z H, inverse
rule for ε acting on ∂_y H). Continuity of E_x and ∂_z E_x (TE) or of H_x and ε⁻¹∂_z H_x (TM) at
both faces gives one linear system (moharam1995, formulation without the transmittance matrix;
the layer modes are referred to the face they decay away from, so no exponential overflows).
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from bwpm._types import ComplexArray


def _toeplitz(
    inside: complex, outside: complex, y0: float, y1: float, period: float, n: int
) -> ComplexArray:
    """Toeplitz matrix [[f]]_{pq} = f̂_{p−q} of the piecewise constant f, orders −n…n."""
    kg = 2.0 * np.pi / period
    q = np.arange(-2 * n, 2 * n + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        coeff = (
            (inside - outside)
            * (np.exp(-1j * q * kg * y0) - np.exp(-1j * q * kg * y1))
            / (1j * q * kg * period)
        )
    coeff = np.where(q == 0, outside + (inside - outside) * (y1 - y0) / period, coeff)
    orders = np.arange(-n, n + 1)
    return np.asarray(coeff[orders[:, None] - orders[None, :] + 2 * n], dtype=np.complex128)


def lamellar_transmission(
    index_background: complex,
    index_outside: complex,
    index_inside: complex,
    period: float,
    y0: float,
    y1: float,
    thickness: float,
    wavelength: float,
    n_orders: int,
    polarization: Literal["TE", "TM"],
) -> dict[int, complex]:
    """Transmitted transverse field E_x (TE) or E_y (TM) per diffraction order at z = L.

    Lengths in one unit (e.g. mm); order m has k_y = 2πm/Λ. The amplitude refers to the plane wave
    E_m e^(i(k_y y + k_z (z − L))) in the half-space behind the layer, the incident wave has
    amplitude 1 at z = 0.
    """
    k0 = 2.0 * np.pi / wavelength
    k = complex(index_background) * k0
    orders = np.arange(-n_orders, n_orders + 1)
    ky = orders * 2.0 * np.pi / period
    kz = np.sqrt((k * k - ky * ky).astype(np.complex128))
    kz = np.where(kz.imag < 0, -kz, kz)
    size = orders.size
    unit = np.eye(size)
    ky_mat = np.diag(ky).astype(np.complex128)
    eps = _toeplitz(index_inside**2, index_outside**2, y0, y1, period, n_orders)
    operator: ComplexArray
    coupling: ComplexArray
    if polarization == "TE":
        operator = k0 * k0 * eps - ky_mat @ ky_mat
        coupling = unit.astype(np.complex128)
        h_factor = 1.0 + 0j
    else:
        inv_eps = _toeplitz(index_inside**-2, index_outside**-2, y0, y1, period, n_orders)
        operator = -np.linalg.solve(inv_eps, ky_mat @ np.linalg.solve(eps, ky_mat) - k0 * k0 * unit)
        coupling = inv_eps
        h_factor = complex(index_background) ** -2
    gamma_sq, modes = np.linalg.eig(operator)
    gamma = np.sqrt(gamma_sq.astype(np.complex128))
    gamma = np.where(gamma.imag < 0, -gamma, gamma)
    decay = np.diag(np.exp(1j * gamma * thickness))
    deriv = coupling @ modes @ np.diag(gamma)
    zero = np.zeros((size, size))
    kz_mat = np.diag(kz)
    # unknowns: reflected r, layer modes c⁺ (from z = 0) and c⁻ (from z = L), transmitted t
    system = np.block(
        [
            [unit, -modes, -modes @ decay, zero],
            [h_factor * kz_mat, deriv, -deriv @ decay, zero],
            [zero, modes @ decay, modes, -unit],
            [zero, deriv @ decay, -deriv, -h_factor * kz_mat],
        ]
    )
    incident = np.zeros(size, dtype=np.complex128)
    incident[n_orders] = 1.0
    rhs = np.concatenate([-incident, h_factor * kz * incident, np.zeros(2 * size)])
    t = np.linalg.solve(system, rhs)[3 * size :]
    if polarization == "TM":
        t = t * kz / k  # H_x → E_y of a +z wave: E_y = −k_z H_x/(ω ε₀ ñ₀²), same factor for all m
    return {int(m): complex(v) for m, v in zip(orders, t, strict=True)}
