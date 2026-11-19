"""B5a: analysis helpers for the 3D fibre-facet pilot (test case K6 in 3D).

Used by `experiments/b3_scaling.py`; they turn the planes returned by `Bidirectional3D` into
the quantities of K6: flux, reflectance against an independent Fresnel reference, and the axial
position of a focus. Conventions as everywhere: lengths in mm, k in 1/mm, fields ∝ exp(i(k·r − ωt)),
transverse field (E_x, E_y) of shape (ny, nx, 2) on the centred grid of `Grid3D`
(x_i = (i − nx//2)·dx). Specification: `docs/spezifikation.md`; pilot: `docs/pilot-imc-2d.md`.
"""

from __future__ import annotations

import numpy as np

from bwpm._types import ComplexArray, FloatArray
from bwpm.volume3d import Grid3D

_GRAZING = 1e-9  # components with cos θ below this are treated as non-propagating


def _spectrum(field: ComplexArray, grid: Grid3D) -> tuple[ComplexArray, FloatArray, FloatArray]:
    """Angular spectrum of both components (FFT, no shift) and the (k_x, k_y) grids."""
    kx, ky = grid.k_grid()
    a = np.fft.fft2(field, axes=(0, 1))
    return np.asarray(a), kx, ky


def _polarisation_split(
    a: ComplexArray, kx: FloatArray, ky: FloatArray, cos_t: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """|E_s|² and |E_p|² per mode from the transverse spectrum a (ny, nx, 2).

    s along ẑ × k̂ = (−k_y, k_x, 0)/k_⊥, p in the plane of k and ẑ. With E ⟂ k:
    E_s = E_t · ŝ and E_p = (E_t · k̂_⊥)/cos θ, so |E|² = |E_s|² + |E_p|².
    """
    kp = np.hypot(kx, ky)
    safe = np.where(kp > 0, kp, 1.0)
    # at k_⊥ = 0 s and p are degenerate; take ŝ = ŷ, k̂_⊥ = x̂ (the power sum is unchanged)
    sx = np.where(kp > 0, -ky / safe, 0.0)
    sy = np.where(kp > 0, kx / safe, 1.0)
    px = np.where(kp > 0, kx / safe, 1.0)
    py = np.where(kp > 0, ky / safe, 0.0)
    e_s = a[..., 0] * sx + a[..., 1] * sy
    e_p = (a[..., 0] * px + a[..., 1] * py) / np.where(cos_t > _GRAZING, cos_t, 1.0)
    return np.abs(e_s) ** 2, np.abs(e_p) ** 2


def _cos_theta(grid: Grid3D, index: float) -> tuple[FloatArray, FloatArray, FloatArray]:
    kx, ky = grid.k_grid()
    k = index * grid.k0
    cos_sq = 1.0 - (kx**2 + ky**2) / k**2
    cos_t = np.sqrt(np.maximum(cos_sq, 0.0))
    return kx, ky, np.where(cos_t > _GRAZING, cos_t, 0.0)


def plane_power(field: ComplexArray, grid: Grid3D, index: float) -> float:
    """Flux through the plane in arbitrary but fixed units (the same for every call).

    P = n Σ cos θ (|E_s|² + |E_p|²) over the propagating modes (Parseval, no ΔA·Z₀ factors):
    a plane wave carries n cos θ |E|². Evanescent modes carry no flux. Direction-independent, so
    forward and backward fields of one medium (or fields in different media, via n) compare.
    """
    a, _, _ = _spectrum(field, grid)
    kx, ky, cos_t = _cos_theta(grid, index)
    s2, p2 = _polarisation_split(a, kx, ky, cos_t)
    return float(index * np.sum(cos_t * (s2 + p2)))


def fresnel_reflectances(n1: float, n2: float, cos_i: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Power reflectances (R_s, R_p) of the interface n1 → n2, lossless, cos of the incidence angle.

    Fresnel amplitude coefficients with Snell's law; total internal reflection gives 1.
    """
    cos_i = np.asarray(cos_i, dtype=np.float64)
    sin_t2 = (n1 / n2) ** 2 * (1.0 - cos_i**2)
    tir = sin_t2 >= 1.0
    cos_t = np.sqrt(np.where(tir, 0.0, 1.0 - sin_t2))
    rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return np.where(tir, 1.0, rs**2), np.where(tir, 1.0, rp**2)


def reference_reflectance(field: ComplexArray, grid: Grid3D, n1: float, n2: float) -> float:
    """Flux-weighted Fresnel reflectance of one interface n1 → n2 for the incident field.

    Independent reference for the reflected power of K6: every plane-wave mode of the incident
    field is weighted with R_s or R_p of its own angle. It ignores multiple reflections, the
    thickness of the slab and the transmitted-side interfaces.
    """
    a, _, _ = _spectrum(field, grid)
    kx, ky, cos_t = _cos_theta(grid, n1)
    s2, p2 = _polarisation_split(a, kx, ky, cos_t)
    rs, rp = fresnel_reflectances(n1, n2, cos_t)
    flux = cos_t * (s2 + p2)
    total = float(np.sum(flux))
    if total == 0.0:
        raise ValueError("the field carries no propagating flux")
    return float(np.sum(cos_t * (s2 * rs + p2 * rp)) / total)


def axial_focus(
    field: ComplexArray,
    grid: Grid3D,
    index: float,
    offsets: FloatArray,
    na: float | None = None,
) -> tuple[float, float]:
    """Distance along the propagation direction (mm) at which the on-axis intensity peaks.

    The field at x = y = 0 after a distance d is Σ_k A_k exp(i k_z d) with A_k the angular
    spectrum referred to the grid centre; this holds for forward and backward waves alike
    (d counts in the direction of travel). Intensity |E_x|² + |E_y|² + |E_z|², E_z = −k_⊥·E_t/k_z.
    Returns (offset of the maximum, maximum); positive d = further along the direction of travel.
    Only propagating modes enter, with `na` only those with sin θ ≤ na / index (a speed-up and a
    guard against noise outside the pupil). The resolution is the step of `offsets`.
    """
    a, kx, ky = _spectrum(field, grid)
    a = a * np.exp(1j * (kx * (grid.nx // 2) * grid.dx + ky * (grid.ny // 2) * grid.dy))[..., None]
    _, _, cos_t = _cos_theta(grid, index)
    prop = cos_t > 0
    if na is not None:
        prop &= np.hypot(kx, ky) <= na * grid.k0
    kz = (index * grid.k0 * cos_t)[prop]
    ax, ay = a[..., 0][prop], a[..., 1][prop]
    az = -(kx[prop] * ax + ky[prop] * ay) / kz
    best, peak = 0.0, -1.0
    for d in np.asarray(offsets, dtype=np.float64):
        ph = np.exp(1j * kz * d)
        intensity = sum(abs(np.sum(c * ph)) ** 2 for c in (ax, ay, az)) / (grid.nx * grid.ny) ** 2
        if intensity > peak:
            best, peak = float(d), float(intensity)
    return best, peak
