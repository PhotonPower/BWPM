"""First Born approximation of the forward-scattered field of a weak volume (reference for K17).

Specification: docs/spezifikation.md, section 9; sources bornwolf1999 (first Born approximation),
chew1995 (spectral form of the scalar and dyadic Green's function). Both are cited
bibliographically, so no equation numbers are given here.

Geometry as in `bwpm.fpwpm` and `bwpm.volume3d`: slab j of the labelled volume occupies
z ∈ [j·dz, (j + 1)·dz] (z = 0 at the input plane), pixel (j, n) is the cell |y − y_n| ≤ dy/2 (and
|x − x_m| ≤ dx/2 in 3D) with constant index. The volume is embedded in a homogeneous background
ñ₀ and periodic transversally. A plane wave of unit amplitude, E = e^(i k z) with k = ñ₀k₀, falls
in along +z. Fields ∝ exp(i(k·r − ωt)).

With the potential V = k₀²(ñ² − ñ₀²), the scalar scattered amplitude with transverse frequency k⊥
behind all scatterers (z ≥ z′) is

    S(k⊥, z) = i/(2k_z) · Σ_j V̂_j(k⊥) · ∫_slab j e^(i k_z (z − z′)) e^(i k z′) dz′,

k_z = sqrt(k² − |k⊥|²) with Im k_z ≥ 0. V̂_j is the exact Fourier coefficient of the pixel
staircase: the DFT of the pixel values times sinc(k_x dx/2)·sinc(k_y dy/2). The dyadic Green's
function turns the incident polarisation ê into (I − k k/k²)·ê; its transverse block is

    E⊥ = S · [[1 − k_x²/k², −k_x k_y/k²], [−k_x k_y/k², 1 − k_y²/k²]] · (e_x, e_y).

In 2D (k_x = 0) this is S for TE (E_x) and S·k_z²/k² for TM (E_y).
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from bwpm._types import ComplexArray, FloatArray, IntArray
from bwpm.fpwpm import Grid2D
from bwpm.interface import longitudinal_wavenumber
from bwpm.volume3d import Grid3D


def _scalar_amplitude(v_hat: ComplexArray, kz: ComplexArray, k: complex, dz: float) -> ComplexArray:
    """S(k⊥) at the output plane from V̂_j(k⊥), shape (nz, …) → (…)."""
    nz = v_hat.shape[0]
    z_out = nz * dz
    # ∫_slab j e^(i k_z (z_out − z′)) e^(i k z′) dz′
    #   = e^(i k z_j) · e^(i k_z (z_out − z_j − dz)) · (e^(i k dz) − e^(i k_z dz)) / (i (k − k_z)),
    # written so that every factor stays bounded for evanescent orders (Im k_z > 0)
    z_left = np.arange(nz) * dz
    shape = (nz,) + (1,) * kz.ndim
    z_left = z_left.reshape(shape)
    mismatch = k - kz
    small = np.abs(mismatch * dz) < 1e-8
    safe = np.where(small, 1.0, mismatch)
    within = np.where(
        small,
        dz * np.exp(0.5j * (k + kz) * dz),
        (np.exp(1j * k * dz) - np.exp(1j * kz * dz)) / (1j * safe),
    )
    travel = np.exp(1j * kz[None, ...] * (z_out - z_left - dz))
    slab_integral = np.exp(1j * k * z_left) * travel * within[None, ...]
    return np.asarray(1j / (2.0 * kz) * np.sum(v_hat * slab_integral, axis=0))


def _potential(
    indices: list[complex], labels: IntArray, index_background: complex, k0: float
) -> ComplexArray:
    """V = k₀²(ñ² − ñ₀²) per pixel in 1/mm²."""
    index_map = np.asarray(indices, dtype=np.complex128)[labels]
    return np.asarray(k0 * k0 * (index_map**2 - complex(index_background) ** 2))


def born_forward_2d(
    grid: Grid2D,
    indices: list[complex],
    labels: IntArray,
    index_background: complex,
    polarization: Literal["TE", "TM"],
) -> ComplexArray:
    """Plane-wave amplitudes of the first-order scattered transverse field at the output plane, 2D.

    Args:
        grid: transverse grid and slab thickness dz (mm), wavelength (µm).
        indices: complex refractive index per material label.
        labels: material label per pixel, shape (nz, ny); slab j spans z ∈ [j·dz, (j + 1)·dz].
        index_background: index ñ₀ of the background in which the incident wave travels.
        polarization: "TE" (E_x) or "TM" (E_y) of the normally incident plane wave.

    Returns:
        Amplitudes a_m = fft(E_s)[m] / ny (numpy order, k_y,m = grid.ky[m]) of the scattered
        transverse component, E_x (TE) or E_y (TM), sampled on `grid.y` at z_out = nz·dz,
        shape (ny,). The incident wave e^(i k z) itself is not included.
    """
    labels = np.asarray(labels)
    if labels.ndim != 2 or labels.shape[1] != grid.ny:
        raise ValueError("labels must have shape (nz, ny)")
    if polarization not in ("TE", "TM"):
        raise ValueError("polarization must be 'TE' or 'TM'")
    k = complex(index_background) * grid.k0
    ky = grid.ky
    kz = longitudinal_wavenumber(index_background, grid.k0, ky * ky)
    potential = _potential(indices, labels, index_background, grid.k0)
    # the samples sit at y_n = (n − ny//2)·dy; the shift is the same phase factor for the
    # incident-field product and the field, so the DFT in sample order matches numpy.fft of the
    # field samples
    v_hat = np.fft.fft(potential, axis=1) / grid.ny * np.sinc(ky * grid.dy / (2.0 * np.pi))
    amplitude = _scalar_amplitude(v_hat, kz, k, grid.dz)
    if polarization == "TM":
        amplitude = amplitude * kz * kz / (k * k)
    return np.asarray(amplitude, dtype=np.complex128)


def born_forward_3d(
    grid: Grid3D,
    indices: list[complex],
    labels: IntArray,
    index_background: complex,
    polarization: tuple[complex, complex],
) -> ComplexArray:
    """Plane-wave amplitudes of the first-order scattered transverse field at the output plane, 3D.

    Args:
        grid: transverse grid and slab thickness dz (mm), wavelength (µm).
        indices: complex refractive index per material label.
        labels: material label per pixel, shape (nz, ny, nx).
        index_background: index ñ₀ of the background in which the incident wave travels.
        polarization: transverse amplitudes (e_x, e_y) of the normally incident plane wave.

    Returns:
        Amplitudes fft2(E_s)/(nx·ny) of (E_x, E_y), shape (ny, nx, 2), numpy order with
        (k_x, k_y) = grid.k_grid(), at z_out = nz·dz. The incident wave is not included.
    """
    labels = np.asarray(labels)
    if labels.ndim != 3 or labels.shape[1:] != (grid.ny, grid.nx):
        raise ValueError("labels must have shape (nz, ny, nx)")
    k = complex(index_background) * grid.k0
    kx, ky = grid.k_grid()
    kz = longitudinal_wavenumber(index_background, grid.k0, kx * kx + ky * ky)
    potential = _potential(indices, labels, index_background, grid.k0)
    cell: FloatArray = np.sinc(kx * grid.dx / (2.0 * np.pi)) * np.sinc(ky * grid.dy / (2.0 * np.pi))
    v_hat = np.fft.fft2(potential, axes=(1, 2)) / (grid.nx * grid.ny) * cell
    scalar = _scalar_amplitude(v_hat, kz, k, grid.dz)
    ex, ey = complex(polarization[0]), complex(polarization[1])
    kk = k * k
    out = np.empty((grid.ny, grid.nx, 2), dtype=np.complex128)
    out[..., 0] = scalar * ((1 - kx * kx / kk) * ex - kx * ky / kk * ey)
    out[..., 1] = scalar * (-kx * ky / kk * ex + (1 - ky * ky / kk) * ey)
    return out
