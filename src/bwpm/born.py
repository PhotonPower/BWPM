"""First Born approximation of the forward-scattered field of a weak 2D volume (reference for K17).

Specification: docs/spezifikation.md, section 9; sources bornwolf1999 (first Born approximation),
chew1995 (spectral form of the scalar and dyadic Green's function). Both are cited
bibliographically, so no equation numbers are given here.

Geometry as in `bwpm.fpwpm`: slab j of the labelled volume occupies z ∈ [j·dz, (j + 1)·dz]
(z = 0 at the input plane), pixel (j, n) is the cell |y − y_n| ≤ dy/2 with constant index. The
volume is embedded in a homogeneous background ñ₀ and periodic in y with period ny·dy. A plane
wave of unit amplitude, E = e^(i k z) with k = ñ₀k₀, falls in along +z; TE means E = E_x, TM means
E = E_y (2D, k_x = 0). Fields ∝ exp(i(k·r − ωt)).

With the potential V = k₀²(ñ² − ñ₀²), the scattered field obeys (∇² + k²)E_s = −V·E_inc (TE). Its
plane-wave amplitude with k_y = k_y,m behind all scatterers (z ≥ z′) is

    E_s,m(z) = i/(2k_z) · Σ_j V̂_j,m · ∫_slab j e^(i k_z (z − z′)) e^(i k z′) dz′,

k_z = sqrt(k² − k_y²) with Im k_z ≥ 0. V̂_j,m is the exact Fourier coefficient of the pixel
staircase: the DFT of the pixel values times sinc(k_y dy/2). For TM the dyadic Green's function
adds the factor 1 − k_y²/k² = k_z²/k² to the transverse component E_y.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from bwpm._types import ComplexArray, IntArray
from bwpm.fpwpm import Grid2D
from bwpm.interface import longitudinal_wavenumber


def born_forward_2d(
    grid: Grid2D,
    indices: list[complex],
    labels: IntArray,
    index_background: complex,
    polarization: Literal["TE", "TM"],
) -> ComplexArray:
    """Plane-wave amplitudes of the first-order scattered transverse field at the output plane.

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
    k0 = grid.k0
    k = complex(index_background) * k0
    ky = grid.ky
    kz = longitudinal_wavenumber(index_background, k0, ky * ky)
    nz = labels.shape[0]
    z_out = nz * grid.dz

    index_map = np.asarray(indices, dtype=np.complex128)[labels]
    potential = k0 * k0 * (index_map**2 - complex(index_background) ** 2)  # 1/mm², (nz, ny)
    # Fourier coefficient of the pixel staircase; the samples sit at y_n = (n − ny//2)·dy, the
    # shift is the same phase factor for the incident-field product and the field, so the DFT in
    # sample index order matches numpy.fft of the field samples
    v_hat = np.fft.fft(potential, axis=1) / grid.ny * np.sinc(ky * grid.dy / (2.0 * np.pi))

    # ∫_slab j e^(i k_z (z_out − z′)) e^(i k z′) dz′
    #   = e^(i k z_j) · e^(i k_z (z_out − z_j − dz)) · (e^(i k dz) − e^(i k_z dz)) / (i (k − k_z)),
    # written so that every factor stays bounded for evanescent orders (Im k_z > 0)
    z_left = np.arange(nz) * grid.dz
    mismatch = k - kz
    small = np.abs(mismatch * grid.dz) < 1e-8
    safe = np.where(small, 1.0, mismatch)
    within = np.where(
        small,
        grid.dz * np.exp(0.5j * (k + kz) * grid.dz),
        (np.exp(1j * k * grid.dz) - np.exp(1j * kz * grid.dz)) / (1j * safe),
    )
    travel = np.exp(1j * kz[None, :] * (z_out - z_left[:, None] - grid.dz))
    slab_integral = np.exp(1j * k * z_left)[:, None] * travel * within[None, :]
    amplitude = 1j / (2.0 * kz) * np.sum(v_hat * slab_integral, axis=0)
    if polarization == "TM":
        amplitude = amplitude * kz * kz / (k * k)
    return np.asarray(amplitude, dtype=np.complex128)
