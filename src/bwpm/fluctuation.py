"""Born increment for weak index fluctuations δn on top of the structure (stage B7, decision E10).

Specification: docs/spezifikation.md §10; source chen2020, eqs. (3)–(9) and algorithm 1 (scalar,
sinc ≈ 1, fields at slab midpoints) with the linearised potential of decision E11 (#46), vector
projector from the dyadic Green's function (chew1995,
checked in K17). Used by `bwpm.fpwpm.Propagator2D` (2D, k_x = 0) and
`bwpm.volume3d.Bidirectional3D` (3D); verified in K18.
"""

from __future__ import annotations

import numpy as np
import scipy.fft as sfft

from bwpm._types import ComplexArray, FloatArray, IntArray
from bwpm.interface import longitudinal_wavenumber


def slab_hosts(
    indices: list[complex], labels: IntArray, fluctuation: FloatArray | ComplexArray
) -> list[complex | None]:
    """Host index ñ₀ per slab: the one structure index under all pixels with δn ≠ 0 (§10).

    Args:
        indices: complex index per material label.
        labels: label per pixel, shape (nz, …).
        fluctuation: δn per pixel, same shape as `labels`.

    Returns:
        per slab the host index, or None if δn ≡ 0 in that slab.
    """
    delta = np.asarray(fluctuation)
    if delta.shape != labels.shape:
        raise ValueError("fluctuation must have the shape of labels")
    hosts: list[complex | None] = []
    for j in range(labels.shape[0]):
        found = {complex(indices[int(label)]) for label in labels[j][delta[j] != 0]}
        if len(found) > 1:
            raise ValueError(f"slab {j}: fluctuation must lie in one host index")
        hosts.append(found.pop() if found else None)
    return hosts


def born_increment(
    spectrum: ComplexArray,
    delta_n: ComplexArray,
    host: complex,
    k0: float,
    dz: float,
    kx: FloatArray,
    ky: FloatArray,
    workers: int | None = None,
) -> ComplexArray:
    """Vector Born increment of one slab in real space, §10.

    ΔE⊥ = F⁻¹{e^(i k_z dz/2)·i dz/(2k_z)·[P(k)·F{V·E_m}]⊥} with V = 2ñ₀k₀²δn (E11, #46; chen2020
    eq. (1) has k₀²((ñ₀ + δn)² − ñ₀²), equal to first order),
    E_m = F⁻¹{e^(i k_z dz/2)(Ê⊥, Ê_z)} the full field at the slab midpoint in the host,
    Ê_z = −k⊥·Ê⊥/k_z (wende2023 eq. (5)), P = I − kkᵀ/k²; propagating waves only (|k⊥| < k).
    Fields ∝ exp(i(k·r − ωt)).

    Args:
        spectrum: F{E⊥} at the slab entry over the transverse axes, shape (…, 2) = (E_x, E_y).
        delta_n: δn per pixel, shape (…).
        host: host index ñ₀ of the slab.
        k0: vacuum wavenumber (1/mm).
        dz: slab thickness (mm).
        kx, ky: transverse frequencies (1/mm) broadcastable to the shape of `delta_n`.
        workers: threads for scipy.fft.

    Returns:
        ΔE⊥ in real space, shape of `spectrum`.
    """
    axes = tuple(range(spectrum.ndim - 1))
    k = host * k0
    kp2 = kx * kx + ky * ky
    kz = longitudinal_wavenumber(host, k0, kp2)
    half = np.exp(1j * kz * dz / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ez = np.where(kz != 0, -(kx * spectrum[..., 0] + ky * spectrum[..., 1]) / kz, 0.0)
    full = np.stack([spectrum[..., 0], spectrum[..., 1], ez], axis=-1) * half[..., None]
    middle = sfft.ifftn(full, axes=axes, workers=workers)
    potential = 2 * host * k0**2 * delta_n
    source = sfft.fftn(potential[..., None] * middle, axes=axes, workers=workers)
    propagating = kp2 < np.real(k * np.conj(k))
    green = np.where(propagating, half * 1j * dz / (2 * np.where(propagating, kz, 1.0)), 0.0)
    longitudinal = (kx * source[..., 0] + ky * source[..., 1] + kz * source[..., 2]) / k**2
    projected = np.stack(
        [source[..., 0] - kx * longitudinal, source[..., 1] - ky * longitudinal], axis=-1
    )
    return np.asarray(sfft.ifftn(green[..., None] * projected, axes=axes, workers=workers))
