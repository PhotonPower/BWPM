"""Born scattering at weak index fluctuations δn on top of the structure (stages B7, B8).

Specification: docs/spezifikation.md §10 (forward) and §11 (back scattering); source chen2020,
eqs. (3)–(9), algorithm 1 and supplement eqs. (S1)–(S7), with the linearised potential of E11
(#46) and Gauss–Legendre points for the slab integral (F1, #53) instead of sinc ≈ 1 resp. the
Taylor series; vector projector from the dyadic Green's function (chew1995, checked in K17).
Used by `bwpm.fpwpm.Propagator2D`, `bwpm.bidirectional.Bidirectional2D` and
`bwpm.volume3d.Bidirectional3D`; verified in K18 and K19.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import scipy.fft as sfft
from scipy.sparse.linalg import LinearOperator, gmres

from bwpm._types import ComplexArray, FloatArray, IntArray
from bwpm.interface import GRAZING_TOLERANCE, longitudinal_wavenumber


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


GAUSS_POINTS = ((1 - 1 / np.sqrt(3)) / 2, (1 + 1 / np.sqrt(3)) / 2)  # ζ/dz, weights 1/2 (F1, #53)


def _wave_numbers(
    host: complex, k0: float, kx: FloatArray, ky: FloatArray
) -> tuple[ComplexArray, np.ndarray]:  # type: ignore[type-arg]
    """k_z in the host and the mask of propagating, non-grazing waves (§10, rule of §3, #25)."""
    kp2 = kx * kx + ky * ky
    kz = longitudinal_wavenumber(host, k0, kp2)
    k2 = complex(host) ** 2 * k0 * k0
    propagating = (kp2 < abs(k2)) & (np.abs(k2 - kp2) > GRAZING_TOLERANCE * abs(k2))
    return kz, propagating


def born_sources(
    spectrum: ComplexArray,
    direction: int,
    delta_n: ComplexArray,
    host: complex,
    k0: float,
    dz: float,
    kx: FloatArray,
    ky: FloatArray,
    workers: int | None = None,
) -> ComplexArray:
    """Born sources Ŝ_m = F{V·E_m} of one slab at the Gauss points ζ_m, §11.

    The incident field of direction d (+1: given at the left slab edge, −1: at the right edge) is
    propagated in the host to ζ_m (distance ζ_m resp. dz − ζ_m) with Ê_z = −d·k⊥·Ê⊥/k_z
    (wende2023 eq. (5)); propagating, non-grazing waves only. V = 2ñ₀k₀²δn (E11, #46).
    Fields ∝ exp(i(k·r − ωt)).

    Args:
        spectrum: F{E⊥} of the incident field, shape (…, 2) = (E_x, E_y).
        direction: +1 or −1.
        delta_n: δn per pixel, shape (…).
        host: host index ñ₀ of the slab.
        k0: vacuum wavenumber (1/mm).
        dz: slab thickness (mm).
        kx, ky: transverse frequencies (1/mm) broadcastable to the shape of `delta_n`.
        workers: threads for scipy.fft.

    Returns:
        Ŝ_m, shape (2, …, 3) = (Gauss point, transverse grid, (x, y, z)).
    """
    axes = tuple(range(spectrum.ndim - 1))
    kz, propagating = _wave_numbers(host, k0, kx, ky)
    safe_kz = np.where(propagating, kz, 1.0)
    ez = np.where(
        propagating, -direction * (kx * spectrum[..., 0] + ky * spectrum[..., 1]) / safe_kz, 0.0
    )
    full = np.stack([spectrum[..., 0], spectrum[..., 1], ez], axis=-1)
    potential = 2 * host * k0**2 * delta_n
    sources = []
    for zeta in GAUSS_POINTS:
        distance = zeta * dz if direction > 0 else (1 - zeta) * dz
        phase = np.where(propagating, np.exp(1j * kz * distance), 0.0)
        field = sfft.ifftn(full * phase[..., None], axes=axes, workers=workers)
        sources.append(sfft.fftn(potential[..., None] * field, axes=axes, workers=workers))
    return np.asarray(sources)


def born_emission(
    sources: ComplexArray,
    direction: int,
    host: complex,
    k0: float,
    dz: float,
    kx: FloatArray,
    ky: FloatArray,
    workers: int | None = None,
) -> ComplexArray:
    """Field radiated by the Born sources of one slab in direction d′, at the right slab edge, §11.

    ΔE⊥ = F⁻¹{Σ_m ½·e^(i d′ k_z (dz − ζ_m))·i dz/(2k_z)·[P(k_d′)·Ŝ_m]⊥}, k_d′ = (k_x, k_y, d′k_z),
    P = I − kkᵀ/k² (chew1995); for d′ = −1 the field leaving the left edge, referred back to the
    right edge. Propagating, non-grazing waves only.

    Args:
        sources: Ŝ_m from `born_sources`, shape (2, …, 3).
        direction: d′ = +1 or −1.
        host, k0, dz, kx, ky, workers: as in `born_sources`.

    Returns:
        ΔE⊥ in real space, shape (…, 2).
    """
    axes = tuple(range(sources.ndim - 2))
    k = host * k0
    kz, propagating = _wave_numbers(host, k0, kx, ky)
    safe_kz = np.where(propagating, kz, 1.0)
    signed_kz = direction * kz
    total = np.zeros(sources.shape[1:-1] + (2,), dtype=np.complex128)
    for m, zeta in enumerate(GAUSS_POINTS):
        source: ComplexArray = sources[m]
        green = np.where(
            propagating,
            0.5 * np.exp(1j * signed_kz * (1 - zeta) * dz) * 1j * dz / (2 * safe_kz),
            0.0,
        )
        longitudinal = (
            kx * source[..., 0] + ky * source[..., 1] + signed_kz * source[..., 2]
        ) / k**2
        total += green[..., None] * np.stack(
            [source[..., 0] - kx * longitudinal, source[..., 1] - ky * longitudinal], axis=-1
        )
    return np.asarray(sfft.ifftn(total, axes=axes, workers=workers))


GMRES_TOLERANCE = 1e-12  # relative residual of the Cayley system (G2, tightened from 1e-10)


def _scattering_operator(
    direction: int,
    delta_n: ComplexArray,
    host: complex,
    k0: float,
    dz: float,
    kx: FloatArray,
    ky: FloatArray,
    workers: int | None,
) -> Callable[[ComplexArray], ComplexArray]:
    """G_d of §12 on spectra F{E⊥} (…, 2): i dz/(2k_z)·[P(k_d)·F{V·e}]⊥ at the slab midpoint."""
    k = host * k0
    kz, propagating = _wave_numbers(host, k0, kx, ky)
    safe_kz = np.where(propagating, kz, 1.0)
    green = np.where(propagating, 1j * dz / (2 * safe_kz), 0.0)
    signed_kz = direction * kz
    potential = 2 * host * k0**2 * delta_n

    def apply(spectrum: ComplexArray) -> ComplexArray:
        axes = tuple(range(spectrum.ndim - 1))
        ex = np.where(propagating, spectrum[..., 0], 0.0)
        ey = np.where(propagating, spectrum[..., 1], 0.0)
        ez = np.where(propagating, -direction * (kx * ex + ky * ey) / safe_kz, 0.0)
        field = sfft.ifftn(np.stack([ex, ey, ez], axis=-1), axes=axes, workers=workers)
        source = sfft.fftn(potential[..., None] * field, axes=axes, workers=workers)
        longitudinal = (
            kx * source[..., 0] + ky * source[..., 1] + signed_kz * source[..., 2]
        ) / k**2
        return np.asarray(
            green[..., None]
            * np.stack(
                [source[..., 0] - kx * longitudinal, source[..., 1] - ky * longitudinal], axis=-1
            )
        )

    return apply


def born_step(
    spectrum: ComplexArray,
    direction: int,
    delta_n: ComplexArray,
    host: complex,
    k0: float,
    dz: float,
    kx: FloatArray,
    ky: FloatArray,
    workers: int | None = None,
) -> ComplexArray:
    """Flux-conserving Born increment of one slab (Cayley form at the slab midpoint), §12.

    x = P_½·Ê, (1 − G_d/2)·y′ = G_d·x (GMRES, G2), ΔÊ = P_½^(±)·y′; with the WPM step this gives
    E_out = P_½(1 − G/2)⁻¹(1 + G/2)P_½·E_in in a slab of host, which conserves the flux
    Σ k_z(|E⊥|² + |E_z|²) for real δn (E13, docs/b9-entwurf.md). For direction −1 the increment
    leaves the left edge and is referred back to the right edge, as R2 in §11.

    Args:
        spectrum: F{E⊥} of the incident field (+1: at the left slab edge, −1: at the right edge),
            shape (…, 2).
        direction: +1 or −1.
        delta_n, host, k0, dz, kx, ky, workers: as in `born_sources`.

    Returns:
        ΔE⊥ in real space at the right slab edge, shape (…, 2).
    """
    kz, propagating = _wave_numbers(host, k0, kx, ky)
    half = np.where(propagating, np.exp(1j * kz * dz / 2), 0.0)[..., None]
    operator = _scattering_operator(direction, delta_n, host, k0, dz, kx, ky, workers)
    x = spectrum * half
    rhs = operator(x)
    shape = rhs.shape
    system = LinearOperator(
        (rhs.size, rhs.size),
        matvec=lambda v: v - 0.5 * operator(v.reshape(shape)).ravel(),
        dtype=np.complex128,
    )
    solution, info = gmres(
        system, rhs.ravel(), x0=rhs.ravel(), rtol=GMRES_TOLERANCE, atol=0.0, restart=30,
        maxiter=100,
    )  # fmt: skip
    if info != 0:
        raise RuntimeError(f"GMRES did not converge in the Born step (info {info})")
    increment = solution.reshape(shape)
    if direction > 0:
        increment = increment * half
    else:
        increment = increment * np.where(propagating, np.exp(-1j * kz * dz / 2), 0.0)[..., None]
    axes = tuple(range(spectrum.ndim - 1))
    return np.asarray(sfft.ifftn(increment, axes=axes, workers=workers))


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
    """Forward Born increment of one slab (§10 with the Gauss points of §11): the field entering at
    the left edge (`spectrum` = F{E⊥}, shape (…, 2)) radiates forward to the right edge."""
    sources = born_sources(spectrum, +1, delta_n, host, k0, dz, kx, ky, workers)
    return born_emission(sources, +1, host, k0, dz, kx, ky, workers)
