"""K20: flux-conserving Born step (B9, E13, docs/spezifikation.md §12, docs/b9-entwurf.md).

References independent of the step: the flux Σ k_z(|E⊥|² + |E_z|²) of the plane-wave spectrum
(E_z from wende2023 eq. (5), written out here), the exact attenuation e^(−2k₀κL) of a plane wave in
a medium with Im ñ = κ, and the instability measured in Scattuccino (#29) for the explicit step.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm.volume3d import Bidirectional3D, Grid3D

N0 = 1.38


def _medium(n: int, nz: int, h_um: float, sigma: float, seed: int) -> np.ndarray:  # type: ignore[type-arg]
    """Whittle–Matérn field, D = 4, l_c = 0.5 µm, as in #29 (Scattuccino S9b-8)."""
    rng = np.random.default_rng(seed)
    fz, fy, fx = (2 * np.pi * np.fft.fftfreq(m, h_um) for m in (nz, n, n))
    q2 = fz[:, None, None] ** 2 + fy[None, :, None] ** 2 + fx[None, None, :] ** 2
    values = np.real(np.fft.ifftn(np.fft.fftn(rng.standard_normal((nz, n, n))) / (1 + 0.25 * q2)))
    return np.asarray(sigma * values / values.std())


def _medium_z(n: int, nz: int, h_um: float, dz_um: float, sigma: float, seed: int) -> np.ndarray:  # type: ignore[type-arg]
    """As `_medium` with a different sampling in z (dz_um) than across (h_um)."""
    rng = np.random.default_rng(seed)
    fz = 2 * np.pi * np.fft.fftfreq(nz, dz_um)
    fy, fx = (2 * np.pi * np.fft.fftfreq(n, h_um) for _ in range(2))
    q2 = fz[:, None, None] ** 2 + fy[None, :, None] ** 2 + fx[None, None, :] ** 2
    values = np.real(np.fft.ifftn(np.fft.fftn(rng.standard_normal((nz, n, n))) / (1 + 0.25 * q2)))
    return np.asarray(sigma * values / values.std())


def _flux(field: np.ndarray, grid: Grid3D) -> float:  # type: ignore[type-arg]
    """Σ Re(k_z)(|Ê⊥|² + |Ê_z|²) over propagating plane waves in the host, E_z = −k⊥·Ê⊥/k_z."""
    kx, ky = grid.k_grid()
    k2 = (N0 * grid.k0) ** 2
    propagating = kx**2 + ky**2 < k2 * (1 - 1e-12)
    kz = np.sqrt(np.where(propagating, k2 - kx**2 - ky**2, 1.0))
    spectrum = np.fft.fft2(field, axes=(0, 1))
    ez = -(kx * spectrum[..., 0] + ky * spectrum[..., 1]) / kz
    power = np.sum(np.abs(spectrum) ** 2, axis=-1) + np.abs(ez) ** 2
    return float(np.sum(np.where(propagating, kz * power, 0.0)))


def _forward(
    sigma: float, h_um: float, depth_um: float, width_um: float = 4.8
) -> tuple[Grid3D, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    n, nz = round(width_um / h_um), round(depth_um / h_um)
    grid = Grid3D(n, n, h_um * 1e-3, h_um * 1e-3, h_um * 1e-3, 0.8)
    field = np.zeros((n, n, 2), dtype=np.complex128)
    field[..., 0] = 1.0
    solver = Bidirectional3D(
        grid,
        [complex(N0)],
        np.zeros((nz, n, n), dtype=np.int64),
        complex(N0),
        fluctuation=_medium(n, nz, h_um, sigma, seed=1),
    )
    return grid, field, solver.run(field, iterations=0).transmitted


@pytest.mark.parametrize(("sigma", "h_um"), [(0.0572, 0.1), (0.0572, 0.05), (0.1, 0.1)])
def test_flux_is_conserved_for_strong_real_fluctuations(sigma: float, h_um: float) -> None:
    """Medium of #29 (σ_n = 0.057 and 0.1, l_c = 0.5 µm), 3 µm, forward only: the flux stays at
    1 within 10⁻⁸ (limited by the GMRES tolerance 10⁻¹⁰ per slab). The explicit step of B7/B8 grew
    it by 4–19 % here (docs/b9-entwurf.md); Σ|E⊥|² itself is not conserved."""
    grid, field, out = _forward(sigma, h_um, 3.0)
    assert abs(_flux(out, grid) / _flux(field, grid) - 1) < 1e-8


def test_scattering_near_grazing_does_not_depend_on_the_step() -> None:
    """σ_n = 0.057 over 3 µm, one medium (lateral grid 0.1 µm, generated at dz = 0.05 µm and
    averaged over pairs of slabs for dz = 0.1 µm): the share of |E⊥|² at cos θ < 0.3 agrees for
    both steps within 30 % (the explicit step of B7/B8 diverged with depth at dz = 0.1 µm)."""
    n, h_um, depth_um = 48, 0.1, 3.0
    fine = _medium_z(n, round(depth_um / 0.05), h_um, 0.05, 0.0572, seed=1)
    shares = []
    for dz_um, delta in ((0.1, 0.5 * (fine[0::2] + fine[1::2])), (0.05, fine)):
        grid = Grid3D(n, n, h_um * 1e-3, h_um * 1e-3, dz_um * 1e-3, 0.8)
        field = np.zeros((n, n, 2), dtype=np.complex128)
        field[..., 0] = 1.0
        out = (
            Bidirectional3D(
                grid,
                [complex(N0)],
                np.zeros(delta.shape, dtype=np.int64),
                complex(N0),
                fluctuation=delta,
            )
            .run(field, iterations=0)
            .transmitted
        )
        kx, ky = grid.k_grid()
        cos_t = np.sqrt(np.clip(1 - (kx**2 + ky**2) / (N0 * grid.k0) ** 2, 0, 1))
        power = np.sum(np.abs(np.fft.fft2(out, axes=(0, 1))) ** 2, axis=-1)
        shares.append(float(np.sum(power[cos_t < 0.3]) / np.sum(power)))
    assert 0.0 < shares[0] < 0.1
    assert abs(shares[0] / shares[1] - 1) < 0.3


def test_absorbing_fluctuation_attenuates_exactly() -> None:
    """Uniform δn = iκ (κ = 10⁻³) over 5 µm, normal incidence: the flux falls as e^(−2k₀κL) to
    10⁻⁶ relative (Cayley factor (1 − b)/(1 + b), b = k₀κ·dz/2, error O(b³) per slab)."""
    n, nz, h_um, kappa = 8, 100, 0.05, 1e-3
    grid = Grid3D(n, n, h_um * 1e-3, h_um * 1e-3, h_um * 1e-3, 0.8)
    field = np.zeros((n, n, 2), dtype=np.complex128)
    field[..., 0] = 1.0
    solver = Bidirectional3D(
        grid,
        [complex(N0)],
        np.zeros((nz, n, n), dtype=np.int64),
        complex(N0),
        fluctuation=np.full((nz, n, n), 1j * kappa),
    )
    out = solver.run(field, iterations=0).transmitted
    expected = np.exp(-2 * grid.k0 * kappa * nz * grid.dz)
    assert abs(_flux(out, grid) / _flux(field, grid) / expected - 1) < 1e-6
