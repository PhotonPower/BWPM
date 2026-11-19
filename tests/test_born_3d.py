"""K17 in 3D: vector Born reference and the FPWPM in 3D (#29).

The 3D reference `bwpm.born.born_forward_3d` is checked against the 2D reference (itself checked
against Airy and RCWA, `tests/test_born.py`) for gratings invariant along x or y, and its
transverse projector against the 2D factors in the frame of each scattered wave (eigenvalue
k_z²/k² along k⊥, 1 across). Then the forward FPWPM in 3D is compared on a thin random layer.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm._types import IntArray
from bwpm.born import born_forward_2d, born_forward_3d
from bwpm.fpwpm import Grid2D
from bwpm.volume3d import Bidirectional3D, Grid3D

LAMBDA_UM = 0.8
N0 = 1.38
DELTA_N = 1e-5
STEP_UM = 0.04


def _random_layer(
    width_um: float, thickness_um: float, corr_um: float, levels: int, seed: int
) -> tuple[Grid3D, list[complex], IntArray]:
    """Gaussian random field (Gaussian correlation), ±2σ on `levels` labels, 4 background slabs."""
    n, nz = round(width_um / STEP_UM), round(thickness_um / STEP_UM)
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((nz, n, n))
    fz, fy, fx = (2 * np.pi * np.fft.fftfreq(m, STEP_UM) for m in (nz, n, n))
    filt = np.exp(
        -(fz[:, None, None] ** 2 + fy[None, :, None] ** 2 + fx[None, None, :] ** 2) * corr_um**2 / 4
    )
    values = np.real(np.fft.ifftn(np.fft.fftn(noise) * filt))
    values /= values.std()
    lab = np.clip(np.round((values + 2) / 4 * (levels - 1)), 0, levels - 1).astype(int) + 1
    indices = [complex(N0)] + [
        complex(N0 + DELTA_N * (-2 + 4 * m / (levels - 1))) for m in range(levels)
    ]
    labels = np.concatenate([lab, np.zeros((4, n, n), dtype=int)])
    grid = Grid3D(n, n, STEP_UM * 1e-3, STEP_UM * 1e-3, STEP_UM * 1e-3, LAMBDA_UM)
    return grid, indices, labels


@pytest.mark.parametrize("axis", ["y", "x"])
def test_born_3d_reduces_to_2d(axis: str) -> None:
    """A grating varying only along y (or x) gives the 2D result: TE for the field across the
    grating vector, TM (factor k_z²/k²) for the field along it."""
    period_px, nx_other, nz = 40, 4, 50
    step = 0.02e-3
    row = np.zeros(period_px, dtype=int)
    row[: period_px // 2] = 1
    labels_2d = np.tile(row, (nz, 1))
    indices = [complex(N0), complex(N0 + DELTA_N)]
    grid_2d = Grid2D(ny=period_px, dy=step, dz=step, wavelength_um=LAMBDA_UM)
    te = born_forward_2d(grid_2d, indices, labels_2d, N0, "TE")
    tm = born_forward_2d(grid_2d, indices, labels_2d, N0, "TM")
    if axis == "y":
        grid = Grid3D(nx_other, period_px, step, step, step, LAMBDA_UM)
        labels = np.repeat(labels_2d[:, :, None], nx_other, axis=2)
        along_x = born_forward_3d(grid, indices, labels, N0, (1, 0))[:, 0, :]  # k_x = 0 column
        along_y = born_forward_3d(grid, indices, labels, N0, (0, 1))[:, 0, :]
        np.testing.assert_allclose(along_x[:, 0], te, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(along_y[:, 1], tm, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(along_x[:, 1], 0, atol=1e-18)
    else:
        grid = Grid3D(period_px, nx_other, step, step, step, LAMBDA_UM)
        labels = np.repeat(labels_2d[:, None, :], nx_other, axis=1)
        across = born_forward_3d(grid, indices, labels, N0, (0, 1))[0, :, :]  # k_y = 0 row
        along = born_forward_3d(grid, indices, labels, N0, (1, 0))[0, :, :]
        np.testing.assert_allclose(across[:, 1], te, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(along[:, 0], tm, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(along[:, 1], 0, atol=1e-18)


def test_born_3d_projector_in_the_frame_of_each_wave() -> None:
    """For every propagating k⊥: the component along û = k⊥/|k⊥| carries k_z²/k² (TM), the one
    along v̂ = ẑ × û carries 1 (TE), both times the same scalar amplitude (x-polarised input)."""
    grid, indices, labels = _random_layer(3.2, 0.4, 0.3, 4, seed=5)
    field = born_forward_3d(grid, indices, labels, N0, (1, 0))
    kx, ky = grid.k_grid()
    k = N0 * grid.k0
    kp = np.hypot(kx, ky)
    sel = (kp > 0) & (kp < 0.98 * k) & (np.abs(ky) > 0.1 * kp)
    ux, uy = kx[sel] / kp[sel], ky[sel] / kp[sel]
    e_u = ux * field[sel, 0] + uy * field[sel, 1]
    e_v = -uy * field[sel, 0] + ux * field[sel, 1]
    scalar = e_v / (-uy)  # ê·v̂ = −u_y for ê = x̂
    np.testing.assert_allclose(e_u, (1 - kp[sel] ** 2 / k**2) * ux * scalar, rtol=1e-10)


def test_fpwpm_3d_is_componentwise_scalar_with_k17_factor() -> None:
    """Thin random layer, x-polarised: the FPWPM scatters into E_x only (E_y < 10⁻⁶ in power) with
    the amplitude r(θ)·S, r = cos θ (3 − cos θ)/2 from the 2D analysis (§9), S the scalar Born
    amplitude; power per angle bin within 2 % of the prediction up to 60° (dy = dz = 0.04 µm)."""
    grid, indices, labels = _random_layer(6.4, 2.0, 0.3, 8, seed=3)
    field = np.zeros((grid.ny, grid.nx, 2), dtype=np.complex128)
    field[..., 0] = 1.0
    out = Bidirectional3D(grid, indices, labels, complex(N0)).run(field, iterations=0).transmitted
    spectrum = np.fft.fft2(out, axes=(0, 1)) / (grid.nx * grid.ny)
    born = born_forward_3d(grid, indices, labels, N0, (1, 0))
    spectrum[0, 0] = born[0, 0] = 0
    kx, ky = grid.k_grid()
    k = N0 * grid.k0
    cos_t = np.sqrt(np.clip(1 - (kx**2 + ky**2) / k**2, 0, 1))
    angle = np.degrees(np.arccos(cos_t))
    scalar = born[..., 0] / (1 - kx**2 / k**2)
    r = cos_t * (3 - cos_t) / 2
    for lo, hi in [(2, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60)]:
        sel = (kx**2 + ky**2 < k**2) & (angle >= lo) & (angle < hi)
        measured = np.sum(np.abs(spectrum[sel, 0]) ** 2) / np.sum(np.abs(scalar[sel]) ** 2)
        predicted = np.sum(r[sel] ** 2 * np.abs(scalar[sel]) ** 2) / np.sum(
            np.abs(scalar[sel]) ** 2
        )
        assert abs(measured / predicted - 1) < 0.02, (lo, hi, measured, predicted)
        assert np.sum(np.abs(spectrum[sel, 1]) ** 2) < 1e-6 * np.sum(np.abs(spectrum[sel, 0]) ** 2)
