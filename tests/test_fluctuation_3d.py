"""K18 in 3D: Born increment for weak fluctuations in `Bidirectional3D` (B7a, E10, #37).

Specification: docs/spezifikation.md §10. References: the 2D increment (`test_fluctuation.py`,
itself checked against Born and RCWA) for structures invariant along x, and the vector first Born
approximation `bwpm.born.born_forward_3d` (K17, checked against 2D) on a thin random layer.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pytest

from bwpm._types import FloatArray, IntArray
from bwpm.born import born_forward_3d
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.volume3d import Bidirectional3D, Grid3D

LAMBDA_UM = 0.8
N0 = 1.38
DELTA_N = 1e-5


def _random_staircase(
    width_um: float, thickness_um: float, step_um: float, corr_um: float, levels: int, seed: int
) -> tuple[Grid3D, list[complex], IntArray, FloatArray]:
    """Gaussian random field on `levels` index steps (±2σ = ±2·DELTA_N), 4 background slabs.

    Returns the grid, the labelled staircase for the Born reference (label 0 = ñ₀) and the same
    staircase as δn on a homogeneous host for the increment.
    """
    n, nz = round(width_um / step_um), round(thickness_um / step_um)
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((nz, n, n))
    fz, fy, fx = (2 * np.pi * np.fft.fftfreq(m, step_um) for m in (nz, n, n))
    filt = np.exp(
        -(fz[:, None, None] ** 2 + fy[None, :, None] ** 2 + fx[None, None, :] ** 2) * corr_um**2 / 4
    )
    values = np.real(np.fft.ifftn(np.fft.fftn(noise) * filt))
    values /= values.std()
    lab = np.clip(np.round((values + 2) / 4 * (levels - 1)), 0, levels - 1).astype(int) + 1
    steps = np.array([0.0] + [DELTA_N * (-2 + 4 * m / (levels - 1)) for m in range(levels)])
    labels = np.concatenate([lab, np.zeros((4, n, n), dtype=int)])
    grid = Grid3D(n, n, step_um * 1e-3, step_um * 1e-3, step_um * 1e-3, LAMBDA_UM)
    return grid, [complex(N0 + s) for s in steps], labels, steps[labels]


def _x_polarised(grid: Grid3D) -> np.ndarray:  # type: ignore[type-arg]
    field = np.zeros((grid.ny, grid.nx, 2), dtype=np.complex128)
    field[..., 0] = 1.0
    return field


def test_zero_fluctuation_is_bit_identical_3d() -> None:
    """δn ≡ 0 gives exactly the unchanged BWPM, also with iterations (K12–K16 untouched)."""
    rng = np.random.default_rng(18)
    grid = Grid3D(16, 16, 0.05e-3, 0.05e-3, 0.05e-3, LAMBDA_UM)
    labels = rng.integers(0, 3, size=(6, 16, 16))
    indices = [complex(N0), complex(1.5), complex(1.45, 0.01)]
    field = _x_polarised(grid)
    reference = Bidirectional3D(grid, indices, labels, complex(N0)).run(field, iterations=2)
    zero = Bidirectional3D(
        grid, indices, labels, complex(N0), fluctuation=np.zeros(labels.shape)
    ).run(field, iterations=2)
    assert np.array_equal(zero.transmitted, reference.transmitted)
    assert np.array_equal(zero.reflected, reference.reflected)


def test_fluctuation_requires_one_host_per_slab_3d() -> None:
    grid = Grid3D(4, 4, 0.1e-3, 0.1e-3, 0.02e-3, LAMBDA_UM)
    labels = np.zeros((2, 4, 4), dtype=int)
    labels[0, :2] = 1
    fluctuation = np.zeros(labels.shape)
    fluctuation[0, 1:3] = 1e-3
    with pytest.raises(ValueError, match="host"):
        Bidirectional3D(
            grid, [complex(N0), complex(1.5)], labels, complex(N0), fluctuation=fluctuation
        )
    with pytest.raises(ValueError, match="shape"):
        Bidirectional3D(
            grid, [complex(N0), complex(1.5)], labels, complex(N0), fluctuation=np.zeros((2, 4, 5))
        )


@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_increment_3d_reduces_to_2d(pol: Literal["TE", "TM"]) -> None:
    """Grating along y, invariant along x: the 3D increment equals the 2D one (k_x = 0)."""
    px, nx, step = 40, 4, 0.02e-3
    nz = 50
    row = np.zeros(px)
    row[: px // 2] = 0.01
    delta_2d = np.vstack([np.tile(row, (nz, 1)), np.zeros((4, px))])
    labels_2d = np.zeros(delta_2d.shape, dtype=int)
    component = 0 if pol == "TE" else 1
    field_2d = np.zeros((px, 2), dtype=np.complex128)
    field_2d[:, component] = 1.0
    out_2d = Propagator2D(
        grid=Grid2D(ny=px, dy=step, dz=step, wavelength_um=LAMBDA_UM),
        indices=[complex(N0)],
        labels=labels_2d,
        index_before=N0,
        fluctuation=delta_2d,
    ).run(field_2d)
    grid = Grid3D(nx, px, step, step, step, LAMBDA_UM)
    field = np.zeros((px, nx, 2), dtype=np.complex128)
    field[..., component] = 1.0
    out = (
        Bidirectional3D(
            grid,
            [complex(N0)],
            np.repeat(labels_2d[:, :, None], nx, axis=2),
            complex(N0),
            fluctuation=np.repeat(delta_2d[:, :, None], nx, axis=2),
        )
        .run(field, iterations=0)
        .transmitted
    )
    for column in range(nx):
        np.testing.assert_allclose(out[:, column, :], out_2d, rtol=0, atol=1e-12)


def _relative_error_per_bin(step_um: float) -> list[float]:
    """Thin random layer, x-polarised: increment against vector Born per angle bin,
    Σ|E − E_Born|²/Σ|E_Born|² over both transverse components."""
    grid, indices, labels, delta = _random_staircase(3.2, 1.0, step_um, 0.3, 8, seed=3)
    host = np.zeros(labels.shape, dtype=int)
    out = (
        Bidirectional3D(grid, [complex(N0)], host, complex(N0), fluctuation=delta)
        .run(_x_polarised(grid), iterations=0)
        .transmitted
    )
    spectrum = np.fft.fft2(out, axes=(0, 1)) / (grid.nx * grid.ny)
    born = born_forward_3d(grid, indices, labels, N0, (1, 0))
    kx, ky = grid.k_grid()
    k = N0 * grid.k0
    angle = np.degrees(np.arcsin(np.sqrt(np.clip((kx**2 + ky**2) / k**2, 0, 1))))
    errors = []
    for lo, hi in [(2, 20), (20, 40), (40, 60), (60, 75)]:
        sel = (kx**2 + ky**2 < k**2) & (angle >= lo) & (angle < hi)
        diff = np.sum(np.abs(spectrum[sel] - born[sel]) ** 2)
        errors.append(float(np.sqrt(diff / np.sum(np.abs(born[sel]) ** 2))))
    return errors


def test_increment_3d_matches_vector_born() -> None:
    """Up to 75° the field per angle bin matches vector Born within 1 % in amplitude (K17
    tolerance; the pure FPWPM misses by up to 70 % in power), falling when dy = dz is halved."""
    coarse = _relative_error_per_bin(0.04)
    fine = _relative_error_per_bin(0.02)
    for c, f in zip(coarse, fine, strict=True):
        assert c < 1e-2
        assert f < 0.5 * c


def test_grazing_grid_points_stay_bounded() -> None:
    """§10 (found in B7b): on a grid with exactly grazing points (|k⊥| = ñ₀k₀, k_z = i·10⁻¹²·k₀ by
    the rule of §3, #25) a random field must not be amplified; with E_z = −k⊥·Ê⊥/k_z taken from
    those points it grew to 10²⁶ within 11 slabs."""
    rng = np.random.default_rng(2)
    nx, ny, nz = 8, 16, 11
    grid = Grid3D(nx, ny, 0.04e-3, 0.05e-3, 0.03e-3, 0.6)  # k_y = 2·2π/(16·0.05 µm) = 1.5·k₀
    labels = np.ones((nz, ny, nx), dtype=np.int64)
    fluctuation = 0.03 * rng.normal(size=(nz, ny, nx))
    field = rng.normal(size=(ny, nx, 2)) + 1j * rng.normal(size=(ny, nx, 2))
    out = (
        Bidirectional3D(grid, [1.33 + 0j, 1.5 + 0j], labels, 1.5 + 0j, fluctuation=fluctuation)
        .run(field, iterations=0)
        .transmitted
    )
    assert np.sum(np.abs(out) ** 2) < 1.1 * np.sum(np.abs(field) ** 2)
