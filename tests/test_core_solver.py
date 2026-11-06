"""B4c (K14): CPU core (bwpm._core.Solver) against the reference files and the Python reference.

Skipped without the core unless BWPM_REQUIRE_CORE=1 (CI).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from bwpm import fast
from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.results import load_json
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d

ROOT = Path(__file__).resolve().parents[1]

if not fast.CORE_AVAILABLE:
    if os.environ.get("BWPM_REQUIRE_CORE") == "1":
        raise RuntimeError("BWPM_REQUIRE_CORE=1 but bwpm._core is not built")
    pytest.skip("bwpm._core not built", allow_module_level=True)


def _load(name: str) -> dict:  # type: ignore[type-arg]
    return load_json((ROOT / "reference" / f"{name}.result.json").read_text("utf-8"))[1]


def _close(a: np.ndarray, b: np.ndarray, tol: dict) -> None:  # type: ignore[type-arg]
    assert a.shape == b.shape
    assert np.allclose(a, b, rtol=tol["rel"], atol=tol["abs"])


def _propagator_2d(d: dict) -> Propagator2D:  # type: ignore[type-arg]
    p = d["parameters"]
    grid = Grid2D(p["ny"], p["dy_mm"], p["dz_mm"], p["wavelength_um"])
    mask = d["inputs"].get("mask")
    return Propagator2D(
        grid, list(d["inputs"]["indices"]), d["inputs"]["labels"], d["inputs"]["index_before"], mask
    )


def test_fpwpm2d_interface() -> None:
    d = _load("fpwpm2d_interface")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], iterations=0, keep=True)
    assert res.forward is not None
    _close(res.forward, d["expected"]["planes"], d["tolerance"])


def test_bwpm2d_tilted_silver() -> None:
    d = _load("bwpm2d_tilted_silver")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], iterations=2, keep=True)
    assert res.forward is not None and res.backward is not None
    _close(res.forward, d["expected"]["forward"], d["tolerance"])
    _close(res.backward, d["expected"]["backward"], d["tolerance"])
    assert np.allclose(
        res.energy_residual, d["expected"]["energy_residual"], rtol=1e-10, atol=1e-12
    )


def test_bwpm2d_slab_exact() -> None:
    d = _load("bwpm2d_slab_exact")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], iterations=40)
    _close(res.transmitted, d["expected"]["transmitted"], d["tolerance"])
    _close(res.reflected, d["expected"]["reflected"], d["tolerance"])
    r = complex(np.mean(res.reflected[:, 0] / d["inputs"]["field"][:, 0]))
    assert abs(r - d["reference"]["r_te_airy"]) < d["reference_tolerance"]["abs"]


def test_bwpm3d_tilted() -> None:
    d = _load("bwpm3d_tilted")
    p = d["parameters"]
    grid = Grid3D(p["nx"], p["ny"], p["dx_mm"], p["dy_mm"], p["dz_mm"], p["wavelength_um"])
    solver = Bidirectional3D(
        grid,
        list(d["inputs"]["indices"]),
        d["inputs"]["labels"],
        d["inputs"]["index_before"],
        d["inputs"]["mask"],
    )
    nz = d["inputs"]["labels"].shape[0]
    res = fast.run_3d(solver, d["inputs"]["field"], iterations=2, keep=(0, 10, nz))
    _close(res.transmitted, d["expected"]["transmitted"], d["tolerance"])
    _close(res.reflected, d["expected"]["reflected"], d["tolerance"])
    fwd = np.stack([res.planes[q][0] for q in (0, 10, nz)])
    bwd = np.stack([res.planes[q][1] for q in (0, 10, nz)])
    _close(fwd, d["expected"]["forward_planes"], d["tolerance"])
    _close(bwd, d["expected"]["backward_planes"], d["tolerance"])


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_random_volumes_equal_python(seed: int) -> None:
    # random 3D volumes with four materials (incl. silver), random field, odd and even sizes
    rng = np.random.default_rng(seed)
    nx, ny, nz = (15, 12, 9) if seed % 2 else (8, 16, 11)
    grid = Grid3D(nx, ny, 0.04e-3, 0.05e-3, 0.03e-3, 0.6)
    labels = rng.integers(0, 4, size=(nz, 1, 1)) * np.ones((1, ny, nx), dtype=np.int64)
    labels[nz // 2] = rng.integers(0, 4, size=(ny, nx))
    indices = [1.0 + 0j, 1.5 + 0j, 2.1 + 0.01j, 0.213 + 6.697j]
    field = rng.normal(size=(ny, nx, 2)) + 1j * rng.normal(size=(ny, nx, 2))
    # propagating spectrum only: white noise would excite evanescent waves near the surface
    # plasmon pole of silver, where |r| is large and rounding differences are amplified
    kx, ky = grid.k_grid()
    spectrum = np.fft.fft2(field, axes=(0, 1))
    spectrum[np.sqrt(kx**2 + ky**2) >= 1.2 * grid.k0] = 0.0
    field = np.fft.ifft2(spectrum, axes=(0, 1))
    mask = absorbing_mask_2d(grid, 2 * grid.dx) if seed == 2 else None
    ref = Bidirectional3D(grid, indices, labels, 1.2 + 0j, mask).run(field, 2, keep=(0, 5, nz))
    core = fast.run_3d(Bidirectional3D(grid, indices, labels, 1.2 + 0j, mask), field, 2, (0, 5, nz))
    # dx ≠ dy puts grid points within rounding of |k⊥| = n·k₀, where k_z = sqrt(n²k₀² − k⊥²)
    # cancels; the core therefore rounds k exactly like 2π·numpy.fft.fftfreq (found here: with
    # 2π·m/(N d) the results differed by 2e-10).
    scale = 1e-12 * float(np.max(np.abs(field)))
    assert np.max(np.abs(core.transmitted - ref.transmitted)) < scale
    assert np.max(np.abs(core.reflected - ref.reflected)) < scale
    for q in (0, 5, nz):
        for j in (0, 1):
            assert np.max(np.abs(core.planes[q][j] - ref.planes[q][j])) < scale
    assert core.stored_values == ref.stored_values


def test_threads_do_not_change_the_result() -> None:
    d = _load("bwpm2d_tilted_silver")
    one = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], iterations=2, threads=1)
    many = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], iterations=2, threads=4)
    assert np.max(np.abs(one.reflected - many.reflected)) < 1e-14


def test_2d_core_equals_reference_class() -> None:
    d = _load("bwpm2d_tilted_silver")
    prop = _propagator_2d(d)
    ref = Bidirectional2D(prop).run(d["inputs"]["field"], iterations=3)
    core = fast.run_2d(prop, d["inputs"]["field"], iterations=3)
    assert np.max(np.abs(core.transmitted - ref.transmitted)) < 1e-12
    assert np.allclose(core.energy_residual, ref.energy_residual, rtol=1e-10, atol=1e-13)
