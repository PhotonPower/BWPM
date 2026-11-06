"""B4d (K15): CUDA core (bwpm._core.CudaSolver) against the reference files and the CPU core.

Runs only where the core was built with CUDA (BWPM_ENABLE_CUDA=ON) and a GPU is present; with
BWPM_REQUIRE_CUDA=1 a missing CUDA backend fails instead of skipping (large machine).
GPU tolerance: cuFFT and the device math library round differently from pocketfft/libm; the
reference files are compared with rel 1e-10 and abs 1e-11 of the maximum.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from bwpm import fast
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.results import load_json
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d

ROOT = Path(__file__).resolve().parents[1]

if not fast.cuda_available():
    if os.environ.get("BWPM_REQUIRE_CUDA") == "1":
        raise RuntimeError("BWPM_REQUIRE_CUDA=1 but the CUDA backend is not available")
    pytest.skip("CUDA backend not built or no GPU", allow_module_level=True)


def _load(name: str) -> dict:  # type: ignore[type-arg]
    return load_json((ROOT / "reference" / f"{name}.result.json").read_text("utf-8"))[1]


def _close(a: np.ndarray, b: np.ndarray) -> None:  # type: ignore[type-arg]
    assert a.shape == b.shape
    assert np.allclose(a, b, rtol=1e-10, atol=1e-11 * float(np.max(np.abs(b)) or 1.0))


def _propagator_2d(d: dict) -> Propagator2D:  # type: ignore[type-arg]
    p = d["parameters"]
    grid = Grid2D(p["ny"], p["dy_mm"], p["dz_mm"], p["wavelength_um"])
    return Propagator2D(
        grid,
        list(d["inputs"]["indices"]),
        d["inputs"]["labels"],
        d["inputs"]["index_before"],
        d["inputs"].get("mask"),
    )


def test_fpwpm2d_interface() -> None:
    d = _load("fpwpm2d_interface")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], 0, keep=True, backend="cuda")
    assert res.forward is not None
    _close(res.forward, d["expected"]["planes"])


def test_bwpm2d_tilted_silver() -> None:
    d = _load("bwpm2d_tilted_silver")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], 2, keep=True, backend="cuda")
    assert res.forward is not None and res.backward is not None
    _close(res.forward, d["expected"]["forward"])
    _close(res.backward, d["expected"]["backward"])


def test_bwpm2d_slab_exact() -> None:
    d = _load("bwpm2d_slab_exact")
    res = fast.run_2d(_propagator_2d(d), d["inputs"]["field"], 40, backend="cuda")
    _close(res.transmitted, d["expected"]["transmitted"])
    _close(res.reflected, d["expected"]["reflected"])


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
    res = fast.run_3d(solver, d["inputs"]["field"], 2, (0, 10, nz), backend="cuda")
    _close(res.transmitted, d["expected"]["transmitted"])
    _close(res.reflected, d["expected"]["reflected"])


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_equals_cpu_core(seed: int) -> None:
    rng = np.random.default_rng(seed)
    nx, ny, nz = (15, 12, 9) if seed % 2 else (8, 16, 11)
    grid = Grid3D(nx, ny, 0.04e-3, 0.05e-3, 0.03e-3, 0.6)
    labels = rng.integers(0, 4, size=(nz, 1, 1)) * np.ones((1, ny, nx), dtype=np.int64)
    labels[nz // 2] = rng.integers(0, 4, size=(ny, nx))
    indices = [1.0 + 0j, 1.5 + 0j, 2.1 + 0.01j, 0.213 + 6.697j]
    field = rng.normal(size=(ny, nx, 2)) + 1j * rng.normal(size=(ny, nx, 2))
    kx, ky = grid.k_grid()
    spectrum = np.fft.fft2(field, axes=(0, 1))
    spectrum[np.sqrt(kx**2 + ky**2) >= 1.2 * grid.k0] = 0.0
    field = np.fft.ifft2(spectrum, axes=(0, 1))
    mask = absorbing_mask_2d(grid, 2 * grid.dx) if seed == 2 else None
    solver = Bidirectional3D(grid, indices, labels, 1.2 + 0j, mask)
    cpu = fast.run_3d(solver, field, 2, (0, 5, nz), backend="cpu")
    gpu = fast.run_3d(solver, field, 2, (0, 5, nz), backend="cuda")
    scale = 1e-10 * float(np.max(np.abs(field)))
    assert np.max(np.abs(gpu.transmitted - cpu.transmitted)) < scale
    assert np.max(np.abs(gpu.reflected - cpu.reflected)) < scale
    assert gpu.stored_values == cpu.stored_values
