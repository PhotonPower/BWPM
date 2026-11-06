"""B3 measurement run: 3D BWPM of the fibre facet with index-matching cap (wende2023 §2.2.2).

Measures run time, peak memory and the number of stored reflection values for several grid
sizes, and the 3D Strehl ratio on the facet (wende2023: 0.97) plus the reflected power.

Usage (on the large machine, results go to experiments/out/):

    python experiments/b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2

FFTs use all cores by default (`--workers -1`). `--fractions f` sets dx = dy = dz =
λ/(f·n_max); wende2023 uses f = 5. The objective field is the Richards-Wolf plane-wave
spectrum of an aplanatic objective (NA 0.4, x-polarised pupil),
amplitude ∝ cos^(−1/2) θ per (k_x, k_y), with the transverse vector factors
(cos²φ cos θ + sin²φ, cos φ sin φ (cos θ − 1)).
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from bwpm import fast
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d, memory_estimate

LAMBDA_UM = 0.635
NA = 0.4
N_CORE, N_CLAD, N_IMC = 1.5, float(np.sqrt(1.5**2 - 0.39**2)), 1.511
D_IMC = 50e-3
AIR_FOCUS_DEPTH = 32.42e-3  # from the 2D pilot; the 3D best focus is reported, not imposed


def objective_field(grid: Grid3D, distance: float) -> np.ndarray:  # type: ignore[type-arg]
    """Field (ny, nx, 2) at the input plane of a beam focusing `distance` behind it in air."""
    kx, ky = grid.k_grid()
    k = grid.k0
    kp = np.sqrt(kx * kx + ky * ky)
    inside = kp <= NA * k
    kz = np.sqrt(np.maximum(k * k - kp * kp, 0.0))
    cos_t = kz / k
    phi = np.arctan2(ky, kx)
    # cos θ^(-1/2) only inside the NA (outside, cos θ can be 0 and the value is not used)
    weight = np.zeros_like(cos_t)
    weight[inside] = cos_t[inside] ** -0.5
    amp = weight * np.exp(-1j * kz * distance)
    ex = amp * (np.cos(phi) ** 2 * cos_t + np.sin(phi) ** 2)
    ey = amp * np.cos(phi) * np.sin(phi) * (cos_t - 1.0)
    shift = np.exp(-1j * (kx * (grid.nx // 2) * grid.dx + ky * (grid.ny // 2) * grid.dy))
    field = np.stack([np.fft.ifft2(ex * shift), np.fft.ifft2(ey * shift)], axis=-1)
    return np.asarray(field)


def volume(grid: Grid3D, with_imc: bool) -> tuple[list[complex], np.ndarray, float]:  # type: ignore[type-arg]
    front = -D_IMC if with_imc else 0.0
    z_in = front - 20e-3
    nz = round((15e-3 - z_in) / grid.dz)
    z = z_in + (np.arange(nz) + 0.5) * grid.dz
    rr = np.sqrt(grid.x[None, :] ** 2 + grid.y[:, None] ** 2)
    labels = np.zeros((nz, grid.ny, grid.nx), dtype=np.int8)
    for j, zj in enumerate(z):
        if with_imc and -D_IMC < zj <= 0:
            labels[j] = 1
        elif zj > 0:
            labels[j] = np.where(rr <= 1e-3, 2, 3)
    return [1.0 + 0j, N_IMC + 0j, N_CORE + 0j, N_CLAD + 0j], labels, z_in


def count_pairs(labels: np.ndarray) -> int:  # type: ignore[type-arg]
    """Number of distinct (label before, label after) pairs between neighbouring slabs.

    Plane by plane with bincount (O(N), little memory); np.unique(axis=0) over all ~10⁸–10⁹
    pairs at once took > 40 min on the server (report in PR #23).
    """
    m = int(labels.max()) + 1
    seen = np.zeros(m * m, dtype=bool)
    for j in range(labels.shape[0] - 1):
        code = labels[j].astype(np.int64).ravel() * m + labels[j + 1].astype(np.int64).ravel()
        seen |= np.bincount(code, minlength=m * m).astype(bool)
    return int(seen.sum())


def strehl(spectrum: np.ndarray, grid: Grid3D, k: float) -> tuple[float, float]:  # type: ignore[type-arg]
    kx, ky = grid.k_grid()
    support = np.sqrt(kx * kx + ky * ky) <= NA * grid.k0
    a = spectrum[support]
    kz = np.sqrt(k * k - kx[support] ** 2 - ky[support] ** 2)
    norm = float(np.sum(np.abs(a))) ** 2
    offsets = np.arange(-5e-3, 5e-3, 0.02e-3)
    values = [abs(np.sum(a * np.exp(1j * kz * d))) ** 2 / norm for d in offsets]
    j = int(np.argmax(values))
    return float(values[j]), float(offsets[j])


def run(
    fraction: float,
    window_um: float,
    iterations: int,
    with_imc: bool,
    workers: int,
    core: bool = False,
    backend: str = "auto",
) -> dict[str, float]:
    d = LAMBDA_UM * 1e-3 / (fraction * N_IMC)
    n = round(window_um * 1e-3 / d)
    n += n % 2
    grid = Grid3D(n, n, d, d, d, LAMBDA_UM)
    indices, labels, z_in = volume(grid, with_imc)
    distance = -z_in - D_IMC + AIR_FOCUS_DEPTH if with_imc else -z_in
    field = objective_field(grid, distance)
    plane = round((-1e-3 - z_in) / d)
    solver = Bidirectional3D(
        grid,
        indices,
        labels,
        1.0 + 0j,
        absorbing_mask_2d(grid, 20 * LAMBDA_UM * 1e-3),
        workers=workers,
    )
    start = time.perf_counter()
    if core:
        res = fast.run_3d(
            solver,
            field,
            iterations,
            (plane,),
            threads=workers if workers > 0 else None,
            backend=backend,
        )
    else:
        res = solver.run(field, iterations=iterations, keep=(plane,))
    seconds = time.perf_counter() - start
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # Linux: kB
    p_in = float(np.sum(np.abs(field) ** 2))
    p_r = float(np.sum(np.abs(res.reflected) ** 2))  # ratio of |E|² sums (same medium, rough)
    forward = res.planes[plane][0][..., 0]
    shift = np.exp(1j * (grid.k_grid()[0] * (n // 2) * d + grid.k_grid()[1] * (n // 2) * d))
    s, best = (
        strehl(np.fft.fft2(forward) * shift, grid, grid.k0 * N_IMC) if with_imc else (0.0, 0.0)
    )
    pairs = count_pairs(labels)
    return {
        "fraction": fraction,
        "n": n,
        "nz": labels.shape[0],
        "with_imc": float(with_imc),
        "seconds": seconds,
        "peak_rss_mb": peak_mb,
        "stored_values": float(res.stored_values),
        "model_mb": memory_estimate(n, n, res.stored_values // 2, pairs, nz=labels.shape[0])
        / 2**20,
        "workers": float(workers),
        "core": float(core and fast.CORE_AVAILABLE),
        "cuda": float(  # backend actually used (run_3d raises if CUDA was requested but missing)
            core and backend != "cpu" and fast.cuda_available()
        ),
        "reflected_ratio": p_r / p_in,
        "strehl": s,
        "best_focus_um": (best - 1e-3) * 1e3,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fractions", type=float, nargs="+", default=[2.5, 3.5, 5.0])
    parser.add_argument("--window-um", type=float, default=96.0)
    parser.add_argument("--iterations", type=int, default=2)
    parser.add_argument("--bare", action="store_true", help="fibre without cap")
    parser.add_argument(
        "--core", action="store_true", help="C++ core (bwpm._core) instead of NumPy"
    )
    parser.add_argument("--backend", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument(
        "--workers", type=int, default=-1, help="FFT threads (scipy.fft), -1 = all cores"
    )
    args = parser.parse_args()
    out = Path(__file__).parent / "out"
    out.mkdir(exist_ok=True)
    results = []
    for f in args.fractions:
        r = run(
            f, args.window_um, args.iterations, not args.bare, args.workers, args.core, args.backend
        )
        print(json.dumps(r), flush=True)
        results.append(r)
    name = (
        f"b3_scaling_{'bare' if args.bare else 'imc'}{'_' + args.backend if args.core else ''}.json"
    )
    (out / name).write_text(json.dumps({"argv": sys.argv[1:], "results": results}, indent=2))


if __name__ == "__main__":
    main()
