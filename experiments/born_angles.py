"""K17 / #29: FPWPM against the first Born approximation, per scattering angle.

Prints two tables (Markdown) for docs/testfaelle.md:

1. Lamellar grating, one period, diffraction order +1 at 2–63°: complex ratio FPWPM/Born of the
   transverse field (TE: E_x, TM: E_y) for several steps dy = dz.
2. Thin random layer (Gaussian random field with Gaussian correlation, quantised to M index
   levels, fixed seed): scattered power per angle bin, FPWPM over Born, for several steps.

Run: python experiments/born_angles.py  (about a minute on one core)
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from bwpm._types import ComplexArray, IntArray
from bwpm.born import born_forward_2d
from bwpm.fpwpm import Grid2D, Propagator2D

LAMBDA_UM = 0.8
N0 = 1.38
DELTA_N = 1e-5
THICKNESS_UM = 2.0
BACKGROUND_SLABS = 4
SEED = 29
POLARIZATIONS: tuple[Literal["TE", "TM"], ...] = ("TE", "TM")


def _run(
    grid: Grid2D, indices: list[complex], labels: IntArray, pol: Literal["TE", "TM"]
) -> tuple[ComplexArray, ComplexArray]:
    """Spectrum (fft/ny) of the FPWPM output and the Born amplitudes, transverse component."""
    comp = 0 if pol == "TE" else 1
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, comp] = 1.0
    out = Propagator2D(grid=grid, indices=indices, labels=labels, index_before=N0).run(field)
    spectrum = np.fft.fft(out[:, comp]) / grid.ny
    born = born_forward_2d(grid, indices, labels, N0, pol)
    return spectrum, born


def lamellar_table() -> None:
    print(
        "| Periode (µm) | Winkel | Schritt (µm) | TE: Betrag | TE: Phase (rad) "
        "| TM: Betrag | TM: Phase (rad) |"
    )
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for period in (20.0, 4.0, 1.2, 0.8, 0.65):
        angle = np.degrees(np.arcsin(LAMBDA_UM / (N0 * period)))
        for step in (0.04, 0.02, 0.01):
            px = round(period / step)
            grid = Grid2D(ny=px, dy=period / px * 1e-3, dz=step * 1e-3, wavelength_um=LAMBDA_UM)
            row = np.zeros(px, dtype=int)
            row[: px // 2] = 1
            labels = np.vstack(
                [
                    np.tile(row, (round(THICKNESS_UM / step), 1)),
                    np.zeros((BACKGROUND_SLABS, px), dtype=int),
                ]
            )
            cells = []
            for pol in POLARIZATIONS:
                spectrum, born = _run(grid, [complex(N0), complex(N0 + DELTA_N)], labels, pol)
                ratio = spectrum[1] / born[1]
                cells += [f"{abs(ratio):.4f}", f"{np.angle(ratio):+.4f}"]
            print(f"| {period} | {angle:.1f}° | {step} | " + " | ".join(cells) + " |")


def random_layer(
    step: float, corr_um: float, levels: int, width_um: float = 40.0
) -> tuple[Grid2D, list[complex], IntArray]:
    """Gaussian random field δn (std DELTA_N, correlation length corr_um) on a grid with step."""
    ny, nz = round(width_um / step), round(THICKNESS_UM / step)
    grid = Grid2D(ny=ny, dy=step * 1e-3, dz=step * 1e-3, wavelength_um=LAMBDA_UM)
    # the same continuous realisation for every step: white noise on the finest grid is not
    # comparable across steps, so draw the spectrum on a fixed set of frequencies instead
    rng = np.random.default_rng(SEED)
    nyc, nzc = round(width_um / 0.01), round(THICKNESS_UM / 0.01)
    noise = rng.standard_normal((nzc, nyc)) + 1j * rng.standard_normal((nzc, nyc))
    fy = np.fft.fftfreq(nyc, 0.01)
    fz = np.fft.fftfreq(nzc, 0.01)
    filt = np.exp(-((2 * np.pi) ** 2) * (fy[None, :] ** 2 + fz[:, None] ** 2) * corr_um**2 / 4)
    coeff = noise * filt
    y = (np.arange(ny) + 0.5) * step
    z = (np.arange(nz) + 0.5) * step
    sel_y = np.abs(fy) < 1 / (2 * step)
    sel_z = np.abs(fz) < 1 / (2 * step)
    c = coeff[np.ix_(sel_z, sel_y)]
    ey = np.exp(2j * np.pi * fy[sel_y][:, None] * y[None, :])
    ez = np.exp(2j * np.pi * fz[sel_z][None, :] * z[:, None])
    values = np.real(ez @ c @ ey)
    values /= values.std()
    # M levels between −2σ and +2σ (pixel centre value, E4)
    lab = np.clip(np.round((values + 2) / 4 * (levels - 1)), 0, levels - 1).astype(int)
    lab = lab + 1  # label 0 is the background
    indices = [complex(N0)] + [
        complex(N0 + DELTA_N * (-2 + 4 * m / (levels - 1))) for m in range(levels)
    ]
    labels = np.vstack([lab, np.zeros((BACKGROUND_SLABS, ny), dtype=int)])
    return grid, indices, labels


def random_table() -> None:
    bins = [(2, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 75)]
    print("| l_c (µm) | Schritt (µm) | Pol. | " + " | ".join(f"{a}–{b}°" for a, b in bins) + " |")
    print("| --- | --- | --- | " + " | ".join("---" for _ in bins) + " |")
    for corr in (0.2, 0.5):
        for step in (0.04, 0.02):
            grid, indices, labels = random_layer(step, corr, 16)
            angle = np.degrees(np.arcsin(np.clip(grid.ky / (grid.k0 * N0), -1, 1)))
            propagating = np.abs(grid.ky) < grid.k0 * N0
            for pol in POLARIZATIONS:
                spectrum, born = _run(grid, indices, labels, pol)
                cells = []
                for lo, hi in bins:
                    sel = propagating & (np.abs(angle) >= lo) & (np.abs(angle) < hi)
                    ratio = np.sum(np.abs(spectrum[sel]) ** 2) / np.sum(np.abs(born[sel]) ** 2)
                    cells.append(f"{ratio:.3f}")
                print(f"| {corr} | {step} | {pol} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    print(f"λ = {LAMBDA_UM} µm, ñ₀ = {N0}, δn = {DELTA_N}, Schicht {THICKNESS_UM} µm\n")
    lamellar_table()
    print()
    random_table()
