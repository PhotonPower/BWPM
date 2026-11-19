"""K17: forward scattering of weak index modulations against the first Born approximation (#29).

The reference `bwpm.born.born_forward_2d` is itself checked against two independent exact
solutions: the Airy formula of a homogeneous slab and an RCWA of a lamellar grating
(`tests/rcwa_reference.py`); in both cases the difference must vanish to first order in δn.
Then the FPWPM is compared with it per diffraction order (TE and TM).

Setting as in Scattuccino S9b-2: λ = 0.8 µm, background ñ₀ = 1.38. Behind the modulated layer
follow a few background slabs, so the field is compared in the background medium (inside the
layer the forward field lacks the reflected part and differs from the total field by the Fresnel
factor 1 − δn/(2ñ₀), which is of first order).
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pytest
from rcwa_reference import lamellar_transmission

from bwpm._types import IntArray
from bwpm.born import born_forward_2d
from bwpm.fpwpm import Grid2D, Propagator2D

LAMBDA_UM = 0.8
N0 = 1.38
THICKNESS_UM = 2.0
BACKGROUND_SLABS = 4
POLARIZATIONS: tuple[Literal["TE", "TM"], ...] = ("TE", "TM")


def _lamellar(
    period_um: float, px: int, step_um: float, dy_um: float | None = None
) -> tuple[Grid2D, IntArray]:
    """One period of a lamellar grating: label 1 on pixels 0 … px/2 − 1, then background slabs."""
    dy = period_um / px if dy_um is None else dy_um
    if not np.isclose(px * dy, period_um):
        raise ValueError("period must be a whole number of pixels")
    grid = Grid2D(ny=px, dy=dy * 1e-3, dz=step_um * 1e-3, wavelength_um=LAMBDA_UM)
    row = np.zeros(px, dtype=int)
    row[: px // 2] = 1
    nz = round(THICKNESS_UM / step_um)
    labels = np.vstack([np.tile(row, (nz, 1)), np.zeros((BACKGROUND_SLABS, px), dtype=int)])
    return grid, labels


def _fpwpm_over_born(
    period_um: float, px: int, step_um: float, delta_n: float, pol: Literal["TE", "TM"]
) -> complex:
    """Ratio FPWPM / Born of diffraction order +1 at the output plane."""
    grid, labels = _lamellar(period_um, px, step_um)
    indices = [complex(N0), complex(N0 + delta_n)]
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, 0 if pol == "TE" else 1] = 1.0
    out = Propagator2D(grid=grid, indices=indices, labels=labels, index_before=N0).run(field)
    spectrum = np.fft.fft(out[:, 0 if pol == "TE" else 1]) / grid.ny
    born = born_forward_2d(grid, indices, labels, N0, pol)
    return complex(spectrum[1] / born[1])


def _angle_deg(period_um: float) -> float:
    return float(np.degrees(np.arcsin(LAMBDA_UM / (N0 * period_um))))


@pytest.mark.parametrize("pol", POLARIZATIONS)
def test_born_slab_matches_airy_to_first_order(pol: Literal["TE", "TM"]) -> None:
    """Homogeneous slab ñ₀ + δn: Born = exact (Airy) transmission − incident wave, error O(δn)."""
    grid = Grid2D(ny=8, dy=0.1e-3, dz=0.02e-3, wavelength_um=LAMBDA_UM)
    nz = round(THICKNESS_UM / 0.02)
    k0 = grid.k0
    length = nz * grid.dz
    errors = []
    for delta_n in (1e-5, 1e-6):
        n1 = N0 + delta_n
        r01 = (N0 - n1) / (N0 + n1)  # normal incidence: TE and TM agree for the transverse field
        t01, t10 = 2 * N0 / (N0 + n1), 2 * n1 / (N0 + n1)
        phase = np.exp(1j * n1 * k0 * length)
        exact = t01 * t10 * phase / (1 - r01 * r01 * phase * phase)
        scattered = exact - np.exp(1j * N0 * k0 * length)
        labels = np.ones((nz, grid.ny), dtype=int)
        born = born_forward_2d(grid, [complex(N0), complex(n1)], labels, N0, pol)
        assert np.allclose(born[1:], 0.0, atol=1e-18)
        errors.append(abs(born[0] / scattered - 1))
    assert errors[0] < 1e-3
    assert 5 < errors[0] / errors[1] < 20  # first order in δn


@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("period_um", [4.0, 1.2, 0.8, 0.65])
def test_born_grating_matches_rcwa_to_first_order(
    period_um: float, pol: Literal["TE", "TM"]
) -> None:
    """Lamellar grating, orders ±1 at 8–63°: Born against RCWA, difference O(δn)."""
    px = 40
    dy = period_um / px * 1e-3
    errors = []
    for delta_n in (1e-5, 1e-6):
        grid, labels = _lamellar(period_um, px, 0.02)
        labels = labels[: round(THICKNESS_UM / 0.02)]  # layer only, output plane at z = L
        born = born_forward_2d(grid, [complex(N0), complex(N0 + delta_n)], labels, N0, pol)
        assert np.all(np.isfinite(born))  # evanescent orders included
        # pixel n is the cell |y − n·dy| ≤ dy/2 in the sample order of the FFT
        rcwa = lamellar_transmission(
            N0, N0, N0 + delta_n, period_um * 1e-3, -dy / 2, (px // 2 - 0.5) * dy,
            THICKNESS_UM * 1e-3, LAMBDA_UM * 1e-3, 10, pol,
        )  # fmt: skip
        errors.append(max(abs(born[1] / rcwa[1] - 1), abs(born[-1] / rcwa[-1] - 1)))
    assert errors[0] < 1e-3
    assert 5 < errors[0] / errors[1] < 20


@pytest.mark.parametrize("pol", POLARIZATIONS)
def test_fpwpm_matches_born_paraxially(pol: Literal["TE", "TM"]) -> None:
    """Grating period 40 µm (first order at 0.83°): FPWPM = Born within 10⁻³."""
    ratio = _fpwpm_over_born(40.0, 400, 0.02, 1e-5, pol)
    assert abs(ratio - 1) < 1e-3


@pytest.mark.xfail(
    strict=True,
    reason="limit of the method (#29): FPWPM scatters too little into large angles (TE) and too "
    "much in E_y (TM); docs/testfaelle.md, K17",
)
@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("period_um", [1.2, 0.8, 0.65])
def test_fpwpm_matches_born_at_large_angles(period_um: float, pol: Literal["TE", "TM"]) -> None:
    """First order at 29°, 46° and 63°: FPWPM = Born within 1 %."""
    ratio = _fpwpm_over_born(period_um, 40, 0.02, 1e-5, pol)
    assert abs(ratio - 1) < 1e-2, f"{_angle_deg(period_um):.1f}°: FPWPM/Born = {ratio:.4f}"


def test_large_angle_deviation_converges_in_step() -> None:
    """The deviation at 46° (TE) is a property of the method, not of the sampling.

    dy = dz = 0.04, 0.02, 0.01 µm: the ratio FPWPM/Born converges (successive differences fall
    at least like first order) and stays far from 1.
    """
    ratios = [_fpwpm_over_born(0.8, px, 0.8 / px, 1e-5, "TE") for px in (20, 40, 80)]
    first, second = abs(ratios[1] - ratios[0]), abs(ratios[2] - ratios[1])
    assert second < 0.7 * first
    assert abs(ratios[2] - 1) > 10 * second


@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("period_um", [4.0, 1.2, 0.8, 0.65])
def test_fpwpm_deviation_matches_first_order_analysis(
    period_um: float, pol: Literal["TE", "TM"]
) -> None:
    """Cause of the deviation (docs/spezifikation.md §9): first-order analysis of the FPWPM step.

    TE: FPWPM/Born = cos θ·(3 − cos θ)/2 (obliquity of the incident instead of the scattered wave
    in the bulk, cos θ, plus the Fresnel factors at entry and exit); TM: E_y as TE without the
    factor cos²θ, (3 − cos θ)/(2 cos θ). The phase error is first order in the step, so the ratio
    is extrapolated from dy = dz = 0.02 and 0.01 µm (Richardson).
    """
    c = float(np.cos(np.radians(_angle_deg(period_um))))
    expected = c * (3 - c) / 2 if pol == "TE" else (3 - c) / (2 * c)
    coarse = _fpwpm_over_born(period_um, round(period_um / 0.02), 0.02, 1e-5, pol)
    fine = _fpwpm_over_born(period_um, round(period_um / 0.01), 0.01, 1e-5, pol)
    assert abs((2 * fine - coarse) / expected - 1) < 1e-3
