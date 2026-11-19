"""K18: weak index fluctuations δn with the Born increment (B7a, decision E10, #37).

Specification: docs/spezifikation.md §10. The reference values are independent of the code under
test: the first Born approximation of the pixel staircase (`bwpm.born.born_forward_2d`, itself
checked against Airy and RCWA in K17) for weak δn, and the RCWA of the lamellar grating
(`tests/rcwa_reference.py`) for multiple scattering at δn = 0.01–0.03. The settings and
tolerances of K17 (`test_born.py`) are reused unchanged.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pytest
from rcwa_reference import lamellar_transmission

from bwpm._types import FloatArray, IntArray
from bwpm.born import born_forward_2d
from bwpm.fpwpm import Grid2D, Propagator2D, power_flux

LAMBDA_UM = 0.8
N0 = 1.38
THICKNESS_UM = 2.0
BACKGROUND_SLABS = 4
POLARIZATIONS: tuple[Literal["TE", "TM"], ...] = ("TE", "TM")


def _component(pol: Literal["TE", "TM"]) -> int:
    return 0 if pol == "TE" else 1


def _grating(
    period_um: float, px: int, step_um: float, thickness_um: float, delta_n: float
) -> tuple[Grid2D, IntArray, FloatArray]:
    """Lamellar grating as δn on a homogeneous host (label 0), then background slabs."""
    grid = Grid2D(ny=px, dy=period_um / px * 1e-3, dz=step_um * 1e-3, wavelength_um=LAMBDA_UM)
    nz = round(thickness_um / step_um)
    labels = np.zeros((nz + BACKGROUND_SLABS, px), dtype=int)
    fluctuation = np.zeros((nz + BACKGROUND_SLABS, px))
    fluctuation[:nz, : px // 2] = delta_n
    return grid, labels, fluctuation


def _run(
    grid: Grid2D, labels: IntArray, fluctuation: FloatArray | None, pol: Literal["TE", "TM"]
) -> np.ndarray:  # type: ignore[type-arg]
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, _component(pol)] = 1.0
    prop = Propagator2D(
        grid=grid, indices=[complex(N0)], labels=labels, index_before=N0, fluctuation=fluctuation
    )
    return np.asarray(prop.run(field))


def _increment_over_born(
    period_um: float, px: int, step_um: float, delta_n: float, pol: Literal["TE", "TM"]
) -> complex:
    """Ratio (FPWPM with Born increment) / Born of diffraction order +1 at the output plane."""
    grid, labels, fluctuation = _grating(period_um, px, step_um, THICKNESS_UM, delta_n)
    out = _run(grid, labels, fluctuation, pol)
    spectrum = np.fft.fft(out[:, _component(pol)]) / grid.ny
    structure = (fluctuation != 0).astype(int)  # the same staircase as labels for the reference
    born = born_forward_2d(grid, [complex(N0), complex(N0 + delta_n)], structure, N0, pol)
    return complex(spectrum[1] / born[1])


def test_zero_fluctuation_is_bit_identical() -> None:
    """δn = None and δn ≡ 0 give exactly the result of the unchanged method (K1–K16)."""
    rng = np.random.default_rng(18)
    grid = Grid2D(ny=64, dy=0.05e-3, dz=0.05e-3, wavelength_um=LAMBDA_UM)
    labels = rng.integers(0, 3, size=(20, grid.ny))
    indices = [complex(N0), complex(1.5), complex(1.45, 0.01)]
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((grid.y / 0.6e-3) ** 2))
    field[:, 1] = 0.5j * field[:, 0]
    reference = Propagator2D(grid=grid, indices=indices, labels=labels, index_before=N0).run(field)
    zero = Propagator2D(
        grid=grid,
        indices=indices,
        labels=labels,
        index_before=N0,
        fluctuation=np.zeros(labels.shape),
    ).run(field)
    assert np.array_equal(zero, reference)


def test_fluctuation_requires_one_host_per_slab() -> None:
    """δn only on pixels whose structure index is the same in the whole slab (§10)."""
    grid = Grid2D(ny=8, dy=0.1e-3, dz=0.02e-3, wavelength_um=LAMBDA_UM)
    labels = np.zeros((2, grid.ny), dtype=int)
    labels[0, :4] = 1
    fluctuation = np.zeros(labels.shape)
    fluctuation[0, 2:6] = 1e-3  # spans labels 1 and 0 in slab 0
    with pytest.raises(ValueError, match="host"):
        Propagator2D(
            grid=grid,
            indices=[complex(N0), complex(1.5)],
            labels=labels,
            index_before=N0,
            fluctuation=fluctuation,
        )
    with pytest.raises(ValueError, match="shape"):
        Propagator2D(
            grid=grid,
            indices=[complex(N0), complex(1.5)],
            labels=labels,
            index_before=N0,
            fluctuation=np.zeros((3, grid.ny)),
        )


@pytest.mark.parametrize("pol", POLARIZATIONS)
def test_increment_matches_born_paraxially(pol: Literal["TE", "TM"]) -> None:
    """Grating period 40 µm (0.83°): within 10⁻³ as in K17."""
    ratio = _increment_over_born(40.0, 400, 0.02, 1e-5, pol)
    assert abs(ratio - 1) < 1e-3


@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("period_um", [1.2, 0.8, 0.65])
def test_increment_matches_born_at_large_angles(period_um: float, pol: Literal["TE", "TM"]) -> None:
    """First order at 29°, 46° and 63°: within 1 % (the K17 tolerance that the pure FPWPM misses),
    and the error falls when dy = dz is halved (midpoint rule and pixel sampling, §10)."""
    coarse = _increment_over_born(period_um, round(period_um / 0.02), 0.02, 1e-5, pol)
    fine = _increment_over_born(period_um, round(period_um / 0.01), 0.01, 1e-5, pol)
    assert abs(coarse - 1) < 1e-2
    assert abs(fine - 1) < 0.5 * abs(coarse - 1)


def _orders_against_rcwa(
    delta_n: float, pol: Literal["TE", "TM"], fluctuation_model: bool
) -> float:
    """Max relative error of orders 0, ±1 against the RCWA.

    Grating period 0.8 µm (orders ±1 at 46°), thickness 5 µm, dy = dz = 0.02 µm. With
    `fluctuation_model=False` the grating is a structure (pure WPM, K17) instead of δn.
    """
    period_um, px, step_um, thickness_um = 0.8, 40, 0.02, 5.0
    grid, labels, fluctuation = _grating(period_um, px, step_um, thickness_um, delta_n)
    if fluctuation_model:
        out = _run(grid, labels, fluctuation, pol)
    else:
        field = np.zeros((grid.ny, 2), dtype=np.complex128)
        field[:, _component(pol)] = 1.0
        structure = (fluctuation != 0).astype(int)
        out = Propagator2D(
            grid=grid,
            indices=[complex(N0), complex(N0 + delta_n)],
            labels=structure,
            index_before=N0,
        ).run(field)
    # back to the exit face z = L through the background slabs (homogeneous host)
    kz = np.sqrt((N0 * grid.k0) ** 2 - grid.ky**2 + 0j)
    back = np.exp(-1j * kz * BACKGROUND_SLABS * grid.dz)
    exit_field = np.fft.ifft(np.fft.fft(out, axis=0) * back[:, None], axis=0)
    spectrum = np.fft.fft(exit_field[:, _component(pol)]) / grid.ny
    dy = grid.dy
    rcwa = lamellar_transmission(
        N0, N0, N0 + delta_n, period_um * 1e-3, -dy / 2, (px // 2 - 0.5) * dy,
        thickness_um * 1e-3, LAMBDA_UM * 1e-3, 15, pol,
    )  # fmt: skip
    return float(max(abs(spectrum[m] / rcwa[m] - 1) for m in (0, 1, -1)))


@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("delta_n", [0.01, 0.03])
def test_multiple_scattering_matches_rcwa(delta_n: float, pol: Literal["TE", "TM"]) -> None:
    """Lamellar grating δn = 0.01 and 0.03 over 5 µm against the RCWA.

    Without reflections a forward model misses the Fresnel factors at entry and exit, of order
    δn/ñ₀ together (§10); the tolerance is δn. The pure WPM misses orders ±1 by 20–70 % (K17).
    """
    error = _orders_against_rcwa(delta_n, pol, fluctuation_model=True)
    error_wpm = _orders_against_rcwa(delta_n, pol, fluctuation_model=False)
    assert error < delta_n
    assert error < error_wpm / 5


@pytest.mark.parametrize("pol", POLARIZATIONS)
@pytest.mark.parametrize("delta_n", [0.01, 0.03])
def test_linearised_potential_halves_the_error(delta_n: float, pol: Literal["TE", "TM"]) -> None:
    """Decision E11 (#46): with V = 2ñ₀k₀²δn the mean phase is exact to second order at normal
    incidence, so the error of k₀δn²L/(2ñ₀) in order 0 (§10) disappears and the error against
    the RCWA stays below δn/2 (with V = k₀²((ñ₀ + δn)² − ñ₀²) of chen2020: up to 0.7·δn)."""
    assert _orders_against_rcwa(delta_n, pol, fluctuation_model=True) < delta_n / 2


@pytest.mark.parametrize("pol", POLARIZATIONS)
def test_flux_is_conserved_to_first_order_in_step(pol: Literal["TE", "TM"]) -> None:
    """Lossless δn = 0.01 over 2 µm: the flux with E_z (eq. (5)) is conserved up to the growth
    O((k₀δn dz)²) per slab (§10), which falls like dz."""
    deviations = []
    for px, step_um in ((40, 0.02), (80, 0.01)):
        grid, labels, fluctuation = _grating(0.8, px, step_um, THICKNESS_UM, 0.01)
        out = _run(grid, labels, fluctuation, pol)
        incident = np.zeros((grid.ny, 2), dtype=np.complex128)
        incident[:, _component(pol)] = 1.0
        deviations.append(power_flux(out, grid, N0) / power_flux(incident, grid, N0) - 1)
    assert abs(deviations[0]) < 1e-3
    assert abs(deviations[1]) < 0.6 * abs(deviations[0])
