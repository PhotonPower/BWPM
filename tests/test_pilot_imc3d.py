"""B5a: analysis helpers of the 3D fibre-facet pilot (K6 in 3D).

Reference values are analytic and independent of `bwpm.pilot_imc3d`: the cos θ law of the flux of
a plane wave, Fresnel at normal incidence and at the Brewster angle, and the focus position of a
spectrum with a known defocus phase.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm.pilot_imc3d import (
    axial_focus,
    fresnel_reflectances,
    plane_power,
    reference_reflectance,
)
from bwpm.volume3d import Grid3D

LAMBDA_UM = 0.635
STEP = LAMBDA_UM * 1e-3 / 4
GRID = Grid3D(nx=64, ny=64, dx=STEP, dy=STEP, dz=STEP, wavelength_um=LAMBDA_UM)


def _mode(mx: int, my: int, ex: complex, ey: complex) -> np.ndarray:  # type: ignore[type-arg]
    """Transverse field (ny, nx, 2) of one plane-wave mode with the given amplitudes."""
    phase = np.exp(
        1j
        * 2
        * np.pi
        * (mx * np.arange(GRID.nx)[None, :] / GRID.nx + my * np.arange(GRID.ny)[:, None] / GRID.ny)
    )
    return np.stack([ex * phase, ey * phase], axis=-1)


def _angle(mx: int, my: int, n: float) -> float:
    kp = 2 * np.pi * np.hypot(mx / (GRID.nx * GRID.dx), my / (GRID.ny * GRID.dy))
    return float(np.arcsin(kp / (n * GRID.k0)))


@pytest.mark.parametrize("polarisation", ["s", "p"])
def test_plane_power_follows_cosine_law(polarisation: str) -> None:
    """A plane wave with |E| = 1 carries a flux ∝ n cos θ, for s and p polarisation."""
    n = 1.511
    ref = plane_power(_mode(0, 0, 1.0, 0.0), GRID, n)
    mx, my = 6, 0
    theta = _angle(mx, my, n)
    if polarisation == "s":  # E along y, k in the x-z plane: |E_t| = 1
        field = _mode(mx, my, 0.0, 1.0)
    else:  # E in the x-z plane: |E_t| = cos θ
        field = _mode(mx, my, np.cos(theta), 0.0)
    assert plane_power(field, GRID, n) / ref == pytest.approx(np.cos(theta), rel=1e-12)


def test_plane_power_scales_with_index_and_amplitude() -> None:
    f = _mode(0, 0, 1.0, 0.0)
    assert plane_power(2 * f, GRID, 1.0) / plane_power(f, GRID, 1.0) == pytest.approx(4.0)
    assert plane_power(f, GRID, 1.5) / plane_power(f, GRID, 1.0) == pytest.approx(1.5)


def test_fresnel_normal_incidence_and_brewster() -> None:
    n1, n2 = 1.0, 1.511
    rs, rp = fresnel_reflectances(n1, n2, np.array([1.0]))
    expected = ((n2 - n1) / (n2 + n1)) ** 2
    assert rs[0] == pytest.approx(expected, rel=1e-14)
    assert rp[0] == pytest.approx(expected, rel=1e-14)
    cos_b = np.cos(np.arctan(n2 / n1))
    assert fresnel_reflectances(n1, n2, np.array([cos_b]))[1][0] == pytest.approx(0.0, abs=1e-28)
    assert fresnel_reflectances(n1, n2, np.array([1e-9]))[0][0] == pytest.approx(1.0, abs=1e-8)


def test_reference_reflectance_weights_by_polarisation() -> None:
    """A p wave at the Brewster angle of the mode grid is not reflected; s is."""
    n1, n2 = 1.0, 1.511
    mx = 10
    theta = _angle(mx, 0, n1)
    p_wave = _mode(mx, 0, np.cos(theta), 0.0)
    s_wave = _mode(mx, 0, 0.0, 1.0)
    cos_t = np.cos(theta)
    rs, rp = fresnel_reflectances(n1, n2, np.array([cos_t]))
    assert reference_reflectance(s_wave, GRID, n1, n2) == pytest.approx(rs[0], rel=1e-12)
    assert reference_reflectance(p_wave, GRID, n1, n2) == pytest.approx(rp[0], rel=1e-12)
    mixed = p_wave + s_wave  # powers add (orthogonal polarisations, same k): |E|² = cos²θ·… + 1
    ps, ss = plane_power(p_wave, GRID, n1), plane_power(s_wave, GRID, n1)
    expected = (ps * rp[0] + ss * rs[0]) / (ps + ss)
    assert reference_reflectance(mixed, GRID, n1, n2) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("distance", [0.012, -0.007, 0.0])
def test_axial_focus_finds_known_defocus(distance: float) -> None:
    """Pupil phase exp(−i k_z D) focuses D further along the propagation direction (mm)."""
    n = 1.0
    kx, ky = GRID.k_grid()
    kp = np.hypot(kx, ky)
    kz = np.sqrt(np.maximum((n * GRID.k0) ** 2 - kp**2, 0.0))
    spectrum = np.where(kp <= 0.3 * GRID.k0, np.exp(-1j * kz * distance), 0.0)
    shift = np.exp(-1j * (kx * (GRID.nx // 2) * GRID.dx + ky * (GRID.ny // 2) * GRID.dy))
    field = np.stack(
        [np.fft.ifft2(spectrum * shift), np.zeros((GRID.ny, GRID.nx), dtype=complex)], axis=-1
    )
    offsets = np.arange(-0.03, 0.03, 0.0002)
    best, peak = axial_focus(field, GRID, n, offsets)
    assert best == pytest.approx(distance, abs=2e-4)
    assert peak > 0
