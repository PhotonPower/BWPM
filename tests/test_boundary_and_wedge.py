"""K7 (absorbing mask, decision E2) and K8 (tilted glass wedge, decision E4) for the 2D FPWPM.

K7 reference: the same propagation in a window four times as wide without mask, where the beam
never reaches the border. K8 references: Snell's law and the Fresnel power transmittance at the
tilted surface (angles written out here, independent of bwpm.interface).
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm._types import ComplexArray
from bwpm.fpwpm import Grid2D, Propagator2D, absorbing_mask, power_flux
from bwpm.geometry import labels_area_fraction, labels_at_centres

LAMBDA_UM = 0.5
STEP = LAMBDA_UM * 1e-3 / 10  # λ/10 in mm


def _tilted_gaussian(
    ny: int, width_lambda: float, order: int, theta: float
) -> tuple[Grid2D, ComplexArray]:
    """Planes every 2.5 µm of a TE Gaussian (w0 = 2 µm) tilted by theta, propagated 40 µm in air."""
    grid = Grid2D(ny, STEP, STEP, LAMBDA_UM)
    nz = round(40e-3 / STEP)
    mask = absorbing_mask(grid, width_lambda * LAMBDA_UM * 1e-3, order) if width_lambda else None
    field = np.zeros((ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((grid.y / 2e-3) ** 2) + 1j * grid.k0 * np.sin(theta) * grid.y)
    planes = Propagator2D(grid, [1.0], np.zeros((nz, ny), dtype=np.int64), 1.0, mask).run(
        field, keep=True
    )
    return grid, np.asarray(planes)[::50, :, 0]


def test_k7_mask_removes_wraparound() -> None:
    theta = np.deg2rad(30.0)
    ny, big = 1024, 4096
    _, reference = _tilted_gaussian(big, 0, 8, theta)
    offset = (big - ny) // 2
    errors = {}
    for width in (0, 10, 15, 20, 30):
        grid, planes = _tilted_gaussian(ny, width, 8, theta)
        inner = np.abs(grid.y) <= np.max(grid.y) - width * LAMBDA_UM * 1e-3
        window = reference[:, offset : offset + ny]
        errors[width] = float(np.max(np.abs(planes[:, inner] - window[:, inner])))
    assert errors[0] > 0.1  # periodic wrap-around without mask
    assert errors[10] > errors[15] > errors[20] > errors[30]
    assert errors[20] < 1e-5  # rule: boundary zone ≥ 20 λ, order 8 (docs/testfaelle.md, K7)


def _wedge(
    alpha_deg: float, polarisation: int, rule: str = "centre"
) -> tuple[float, float, float, float]:
    """Glass wedge n = 1.5 with surface z = 10 µm + y·tan α, normal incidence of a Gaussian.

    Returns (direction in glass, expected direction, transmitted power, Fresnel transmittance),
    directions in degrees against z.
    """
    n = 1.5
    grid = Grid2D(2048, STEP, STEP, LAMBDA_UM)
    nz = round(40e-3 / STEP)
    alpha = np.deg2rad(alpha_deg)

    def inside(y: np.ndarray, z: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
        return np.asarray(z > 10e-3 + y * np.tan(alpha))

    if rule == "centre":
        indices, labels = [1.0 + 0j, n + 0j], labels_at_centres(grid, nz, 0.0, inside)
    else:
        indices, labels = labels_area_fraction(grid, nz, 0.0, inside, 1.0, n)
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, polarisation] = np.exp(-((grid.y / 5e-3) ** 2))
    mask = absorbing_mask(grid, 20 * LAMBDA_UM * 1e-3, 8)
    out = Propagator2D(grid, indices, labels, 1.0, mask).run(field)
    spectrum = np.abs(np.fft.fft(out[:, polarisation])) ** 2
    k_centre = float(np.sum(grid.ky * spectrum) / np.sum(spectrum))
    k0 = grid.k0
    k_t = k0 * np.sin(alpha)  # tangential component at the tilted surface
    k_n = np.sqrt((n * k0) ** 2 - k_t**2)
    ky_expected = k_t * np.cos(alpha) - k_n * np.sin(alpha)
    cos_i, cos_t = np.cos(alpha), k_n / (n * k0)
    if polarisation == 0:
        t = 2 * cos_i / (cos_i + n * cos_t)
    else:
        t = 2 * cos_i / (n * cos_i + cos_t)
    transmittance = n * cos_t / cos_i * abs(t) ** 2
    return (
        float(np.degrees(np.arcsin(k_centre / (n * k0)))),
        float(np.degrees(np.arcsin(ky_expected / (n * k0)))),
        power_flux(out, grid, n) / power_flux(field, grid, 1.0),
        float(transmittance),
    )


def test_k8_normal_surface_exact() -> None:
    # z-normal surface: exact up to the angular spread of the Gaussian (w0 = 10 λ): ΔT ≈ 1.4e-5
    _, _, power, fresnel = _wedge(0.0, 0)
    assert abs(power - fresnel) < 5e-5


@pytest.mark.parametrize("alpha", [10.0, 30.0, 45.0])
@pytest.mark.parametrize("polarisation", [0, 1], ids=["TE", "TM"])
def test_k8_wedge_direction_and_power(alpha: float, polarisation: int) -> None:
    direction, expected, power, fresnel = _wedge(alpha, polarisation)
    assert abs(direction - expected) < 0.02  # degrees: refraction direction is right
    # power: the staircase applies Fresnel at z-normal steps; measured deviation up to 7 %
    assert abs(power - fresnel) < 0.08


def test_k8_area_fraction_rule() -> None:
    direction, expected, power, fresnel = _wedge(30.0, 0, rule="area")
    assert abs(direction - expected) < 0.02
    assert abs(power - fresnel) < 0.08
