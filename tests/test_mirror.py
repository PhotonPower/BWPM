"""K5 and K9 (B2): reflecting curved and tilted surfaces with the 2D BWPM.

K5: silver parabolic mirror z = z_v − y²/(4f), f = 50 µm, collimated TE Gaussian (wende2023
§2.1.4 as plausibility): the on-axis intensity maximum of the backward field lies at z_v − f; a
parabola focuses a collimated on-axis beam without aberration.
K9: tilted plane silver mirror: the reflected beam leaves at 2α (law of reflection); the power is
recorded against Fresnel (a limit of the method, docs/testfaelle.md, K9).
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D, absorbing_mask, power_flux
from bwpm.geometry import labels_at_centres

SILVER = 0.213 + 6.697j


@pytest.mark.slow
def test_k5_parabolic_mirror_focus() -> None:
    lam, n_max = 0.78, 1.505
    d = lam * 1e-3 / (10 * n_max)  # wende2023 §2.1.4: λ/(10 n_max)
    f, z_v = 50e-3, 70e-3
    ny = round(160e-3 / d)
    ny += ny % 2
    grid = Grid2D(ny, d, d, lam)
    nz = round((z_v + 6e-3) / d)
    z = (np.arange(nz) + 0.5) * d
    zz, yy = np.meshgrid(z, grid.y, indexing="ij")
    surface = z_v - yy**2 / (4 * f)
    labels = np.zeros((nz, ny), dtype=np.int64)
    labels[zz >= surface] = 1  # silver 0.4 µm
    labels[zz >= surface + 0.4e-3] = 2  # IP-S substrate
    field = np.zeros((ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((grid.y / 30e-3) ** 2))
    prop = Propagator2D(
        grid, [1.0, SILVER, 1.505], labels, 1.0, absorbing_mask(grid, 20 * lam * 1e-3, 8)
    )
    res = Bidirectional2D(prop).run(field, iterations=1, keep=True)
    assert res.backward is not None
    intensity = np.abs(res.backward[:, ny // 2, 0]) ** 2
    planes = np.arange(nz + 1) * d
    intensity[planes > z_v - 5e-3] = 0.0
    k = int(np.argmax(intensity))
    a, b, c = intensity[k - 1 : k + 2]
    z_focus = planes[k] + 0.5 * (a - c) / (a - 2 * b + c) * d
    # measured 49.87 µm from the vertex (λ/10 and λ/20 agree to 0.01 µm); wende2023: 49.85 µm
    assert abs((z_v - z_focus) - f) < 0.3e-3


def _tilted_mirror(alpha_deg: float, fraction: int) -> tuple[float, float, float]:
    """(reflected direction in degrees, reflected power, Fresnel |r|²) for a tilted silver plane."""
    lam = 0.5
    d = lam * 1e-3 / fraction
    grid = Grid2D(round(102.4e-3 / d), d, d, lam)
    nz = round(30e-3 / d)
    alpha = np.deg2rad(alpha_deg)

    def inside(y: np.ndarray, z: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
        return np.asarray(z > 15e-3 + y * np.tan(alpha))

    labels = labels_at_centres(grid, nz, 0.0, inside)
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((grid.y / 8e-3) ** 2))
    mask = absorbing_mask(grid, 20 * lam * 1e-3, 8)
    res = Bidirectional2D(Propagator2D(grid, [1.0, SILVER], labels, 1.0, mask)).run(field, 1)
    spectrum = np.abs(np.fft.fft(res.reflected[:, 0])) ** 2
    k_centre = float(np.sum(grid.ky * spectrum) / np.sum(spectrum))
    power = power_flux(res.reflected, grid, 1.0, direction=-1) / power_flux(field, grid, 1.0)
    kz = np.sqrt(SILVER**2 - np.sin(alpha) ** 2)
    r = (np.cos(alpha) - kz) / (np.cos(alpha) + kz)
    return float(abs(np.degrees(np.arcsin(k_centre / grid.k0)))), power, float(abs(r) ** 2)


@pytest.mark.slow
@pytest.mark.parametrize("alpha", [5.0, 10.0, 20.0])
def test_k9_tilted_mirror_direction_converges(alpha: float) -> None:
    coarse, _, _ = _tilted_mirror(alpha, 10)
    fine, power, fresnel = _tilted_mirror(alpha, 20)
    assert abs(fine - 2 * alpha) < 0.05  # λ/40 gives < 0.01° (docs/testfaelle.md, K9)
    assert abs(fine - 2 * alpha) <= abs(coarse - 2 * alpha)
    # limit of the method: the reflected power falls short of Fresnel and does not converge
    assert power < fresnel
