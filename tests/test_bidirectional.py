"""K4 (B2): BWPM for planar layers against the exact multilayer solution of the discretised stack.

Reference (independent of bwpm): tangential fields with forward/backward amplitudes in every layer,
continuity of the tangential E and H at every interface (H_t ∝ η·(E⁺ − E⁻) with η = k_z for TE and
η = ñ²k₀²/k_z for TM, E_t = E_x resp. E_y), propagation exp(±i k_z d). The discretised stack of
the BWPM: medium `before` up to plane 0, slab j between planes j and j + 1, the last slab's medium
continues behind plane nz. Plane waves e^(i k_y y) on the grid are exact eigenfunctions of every
step, so after enough iterations the BWPM must agree to rounding.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D

LAMBDA_UM = 1.0  # as wende2023 §2.1.3
STEP = LAMBDA_UM * 1e-3 / 20
GRID = Grid2D(ny=64, dy=STEP, dz=STEP, wavelength_um=LAMBDA_UM)
K0 = 2 * np.pi / (LAMBDA_UM * 1e-3)
SILVER = 0.213 + 6.697j


def _kz(n: complex, ky: float) -> complex:
    q = complex(np.sqrt(complex(n * n * K0 * K0 - ky * ky)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _multilayer(
    media: list[complex], before: complex, ky: float, pol: str
) -> tuple[complex, complex]:
    """(r, t) of the tangential field: E⁻ at plane 0 per E⁺ at plane 0; E⁺ at plane nz per E⁺(0).

    `media[j]` is the index of slab j (thickness dz); the last slab continues to +∞.
    """
    d = GRID.dz
    layers = [before, *media]
    kz = [_kz(n, ky) for n in layers]
    eta = [k if pol == "TE" else n * n * K0 * K0 / k for n, k in zip(layers, kz, strict=True)]
    # start in the last medium at plane nz - 1 (its left interface) with E+ = 1, E- = 0
    a_plus, a_minus = 1.0 + 0j, 0.0 + 0j
    for j in range(len(layers) - 1, 0, -1):
        # interface j - 1 at plane j - 1: left medium layers[j - 1], right layers[j]
        e_t = a_plus + a_minus
        h_t = eta[j] * (a_plus - a_minus)
        left_plus = 0.5 * (e_t + h_t / eta[j - 1])
        left_minus = 0.5 * (e_t - h_t / eta[j - 1])
        if j - 1 == 0:
            a_plus, a_minus = left_plus, left_minus
        else:
            # move to the left edge of slab j - 2's right neighbour: through slab (j - 1) of
            # thickness dz, i.e. from plane j - 1 back to plane j - 2
            a_plus = left_plus * np.exp(-1j * kz[j - 1] * d)
            a_minus = left_minus * np.exp(1j * kz[j - 1] * d)
    r = a_minus / a_plus
    t = np.exp(1j * kz[-1] * d) / a_plus  # last slab: from its left interface to plane nz
    return complex(r), complex(t)


def _run(
    media: list[complex], before: complex, m: int, pol: str, iterations: int
) -> tuple[complex, complex, list[float]]:
    labels_list = sorted(set(media), key=lambda c: (c.real, c.imag))
    labels = np.array([[labels_list.index(n)] * GRID.ny for n in media], dtype=np.int64)
    prop = Propagator2D(GRID, labels_list, labels, before)
    ky = 2 * np.pi * m / (GRID.ny * GRID.dy)
    field = np.zeros((GRID.ny, 2), dtype=np.complex128)
    col = 0 if pol == "TE" else 1
    field[:, col] = np.exp(1j * ky * GRID.y)
    res = Bidirectional2D(prop).run(field, iterations=iterations)
    r = complex(np.mean(res.reflected[:, col] / field[:, col]))
    t = complex(np.mean(res.transmitted[:, col] / field[:, col]))
    return r, t, res.energy_residual


def _slab(n: complex, thickness_steps: int) -> list[complex]:
    return [1.0 + 0j] * 5 + [complex(n)] * thickness_steps + [1.0 + 0j] * 5


STACKS = {
    # wende2023 §2.1.3: slab 2.5 µm, n = 1.5 in vacuum, λ = 1 µm
    "slab": _slab(1.5, 50),
    "multilayer": [1.0 + 0j] * 3
    + [1.45 + 0j] * 7
    + [2.1 + 0j] * 4
    + [1.45 + 0j] * 9
    + [1.52 + 0j] * 6,
    "silver": [1.0 + 0j] * 4 + [1.505 + 0j] * 10 + [SILVER] * 8 + [1.505 + 0j] * 5,
}


@pytest.mark.parametrize("name", list(STACKS))
@pytest.mark.parametrize(
    "m", [0, 1, 2, 3, 5]
)  # sin θ = 0, 0.31, 0.63, 0.94 in vacuum; 5: evanescent
@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_k4_converges_to_exact_multilayer(name: str, m: int, pol: str) -> None:
    if name == "multilayer" and m == 5:
        return  # trapped wave, see test_k4_trapped_wave_does_not_converge
    media = STACKS[name]
    ky = 2 * np.pi * m / (GRID.ny * GRID.dy)
    r_exact, t_exact = _multilayer(media, 1.0, ky, pol)
    r, t, _ = _run(media, 1.0, m, pol, iterations=60)
    assert abs(r - r_exact) < 1e-12
    assert abs(t - t_exact) < 1e-12


def test_k4_trapped_wave_does_not_converge() -> None:
    # sin θ = 1.56: evanescent in vacuum, n = 1.45 and n = 1.52, propagating in n = 2.1 between
    # them. The iteration is a Neumann series of multiple reflections; for a wave trapped by total
    # internal reflection its spectral radius is not below 1, and it does not converge (limit of
    # the method, docs/testfaelle.md, K4).
    media = STACKS["multilayer"]
    ky = 2 * np.pi * 5 / (GRID.ny * GRID.dy)
    r_exact, _ = _multilayer(media, 1.0, ky, "TE")
    errors = [abs(_run(media, 1.0, 5, "TE", iterations=n)[0] - r_exact) for n in (60, 200)]
    assert min(errors) > 1e-3


def test_k4_slab_against_airy() -> None:
    # ideal slab (thickness 50 dz) in vacuum, TE: Airy formula with Fresnel coefficients
    n, thickness = 1.5, 50 * GRID.dz
    for m in (0, 1, 2, 3):
        ky = 2 * np.pi * m / (GRID.ny * GRID.dy)
        k1, k2 = _kz(1.0, ky), _kz(n, ky)
        r12 = (k1 - k2) / (k1 + k2)
        t12, t21 = 2 * k1 / (k1 + k2), 2 * k2 / (k1 + k2)
        phase = np.exp(2j * k2 * thickness)
        r_airy = (r12 - r12 * phase) / (1 - r12 * r12 * phase)
        r, _, _ = _run(STACKS["slab"], 1.0, m, "TE", iterations=60)
        # the slab starts at plane 5: the reflection is referred to plane 0 (5 dz of vacuum, twice)
        assert abs(r - r_airy * np.exp(2j * k1 * 5 * GRID.dz)) < 1e-12
        assert t12 * t21 != 0


def test_k4_energy_residual_decreases() -> None:
    # wende2023 Fig. 2(e, f): ε = 1 − |R + T| falls roughly geometrically with the iterations
    for m in (1, 2):
        for pol in ("TE", "TM"):
            _, _, eps = _run(STACKS["slab"], 1.0, m, pol, iterations=12)
            assert abs(eps[0]) > 1e-3
            assert all(abs(b) < abs(a) for a, b in pairwise(eps[:6]))
            assert abs(eps[-1]) < 1e-12
