"""K10–K12 (B3): BWPM in 3D.

K10 homogeneous medium, K11 planar stack against the transfer matrix for arbitrary azimuth and
polarisation (s and p separately with tangential fields, independent of bwpm), K12 consistency
with 2D and of sparse against full storage of the reflections.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm._types import ComplexArray
from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d, memory_estimate

LAMBDA_UM = 1.0
STEP = LAMBDA_UM * 1e-3 / 16
GRID = Grid3D(nx=32, ny=32, dx=STEP, dy=STEP, dz=STEP, wavelength_um=LAMBDA_UM)
K0 = 2 * np.pi / (LAMBDA_UM * 1e-3)
SILVER = 0.213 + 6.697j


def _kz(n: complex, kp2: float) -> complex:
    q = complex(np.sqrt(complex(n * n * K0 * K0 - kp2)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _plane_wave(
    mx: int, my: int, e_perp: tuple[complex, complex]
) -> tuple[ComplexArray, float, float]:
    kx = 2 * np.pi * mx / (GRID.nx * GRID.dx)
    ky = 2 * np.pi * my / (GRID.ny * GRID.dy)
    phase = np.exp(1j * (kx * GRID.x[None, :] + ky * GRID.y[:, None]))
    field = np.stack([e_perp[0] * phase, e_perp[1] * phase], axis=-1)
    return np.asarray(field), kx, ky


def _stack(media: list[complex]) -> tuple[list[complex], np.ndarray]:  # type: ignore[type-arg]
    unique = sorted(set(media), key=lambda c: (c.real, c.imag))
    labels = np.array([np.full((GRID.ny, GRID.nx), unique.index(n)) for n in media])
    return unique, labels


@pytest.mark.parametrize(("mx", "my"), [(0, 0), (3, 1), (-2, 4)])
def test_k10_homogeneous_plane_wave(mx: int, my: int) -> None:
    n, nz = 1.4 + 0.01j, 20
    indices, labels = _stack([n] * nz)
    field, kx, ky = _plane_wave(mx, my, (0.3 + 0.1j, -0.7))
    out = Bidirectional3D(GRID, indices, labels, n).run(field, iterations=0).transmitted
    expected = field * np.exp(1j * _kz(n, kx * kx + ky * ky) * nz * GRID.dz)
    assert np.max(np.abs(out - expected)) < 1e-12


def _transfer(
    media: list[complex], before: complex, kp2: float, pol: str
) -> tuple[complex, complex]:
    """(r, t) of the tangential field component (same scheme as tests/test_bidirectional.py)."""
    layers = [before, *media]
    kz = [_kz(n, kp2) for n in layers]
    eta = [k if pol == "s" else n * n * K0 * K0 / k for n, k in zip(layers, kz, strict=True)]
    a_plus, a_minus = 1.0 + 0j, 0.0 + 0j
    for j in range(len(layers) - 1, 0, -1):
        e_t, h_t = a_plus + a_minus, eta[j] * (a_plus - a_minus)
        left_plus, left_minus = 0.5 * (e_t + h_t / eta[j - 1]), 0.5 * (e_t - h_t / eta[j - 1])
        if j - 1 == 0:
            a_plus, a_minus = left_plus, left_minus
        else:
            a_plus = left_plus * np.exp(-1j * kz[j - 1] * GRID.dz)
            a_minus = left_minus * np.exp(1j * kz[j - 1] * GRID.dz)
    return complex(a_minus / a_plus), complex(np.exp(1j * kz[-1] * GRID.dz) / a_plus)


# |k⊥| = k0·sqrt(mx² + my²)/2: 30° along x, 45° at azimuth 45° and 135°, normal incidence
# (evanescent directions in air that are guided in n = 2.1 would be the trapped case of K4)
@pytest.mark.parametrize(("mx", "my"), [(1, 0), (1, 1), (-1, 1), (0, 0)])
def test_k11_planar_stack_vector_coupling(mx: int, my: int) -> None:
    # 30 iterations: converged to rounding. Many more iterations let rounding noise in waves that
    # are trapped in the n = 2.1 layer grow by about 1.3 per iteration (docs/testfaelle.md, K11).
    media = [1.0 + 0j] * 3 + [1.5 + 0j] * 10 + [2.1 + 0j] * 5 + [SILVER] * 3 + [1.5 + 0j] * 4
    indices, labels = _stack(media)
    e0 = (0.6 + 0.2j, -0.3 + 0.5j)
    field, kx, ky = _plane_wave(mx, my, e0)
    res = Bidirectional3D(GRID, indices, labels, 1.0 + 0j).run(field, iterations=30)
    kp2 = kx * kx + ky * ky
    # tangential basis: s along (−k_y, k_x)/k⊥, p (tangential part) along (k_x, k_y)/k⊥
    if kp2 > 0:
        s_dir = np.array([-ky, kx]) / np.sqrt(kp2)
        p_dir = np.array([kx, ky]) / np.sqrt(kp2)
    else:
        s_dir, p_dir = np.array([0.0, 1.0]), np.array([1.0, 0.0])
    e = np.array(e0)
    r_s, t_s = _transfer(media, 1.0, kp2, "s")
    r_p, t_p = _transfer(media, 1.0, kp2, "p")
    a_s, a_p = e @ s_dir, e @ p_dir
    expected_t = t_s * a_s * s_dir + t_p * a_p * p_dir
    expected_r = r_s * a_s * s_dir + r_p * a_p * p_dir
    centre = (GRID.ny // 2, GRID.nx // 2)  # y = x = 0, phase of the plane wave 1
    assert np.max(np.abs(res.transmitted[centre] - expected_t)) < 1e-12
    assert np.max(np.abs(res.reflected[centre] - expected_r)) < 1e-12


def test_k12_equals_2d_and_sparse_equals_full() -> None:
    # a volume homogeneous in x with a field constant in x is the 2D problem
    nz = 40
    g2 = Grid2D(ny=GRID.ny, dy=GRID.dy, dz=GRID.dz, wavelength_um=LAMBDA_UM)
    z = (np.arange(nz) + 0.5) * GRID.dz
    zz, yy = np.meshgrid(z, g2.y, indexing="ij")
    labels2 = np.where(zz > 10 * GRID.dz + 0.4 * yy, 1, 0).astype(np.int64)
    labels2[zz > 30 * GRID.dz] = 2
    indices = [1.0 + 0j, 1.5 + 0j, SILVER]
    field2 = np.zeros((GRID.ny, 2), dtype=np.complex128)
    field2[:, 0] = np.exp(-((g2.y / (8 * STEP)) ** 2))
    field2[:, 1] = 0.5j * np.exp(-((g2.y / (8 * STEP)) ** 2))
    res2 = Bidirectional2D(Propagator2D(g2, indices, labels2, 1.0)).run(field2, iterations=2)
    labels3 = np.repeat(labels2[:, :, None], GRID.nx, axis=2)
    field3 = np.repeat(field2[:, None, :], GRID.nx, axis=1)
    sparse = Bidirectional3D(GRID, indices, labels3, 1.0).run(field3, iterations=2)
    full = Bidirectional3D(GRID, indices, labels3, 1.0, sparse=False).run(field3, iterations=2)
    for m in (0, 7, GRID.nx - 1):
        assert np.max(np.abs(sparse.transmitted[:, m] - res2.transmitted)) < 1e-12
        assert np.max(np.abs(sparse.reflected[:, m] - res2.reflected)) < 1e-12
    assert np.array_equal(sparse.transmitted, full.transmitted)
    assert np.array_equal(sparse.reflected, full.reflected)
    assert 0 < sparse.stored_values < full.stored_values / 5


def test_mask_and_memory_model() -> None:
    m = absorbing_mask_2d(GRID, 4 * STEP)
    assert m[GRID.ny // 2, GRID.nx // 2] == 1.0 and m[0, 0] < 1e-6
    assert (
        memory_estimate(1000, 1000, 2_000_000, 5)
        == 4 * 5 * 10**6 * 64 + 8 * 10**6 * 32 + 2 * 2_000_000 * 40
    )
