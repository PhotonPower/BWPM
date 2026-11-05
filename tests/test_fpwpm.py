"""K2, K3 (B1): forward FPWPM in 2D against plane-wave solutions and the angular spectrum.

References are independent of bwpm.fpwpm and bwpm.interface: plane waves on the grid with
k_z from the dispersion relation, Fresnel coefficients written with angles (Byrnes eq. (6)),
and for TM the transverse ratio E_y,t/E_y,i = t_p·cos θ₂/cos θ₁.
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm._types import ComplexArray
from bwpm.fpwpm import Grid2D, Propagator2D, power_flux

LAMBDA_UM = 0.5
GRID = Grid2D(ny=256, dy=0.05e-3, dz=0.05e-3, wavelength_um=LAMBDA_UM)  # λ/10 sampling
K0 = 2 * np.pi / (LAMBDA_UM * 1e-3)
SILVER = 0.213 + 6.697j


def _kz(n: complex, ky: float) -> complex:
    q = complex(np.sqrt(complex(n * n * K0 * K0 - ky * ky)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _plane_wave(m: int, polarisation: str) -> tuple[ComplexArray, float]:
    """Plane wave e^(i k_y y) with k_y = 2π m/(ny dy) on the grid, unit transverse component."""
    ky = 2 * np.pi * m / (GRID.ny * GRID.dy)
    field = np.zeros((GRID.ny, 2), dtype=np.complex128)
    field[:, 0 if polarisation == "TE" else 1] = np.exp(1j * ky * GRID.y)
    return field, ky


@pytest.mark.parametrize("index", [1.0, 1.5, 1.5 + 0.01j])
@pytest.mark.parametrize("m", [0, 20, 90])
def test_k2_homogeneous_plane_wave_exact(index: complex, m: int) -> None:
    nz = 40
    prop = Propagator2D(GRID, [index], np.zeros((nz, GRID.ny), dtype=np.int64), index)
    for pol in ("TE", "TM"):
        field, ky = _plane_wave(m, pol)
        out = prop.run(field)
        expected = field * np.exp(1j * _kz(index, ky) * nz * GRID.dz)
        assert np.max(np.abs(out - expected)) < 1e-12


def test_k2_gaussian_equals_single_angular_spectrum_step() -> None:
    # many small steps = one step over the full length (homogeneous medium, no mask)
    nz, index = 200, 1.3
    w0 = 1.5e-3
    field = np.zeros((GRID.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((GRID.y / w0) ** 2))
    field[:, 1] = 0.3j * np.exp(-((GRID.y / w0) ** 2))
    out = Propagator2D(GRID, [index], np.zeros((nz, GRID.ny), dtype=np.int64), index).run(field)
    kz = np.array([_kz(index, k) for k in GRID.ky])
    reference = np.fft.ifft(
        np.fft.fft(field, axis=0) * np.exp(1j * kz * nz * GRID.dz)[:, None], axis=0
    )
    assert np.max(np.abs(out - reference)) < 1e-12 * np.max(np.abs(field))


def _fresnel_transverse(n1: complex, n2: complex, ky: float) -> tuple[complex, complex]:
    """Transverse transmission factors (TE: t_s, TM: t_p cos θ₂/cos θ₁), angles from k_z."""
    cos1, cos2 = _kz(n1, ky) / (n1 * K0), _kz(n2, ky) / (n2 * K0)
    t_s = 2 * n1 * cos1 / (n1 * cos1 + n2 * cos2)
    t_p = 2 * n1 * cos1 / (n2 * cos1 + n1 * cos2)
    return t_s, t_p * cos2 / cos1


@pytest.mark.parametrize(
    ("n1", "n2"), [(1.0, 1.5), (1.5, 1.0), (1.0, 1.5 + 0.02j), (1.505, 1.511), (1.0, SILVER)]
)
@pytest.mark.parametrize("m", [0, 25, 70])
def test_k3_plane_interface_exact(n1: complex, n2: complex, m: int) -> None:
    # interface after slab 9 (at z = 10 dz), propagate to z = 30 dz
    nz, n_before = 30, 10
    labels = np.zeros((nz, GRID.ny), dtype=np.int64)
    labels[n_before:] = 1
    prop = Propagator2D(GRID, [n1, n2], labels, n1)
    for pol, column in (("TE", 0), ("TM", 1)):
        field, ky = _plane_wave(m, pol)
        out = prop.run(field)
        t_te, t_tm = _fresnel_transverse(n1, n2, ky)
        phase = np.exp(1j * (_kz(n1, ky) * n_before + _kz(n2, ky) * (nz - n_before)) * GRID.dz)
        expected = (t_te if pol == "TE" else t_tm) * phase * field[:, column]
        assert np.max(np.abs(out[:, column] - expected)) < 1e-12 * max(1.0, abs(t_te))
        assert np.max(np.abs(out[:, 1 - column])) < 1e-13


def test_k3_power_transmission_of_a_gaussian() -> None:
    # TE Gaussian air -> glass: transmitted power = Σ_k T_s(k)·P_in(k) (angular spectrum)
    n1, n2, w0 = 1.0, 1.5, 1.0e-3
    nz = 20
    labels = np.ones((nz, GRID.ny), dtype=np.int64)
    field = np.zeros((GRID.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((GRID.y / w0) ** 2))
    out = Propagator2D(GRID, [n1, n2], labels, n1).run(field)
    spectrum = np.abs(np.fft.fft(field[:, 0])) ** 2
    weights = []
    for k, power in zip(GRID.ky, spectrum, strict=True):
        kz1, kz2 = _kz(n1, float(k)), _kz(n2, float(k))
        t_s = 2 * kz1 / (kz1 + kz2)
        weights.append(power * abs(t_s) ** 2 * kz2.real)
    expected = float(np.sum(weights)) * GRID.dy / GRID.ny
    assert abs(power_flux(out, GRID, n2) / expected - 1.0) < 1e-12
