"""K19: back scattering at weak fluctuations δn in the BWPM (B8a, E12, docs/spezifikation.md §11).

References independent of bwpm: the exact multilayer solution (transfer matrix of tangential E and
H, as in test_bidirectional.py) of the stack ñ₀ + δn_j, the closed slab integral of the first
Born approximation, and the RCWA of a lamellar grating in reflection (`tests/rcwa_reference.py`,
itself checked here against Airy and energy conservation). The structure is a homogeneous host.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pytest
from rcwa_reference import lamellar_reflection, lamellar_transmission

from bwpm.bidirectional import Bidirectional2D
from bwpm.fluctuation import born_emission, born_sources
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.volume3d import Bidirectional3D, Grid3D

N0 = 1.5


def _kz(n: complex, k0: float, ky: float) -> complex:
    q = complex(np.sqrt(complex(n * n * k0 * k0 - ky * ky)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _stack(
    media: list[complex], before: complex, k0: float, dz: float, ky: float, pol: str
) -> tuple[complex, complex]:
    """Exact (r, t) of the discretised stack: slab j of thickness dz has index media[j], the last
    medium continues; r refers to plane 0, t to plane nz (tangential E, TE: E_x, TM: E_y)."""
    layers = [before, *media]
    kz = [_kz(n, k0, ky) for n in layers]
    eta = [k if pol == "TE" else n * n * k0 * k0 / k for n, k in zip(layers, kz, strict=True)]
    a_plus, a_minus = 1.0 + 0j, 0.0 + 0j
    for j in range(len(layers) - 1, 0, -1):
        e_t = a_plus + a_minus
        h_t = eta[j] * (a_plus - a_minus)
        left_plus = 0.5 * (e_t + h_t / eta[j - 1])
        left_minus = 0.5 * (e_t - h_t / eta[j - 1])
        if j - 1 == 0:
            a_plus, a_minus = left_plus, left_minus
        else:
            a_plus = left_plus * np.exp(-1j * kz[j - 1] * dz)
            a_minus = left_minus * np.exp(1j * kz[j - 1] * dz)
    return complex(a_minus / a_plus), complex(np.exp(1j * kz[-1] * dz) / a_plus)


def _bwpm(
    grid: Grid2D,
    delta: np.ndarray,  # type: ignore[type-arg]
    m: int,
    pol: str,
    iterations: int,
) -> tuple[complex, complex, list[float]]:
    labels = np.zeros((delta.shape[0], grid.ny), dtype=np.int64)
    fluctuation = np.repeat(delta[:, None], grid.ny, axis=1)
    prop = Propagator2D(grid, [complex(N0)], labels, complex(N0), fluctuation=fluctuation)
    ky = 2 * np.pi * m / (grid.ny * grid.dy)
    col = 0 if pol == "TE" else 1
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, col] = np.exp(1j * ky * grid.y)
    res = Bidirectional2D(prop).run(field, iterations=iterations)
    r = complex(np.mean(res.reflected[:, col] / field[:, col]))
    t = complex(np.mean(res.transmitted[:, col] / field[:, col]))
    return r, t, res.energy_residual


def test_single_slab_back_scattering_matches_the_exact_born_integral() -> None:
    """R1 of one slab with uniform δn for a plane wave (normal incidence, TE): the first Born
    amplitude at the left edge is i V/(2k)·∫₀^dz e^(2ikζ) dζ; with two Gauss points the slab
    integral has the error of the quadrature, < 2·10⁻⁴ for 2k·dz = 0.94 (docs/b8-entwurf.md)."""
    grid = Grid2D(ny=16, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.0)
    k = N0 * grid.k0
    delta_n = 1e-4
    spectrum = np.zeros((grid.ny, 2), dtype=np.complex128)
    spectrum[0, 0] = grid.ny  # E_x = 1 everywhere
    kx = np.zeros_like(grid.ky)
    delta = np.full(grid.ny, delta_n, dtype=np.complex128)
    sources = born_sources(spectrum, +1, delta, N0, grid.k0, grid.dz, kx, grid.ky)
    back = born_emission(sources, -1, N0, grid.k0, grid.dz, kx, grid.ky)
    # emission is referred to the right edge: back to the left edge by e^(+ik dz)
    at_left = complex(np.mean(back[:, 0])) * np.exp(1j * k * grid.dz)
    potential = 2 * N0 * grid.k0**2 * delta_n
    exact = 1j * potential / (2 * k) * (np.exp(2j * k * grid.dz) - 1) / (2j * k)
    assert abs(at_left / exact - 1) < 2e-4
    forward = born_emission(sources, +1, N0, grid.k0, grid.dz, kx, grid.ky)
    exact_fwd = 1j * potential / (2 * k) * grid.dz * np.exp(1j * k * grid.dz)
    assert abs(complex(np.mean(forward[:, 0])) / exact_fwd - 1) < 1e-12  # no phase mismatch


@pytest.mark.parametrize("pol", ["TE", "TM"])
@pytest.mark.parametrize("m", [0, 1, 2])
def test_slab_reflection_and_transmission_match_the_exact_stack(m: int, pol: str) -> None:
    """Slab δn = 10⁻³, 2.5 µm (50 dz, λ = 1 µm, ñ₀ = 1.5) at 0°, 12° and 25° in the host:
    converged BWPM against the exact stack. Remaining error of first order relative to the
    reflection: second order in δn of the Fresnel coefficients, δn/ñ₀ ≈ 7·10⁻⁴, plus the
    quadrature (2·10⁻⁴); tolerance 2·10⁻³ relative for r, 10⁻⁴ absolute for t."""
    grid = Grid2D(ny=64, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.0)
    delta = np.zeros(60)
    delta[5:55] = 1e-3
    r, t, _ = _bwpm(grid, delta, m, pol, iterations=20)
    ky = 2 * np.pi * m / (grid.ny * grid.dy)
    media = [complex(N0 + d) for d in delta]
    r_exact, t_exact = _stack(media, complex(N0), grid.k0, grid.dz, ky, pol)
    assert abs(r_exact) > 1e-4  # the reflection is resolved, not rounding
    assert abs(r / r_exact - 1) < 2e-3
    assert abs(t - t_exact) < 1e-4


@pytest.mark.slow
@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_bragg_grating_matches_the_exact_stack(pol: str) -> None:
    """Bragg grating in z: δn = 0.01·sin(2πz/Λ), Λ = 8 dz = 0.4 µm, 40 periods, λ = 2ñ₀Λ = 1.2 µm,
    normal incidence: strong coherent reflection (multiple scattering between slabs) against the
    exact stack; tolerance 1 % relative (second order in δn per slab, accumulated)."""
    grid = Grid2D(ny=8, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.2)
    z = (np.arange(320) + 0.5) * grid.dz
    delta = np.concatenate([np.zeros(4), 0.01 * np.sin(2 * np.pi * z / (8 * grid.dz)), np.zeros(4)])
    r, t, _ = _bwpm(grid, delta, 0, pol, iterations=40)
    media = [complex(N0 + d) for d in delta]
    r_exact, t_exact = _stack(media, complex(N0), grid.k0, grid.dz, 0.0, pol)
    assert abs(r_exact) ** 2 > 0.05
    assert abs(r / r_exact - 1) < 1e-2
    assert abs(t / t_exact - 1) < 1e-2


@pytest.mark.slow
def test_energy_is_conserved_for_lossless_fluctuations() -> None:
    """Bragg grating as above: 1 − (P_R + P_T)/P_in after convergence below 10⁻³ (growth of the
    norm of order (k₀δn dz)² per slab and pass, §10)."""
    grid = Grid2D(ny=8, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.2)
    z = (np.arange(320) + 0.5) * grid.dz
    delta = np.concatenate([np.zeros(4), 0.01 * np.sin(2 * np.pi * z / (8 * grid.dz)), np.zeros(4)])
    _, _, residual = _bwpm(grid, delta, 0, "TE", iterations=40)
    assert abs(residual[-1]) < 1e-3


def test_one_iteration_transmits_like_the_forward_propagator() -> None:
    """Gap of B7: Bidirectional2D ignored δn. The first forward pass must equal Propagator2D."""
    rng = np.random.default_rng(19)
    grid = Grid2D(ny=32, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.0)
    labels = np.zeros((20, grid.ny), dtype=np.int64)
    labels[12:] = 1
    fluctuation = np.zeros(labels.shape)
    fluctuation[2:10] = 0.01 * rng.standard_normal((8, grid.ny))
    prop = Propagator2D(grid, [complex(N0), 1.6 + 0j], labels, complex(N0), fluctuation=fluctuation)
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((grid.y / (6 * grid.dy)) ** 2))
    forward = prop.run(field)
    res = Bidirectional2D(prop).run(field, iterations=1)
    np.testing.assert_allclose(res.transmitted, forward, rtol=0, atol=1e-14)


def test_zero_fluctuation_is_bit_identical_in_bwpm_2d() -> None:
    rng = np.random.default_rng(20)
    grid = Grid2D(ny=32, dy=0.05e-3, dz=0.05e-3, wavelength_um=1.0)
    labels = rng.integers(0, 2, size=(15, grid.ny))
    indices = [complex(N0), 1.6 + 0.01j]
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    field[:, 1] = 1.0
    plain = Bidirectional2D(Propagator2D(grid, indices, labels, complex(N0))).run(field, 3)
    zero = Bidirectional2D(
        Propagator2D(grid, indices, labels, complex(N0), fluctuation=np.zeros(labels.shape))
    ).run(field, 3)
    assert np.array_equal(plain.reflected, zero.reflected)
    assert np.array_equal(plain.transmitted, zero.transmitted)


@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_bwpm_3d_reduces_to_2d(pol: Literal["TE", "TM"]) -> None:
    """Random δn varying along y and z only: 3D BWPM (3 iterations) equals the 2D BWPM."""
    rng = np.random.default_rng(21)
    ny, nx, nz, step = 32, 4, 24, 0.05e-3
    delta_2d = np.zeros((nz, ny))
    delta_2d[3:20] = 0.01 * rng.standard_normal((17, ny))
    labels_2d = np.zeros((nz, ny), dtype=np.int64)
    labels_2d[21:] = 1
    indices = [complex(N0), 1.6 + 0j]
    col = 0 if pol == "TE" else 1
    field_2d = np.zeros((ny, 2), dtype=np.complex128)
    field_2d[:, col] = 1.0
    grid_2d = Grid2D(ny=ny, dy=step, dz=step, wavelength_um=1.0)
    ref = Bidirectional2D(
        Propagator2D(grid_2d, indices, labels_2d, complex(N0), fluctuation=delta_2d)
    ).run(field_2d, 3)
    grid = Grid3D(nx, ny, step, step, step, 1.0)
    field = np.repeat(field_2d[:, None, :], nx, axis=1)
    res = Bidirectional3D(
        grid,
        indices,
        np.repeat(labels_2d[:, :, None], nx, axis=2),
        complex(N0),
        fluctuation=np.repeat(delta_2d[:, :, None], nx, axis=2),
    ).run(field, iterations=3)
    for column in range(nx):
        np.testing.assert_allclose(res.reflected[:, column], ref.reflected, rtol=0, atol=1e-12)
        np.testing.assert_allclose(res.transmitted[:, column], ref.transmitted, rtol=0, atol=1e-12)


@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_rcwa_reflection_is_exact_for_a_slab_and_conserves_energy(
    pol: Literal["TE", "TM"],
) -> None:
    """The reference itself: homogeneous layer 1.40 in 1.38 (2.1 µm, normal incidence) against Airy,
    and a strong lossless grating (1.38/1.6) conserving the flux over all orders (TM with the
    factor −k_z/k of the reflected E_y, docs/testfaelle.md K19)."""
    lam, n0, n1, thick = 0.8e-3, 1.38, 1.40, 2.1e-3
    k0 = 2 * np.pi / lam
    r = lamellar_reflection(n0, n1, n1, 0.8e-3, 0.0, 0.4e-3, thick, lam, 5, pol)[0]
    r01 = (n0 - n1) / (n0 + n1)
    phase = np.exp(2j * n1 * k0 * thick)
    assert abs(r - r01 * (1 - phase) / (1 - r01**2 * phase)) < 1e-14
    period = 0.8e-3
    refl = lamellar_reflection(n0, n0, 1.6, period, -0.01e-3, 0.39e-3, thick, lam, 15, pol)
    trans = lamellar_transmission(n0, n0, 1.6, period, -0.01e-3, 0.39e-3, thick, lam, 15, pol)
    k = n0 * k0
    total = 0.0
    for m in refl:
        ky = 2 * np.pi * m / period
        if ky * ky < k * k:
            kz = np.sqrt(k * k - ky * ky)
            weight = kz / k if pol == "TE" else k / kz
            total += (abs(refl[m]) ** 2 + abs(trans[m]) ** 2) * weight
    assert abs(total - 1) < 1e-12


@pytest.mark.parametrize("pol", ["TE", "TM"])
@pytest.mark.parametrize("delta_n", [0.01, 0.03])
def test_grating_reflection_matches_rcwa(delta_n: float, pol: Literal["TE", "TM"]) -> None:
    """Lamellar δn grating (period 0.8 µm, 2.1 µm thick, λ = 0.8 µm, ñ₀ = 1.38, dy = dz = 0.02 µm):
    reflected orders 0 and ±1 (±1 at 46° backwards) after 10 iterations against the RCWA. The
    remaining error is of second order in δn, of the size δn/ñ₀ (§10, §11); tolerance δn."""
    period_um, px, step_um, thick_um = 0.8, 40, 0.02, 2.1
    grid = Grid2D(ny=px, dy=period_um / px * 1e-3, dz=step_um * 1e-3, wavelength_um=0.8)
    nz = round(thick_um / step_um)
    labels = np.zeros((nz + 4, px), dtype=np.int64)
    fluctuation = np.zeros((nz + 4, px))
    fluctuation[:nz, : px // 2] = delta_n
    prop = Propagator2D(grid, [1.38 + 0j], labels, 1.38 + 0j, fluctuation=fluctuation)
    col = 0 if pol == "TE" else 1
    field = np.zeros((px, 2), dtype=np.complex128)
    field[:, col] = 1.0
    reflected = Bidirectional2D(prop).run(field, iterations=10).reflected
    spectrum = np.fft.fft(reflected[:, col]) / px
    dy = grid.dy
    # pixel n is the cell |y − n·dy| ≤ dy/2 in the sample order of the FFT (as K17, K18)
    rcwa = lamellar_reflection(
        1.38, 1.38, 1.38 + delta_n, period_um * 1e-3, -dy / 2, (px // 2 - 0.5) * dy,
        thick_um * 1e-3, 0.8e-3, 15, pol,
    )  # fmt: skip
    for m in (0, 1, -1):
        assert abs(rcwa[m]) > 1e-4
        assert abs(spectrum[m] / rcwa[m] - 1) < delta_n, (m, spectrum[m], rcwa[m])
