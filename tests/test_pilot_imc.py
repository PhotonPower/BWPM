"""Fusion pilot K6 (2D): fibre facet with index-matching cap, Field Tracing → BWPM → Field Tracing.

Needs the optional package fieldtracing (extra `fusion`); skipped otherwise. References are
independent of the chain:
- reflected focus: a plane mirror images the aberration-free air focus of the objective, so the
  light reflected at the IMC front focuses at front − x (x = depth of the air focus behind it);
- powers: Fresnel at the facet, bounded by core and cladding indices, weighted with the
  aplanatic spectrum (amplitude ∝ cos^(−1/2) θ per k_y) of the objective;
- Strehl ratio: the same spectrum with the analytic phase of the plate,
  exp(i(−k_z,air·x + k_z,IMC·(d − 1 µm))), and t_s of the front surface.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("fieldtracing")

from bwpm.pilot_imc import (  # noqa: E402
    IMC_THICKNESS,
    N_CLAD,
    N_CORE,
    N_IMC,
    NA_OBJECTIVE,
    make_grid,
    run_configuration,
    strehl_from_spectrum,
)

GRID = make_grid()
K = GRID.k0
AIR_FOCUS_DEPTH = 32.42e-3  # puts the best focus behind the cap onto the facet (measured)


def _aplanatic() -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    ky = GRID.ky[np.abs(GRID.ky) <= K * NA_OBJECTIVE]
    kz = np.sqrt(K * K - ky * ky)
    return ky, kz, (kz / K) ** -0.5


def _fresnel_power(n: float) -> float:
    _, kz, a = _aplanatic()
    kz2 = np.sqrt((n * K) ** 2 - (K * K - kz * kz))
    r = (kz - kz2) / (kz + kz2)
    return float(np.sum(a**2 * kz * np.abs(r) ** 2) / np.sum(a**2 * kz))


@pytest.fixture(scope="module")
def bare() -> object:
    return run_configuration(False, 0.0, GRID)


@pytest.fixture(scope="module")
def capped() -> object:
    return run_configuration(True, AIR_FOCUS_DEPTH, GRID)


def test_bare_fibre(bare: object) -> None:
    # reflection from the facet: between Fresnel for cladding and core; refocuses on the facet
    assert _fresnel_power(N_CLAD) < bare.reflected_power < _fresnel_power(N_CORE)  # type: ignore[attr-defined]
    assert abs(bare.reflected_focus) < 0.25e-3  # type: ignore[attr-defined]


def test_cap_moves_the_reflection_out_of_focus(capped: object) -> None:
    r = capped
    front = -IMC_THICKNESS
    # reflection at the IMC front: Fresnel air -> IMC
    assert abs(r.reflected_power - _fresnel_power(N_IMC)) < 2e-4  # type: ignore[attr-defined]
    # its focus is the mirror image of the air focus: 82.4 µm in front of the facet
    assert abs(r.reflected_focus - (front - AIR_FOCUS_DEPTH)) < 0.25e-3  # type: ignore[attr-defined]
    # backward power between cap and facet: reflection IMC -> fibre only (core/cladding bounds)
    zp, prof = r.z_planes, r.backward_profile  # type: ignore[attr-defined]
    inside = prof[np.argmin(np.abs(zp - (-20e-3)))]
    core = ((N_IMC - N_CORE) / (N_IMC + N_CORE)) ** 2
    clad = ((N_IMC - N_CLAD) / (N_IMC + N_CLAD)) ** 2
    assert 0.5 * core < inside < clad


def test_focus_quality_behind_the_cap(capped: object) -> None:
    ky, spectrum = capped.forward_spectrum_in_imc  # type: ignore[attr-defined]
    # k_y is conserved at planar interfaces: the physical spectrum is the objective's NA; the
    # 0.1 % of power outside is leakage from the absorbing mask and is left out of the Strehl ratio
    support = np.abs(ky) <= K * NA_OBJECTIVE
    offsets = np.arange(-5e-3, 5e-3, 0.02e-3)
    strehl, best = strehl_from_spectrum(ky[support], spectrum[support], K * N_IMC, offsets)
    # analytic reference: aplanatic spectrum, t_s of the front, plate phase
    kyr, kz, a = _aplanatic()
    kz_g = np.sqrt((N_IMC * K) ** 2 - kyr * kyr)
    t_s = 2 * kz / (kz + kz_g)
    reference = a * t_s * np.exp(1j * (-kz * AIR_FOCUS_DEPTH + kz_g * (IMC_THICKNESS - 1e-3)))
    strehl_ref, best_ref = strehl_from_spectrum(kyr, reference, K * N_IMC, offsets)
    assert abs(strehl - strehl_ref) < 1e-3
    assert abs(best - best_ref) < 0.05e-3
    assert abs(best - 1e-3) < 0.1e-3  # best focus on the facet (1 µm behind the sampling plane)
