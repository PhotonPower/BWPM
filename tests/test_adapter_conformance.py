"""B4e: Fusion adapter and conformance tests (Fusion plan 6.5, docs/b4-entwurf.md §5).

Independent references: tangential transfer matrix (s and p separately) and Airy formula for
planar stacks, energy conservation for lossless stacks, reciprocity of the transmission, and
determinism. Every test runs for every available backend (Python reference, CPU core, CUDA).
"""

from __future__ import annotations

import numpy as np
import pytest

from bwpm import adapter

BACKENDS = adapter.describe()["backends"]
LAMBDA_UM = 1.0
STEP = LAMBDA_UM * 1e-3 / 16
K0 = 2 * np.pi / (LAMBDA_UM * 1e-3)
SILVER = 0.213 + 6.697j


def _kz(n: complex, kp2: float) -> complex:
    q = complex(np.sqrt(complex(n * n * K0 * K0 - kp2)))
    return -q if (q.imag < 0 or (q.imag == 0 and q.real < 0)) else q


def _structure(media: list[complex], n: int = 32, before: complex = 1.0 + 0j) -> adapter.Structure:
    unique = sorted(set(media), key=lambda c: (c.real, c.imag))
    labels = np.array([np.full((n, n), unique.index(m)) for m in media], dtype=np.int64)
    return adapter.Structure(n, n, STEP, STEP, STEP, tuple(unique), labels, before)


def _transfer(
    media: list[complex], before: complex, kp2: float, pol: str
) -> tuple[complex, complex]:
    """(r, t) of the tangential field component, discretised stack (as tests/test_volume3d.py)."""
    layers = [before, *media]
    kz = [_kz(m, kp2) for m in layers]
    eta = [k if pol == "s" else m * m * K0 * K0 / k for m, k in zip(layers, kz, strict=True)]
    a_plus, a_minus = 1.0 + 0j, 0.0 + 0j
    for j in range(len(layers) - 1, 0, -1):
        e_t, h_t = a_plus + a_minus, eta[j] * (a_plus - a_minus)
        lp, lm = 0.5 * (e_t + h_t / eta[j - 1]), 0.5 * (e_t - h_t / eta[j - 1])
        if j - 1 == 0:
            a_plus, a_minus = lp, lm
        else:
            a_plus, a_minus = (
                lp * np.exp(-1j * kz[j - 1] * STEP),
                lm * np.exp(1j * kz[j - 1] * STEP),
            )
    return complex(a_minus / a_plus), complex(np.exp(1j * kz[-1] * STEP) / a_plus)


STACK = [1.0 + 0j] * 2 + [1.5 + 0j] * 9 + [2.1 + 0.02j] * 4 + [1.45 + 0j] * 6 + [1.0 + 0j] * 3


def test_describe_declares_conventions_and_limits() -> None:
    d = adapter.describe()
    assert d["conventions"]["time"].startswith("exp(-i omega t)")
    assert "Convention A" in d["conventions"]["fresnel"]
    assert d["conventions"]["polarisation"].startswith("s = TE, p = TM")
    assert d["linear"] and d["reciprocal"] and not d["operator"]
    assert len(d["validity"]) >= 4 and "python" in d["backends"]


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize(("mx", "my"), [(0, 0), (1, 0), (1, 1), (-1, 1)])
def test_planar_stack_against_transfer_matrix(backend: str, mx: int, my: int) -> None:
    handle = adapter.prepare(_structure(STACK), [LAMBDA_UM])
    e0 = (0.6 + 0.2j, -0.3 + 0.5j)
    incident = adapter.plane_wave(adapter.front_port(handle, LAMBDA_UM), mx, my, e0)
    out = adapter.apply(handle, incident, iterations=30, backend=backend)
    # 32 points at λ/16: |k⊥| = k0·sqrt(mx² + my²)/2 (30°, 45° at azimuth 45° and 135°)
    kx, ky = 2 * np.pi * mx / (32 * STEP), 2 * np.pi * my / (32 * STEP)
    kp2 = kx * kx + ky * ky
    s_dir = np.array([-ky, kx]) / np.sqrt(kp2) if kp2 else np.array([0.0, 1.0])
    p_dir = np.array([kx, ky]) / np.sqrt(kp2) if kp2 else np.array([1.0, 0.0])
    e = np.array(e0)
    r_s, t_s = _transfer(STACK, 1.0, kp2, "s")
    r_p, t_p = _transfer(STACK, 1.0, kp2, "p")
    expected_r = r_s * (e @ s_dir) * s_dir + r_p * (e @ p_dir) * p_dir
    expected_t = t_s * (e @ s_dir) * s_dir + t_p * (e @ p_dir) * p_dir
    centre = (16, 16)
    assert out["front"].port.direction == -1 and out["back"].port.direction == 1
    assert np.max(np.abs(out["front"].field[centre] - expected_r)) < 1e-12
    assert np.max(np.abs(out["back"].field[centre] - expected_t)) < 1e-12


@pytest.mark.parametrize("backend", BACKENDS)
def test_slab_against_airy(backend: str) -> None:
    # TE (E_x for k⊥ along y), slab n = 1.5 of 40 dz between 4 dz of vacuum on each side
    media = [1.0 + 0j] * 4 + [1.5 + 0j] * 40 + [1.0 + 0j] * 4
    handle = adapter.prepare(_structure(media), [LAMBDA_UM])
    incident = adapter.plane_wave(adapter.front_port(handle, LAMBDA_UM), 0, 1, (1.0, 0.0))
    out = adapter.apply(handle, incident, iterations=40, backend=backend)
    ky = 2 * np.pi / (32 * STEP)
    k1, k2 = _kz(1.0, ky * ky), _kz(1.5, ky * ky)
    r12 = (k1 - k2) / (k1 + k2)
    phase = np.exp(2j * k2 * 40 * STEP)
    r_airy = (r12 - r12 * phase) / (1 - r12 * r12 * phase) * np.exp(2j * k1 * 4 * STEP)
    assert abs(out["front"].field[16, 16, 0] - r_airy) < 1e-12


@pytest.mark.parametrize("backend", BACKENDS)
def test_energy_conservation_lossless(backend: str) -> None:
    media = [1.0 + 0j] * 2 + [1.5 + 0j] * 9 + [2.1 + 0j] * 4 + [1.45 + 0j] * 6 + [1.0 + 0j] * 3
    handle = adapter.prepare(_structure(media), [LAMBDA_UM])
    port = adapter.front_port(handle, LAMBDA_UM)
    # band-limited beam: propagating in vacuum (|k⊥| < 0.9 k0) so nothing is trapped in n = 2.1
    rng = np.random.default_rng(4)
    field = rng.normal(size=(32, 32, 2)) + 1j * rng.normal(size=(32, 32, 2))
    kx = 2 * np.pi * np.fft.fftfreq(32, STEP)
    kxx, kyy = np.meshgrid(kx, kx)
    spectrum = np.fft.fft2(field, axes=(0, 1))
    spectrum[np.sqrt(kxx**2 + kyy**2) >= 0.9 * K0] = 0.0
    incident = adapter.PortField(port, np.fft.ifft2(spectrum, axes=(0, 1)))
    out = adapter.apply(handle, incident, iterations=60, backend=backend)
    p_in = adapter.port_power(incident)
    balance = (adapter.port_power(out["front"]) + adapter.port_power(out["back"])) / p_in
    assert abs(balance - 1.0) < 1e-12


@pytest.mark.parametrize("backend", BACKENDS)
def test_reciprocity_of_transmission(backend: str) -> None:
    # asymmetric absorbing stack with vacuum on both sides: t is the same from both sides,
    # r is not (incidence from the back = the mirrored stack)
    media = [1.0 + 0j] * 2 + [1.5 + 0j] * 5 + [SILVER] * 1 + [2.1 + 0.05j] * 7 + [1.0 + 0j] * 3
    forward = adapter.prepare(_structure(media), [LAMBDA_UM])
    backward = adapter.prepare(_structure(media[::-1]), [LAMBDA_UM])
    results = []
    for handle in (forward, backward):
        incident = adapter.plane_wave(adapter.front_port(handle, LAMBDA_UM), 1, 1, (0.4, 0.7j))
        results.append(adapter.apply(handle, incident, iterations=40, backend=backend))
    t_fwd, t_bwd = results[0]["back"].field[16, 16], results[1]["back"].field[16, 16]
    r_fwd, r_bwd = results[0]["front"].field[16, 16], results[1]["front"].field[16, 16]
    assert np.max(np.abs(t_fwd - t_bwd)) < 1e-12
    assert np.max(np.abs(r_fwd - r_bwd)) > 1e-3


@pytest.mark.parametrize("backend", BACKENDS)
def test_determinism_and_hash(backend: str) -> None:
    structure = _structure(STACK)
    a, b = adapter.prepare(structure, [LAMBDA_UM]), adapter.prepare(structure, [LAMBDA_UM])
    assert a.content_hash == b.content_hash
    other = _structure(STACK[:-1] + [1.5 + 0j])
    assert adapter.prepare(other, [LAMBDA_UM]).content_hash != a.content_hash
    incident = adapter.plane_wave(adapter.front_port(a, LAMBDA_UM), 1, 0, (1.0, 0.5))
    one = adapter.apply(a, incident, 3, backend)
    two = adapter.apply(b, incident, 3, backend)
    tolerance = 0.0 if backend != "cuda" else 1e-10
    assert np.max(np.abs(one["back"].field - two["back"].field)) <= tolerance


@pytest.mark.parametrize("backend", BACKENDS)
def test_port_grid_is_centred(backend: str) -> None:
    # laterally homogeneous stack: shifting the input by whole pixels shifts the output (the
    # spectrum is referred to (0, 0) and nothing depends on the absolute index)
    handle = adapter.prepare(_structure(STACK), [LAMBDA_UM])
    port = adapter.front_port(handle, LAMBDA_UM)
    y = (np.arange(32) - 16) * STEP
    beam = np.exp(-((y[:, None] / (3 * STEP)) ** 2 + (y[None, :] / (3 * STEP)) ** 2))
    field = np.stack([beam, 0.3j * beam], axis=-1)
    a = adapter.apply(handle, adapter.PortField(port, field), 3, backend)
    b = adapter.apply(
        handle, adapter.PortField(port, np.roll(field, (2, -3), axis=(0, 1))), 3, backend
    )
    assert np.max(np.abs(np.roll(a["back"].field, (2, -3), axis=(0, 1)) - b["back"].field)) < 1e-13


def test_estimate_and_errors() -> None:
    handle = adapter.prepare(_structure(STACK), [LAMBDA_UM])
    est = adapter.estimate(handle, iterations=3, backend="cpu")
    assert est["seconds"] > 0 and est["memory_bytes"] > 0
    with pytest.raises(ValueError):
        adapter.prepare(_structure(STACK), [])
    wrong = adapter.Port(0.0, -1, 1.0, LAMBDA_UM, 32, 32, STEP, STEP)
    with pytest.raises(ValueError):
        adapter.apply(handle, adapter.PortField(wrong, np.zeros((32, 32, 2), complex)))
