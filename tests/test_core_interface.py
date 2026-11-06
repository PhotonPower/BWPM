"""B4b: C++ interface matrices (bwpm._core) against the reference file and the Python reference.

Skipped when the core is not built, except when BWPM_REQUIRE_CORE=1 (CI), where a missing core
fails.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from bwpm import fast
from bwpm.interface import reflection_matrix, transmission_matrix
from bwpm.results import load_json

ROOT = Path(__file__).resolve().parents[1]

if not fast.CORE_AVAILABLE:
    if os.environ.get("BWPM_REQUIRE_CORE") == "1":
        raise RuntimeError("BWPM_REQUIRE_CORE=1 but bwpm._core is not built")
    pytest.skip("bwpm._core not built", allow_module_level=True)

from bwpm import _core  # type: ignore[attr-defined]  # noqa: E402


def test_reference_file() -> None:
    _, data = load_json((ROOT / "reference" / "interface_matrices.result.json").read_text("utf-8"))
    kx, ky = data["inputs"]["k_x_per_mm"], data["inputs"]["k_y_per_mm"]
    k0 = data["parameters"]["k0_per_mm"]
    tol = data["tolerance"]
    for j, (n1, n2) in enumerate(data["inputs"]["index_pairs"]):
        scale = n1.real
        t = _core.transmission_matrix(kx * scale, ky * scale, n1, n2, k0)
        r = _core.reflection_matrix(kx * scale, ky * scale, n1, n2, k0)
        assert np.allclose(t, data["expected"]["T"][j], rtol=tol["rel"], atol=tol["abs"])
        assert np.allclose(r, data["expected"]["R"][j], rtol=tol["rel"], atol=tol["abs"])


@pytest.mark.parametrize(("n1", "n2"), [(1.0, 1.5), (1.5, 1.0), (1.0, 0.213 + 6.697j), (1.3, 1.3)])
def test_equals_python_reference_on_a_grid(n1: complex, n2: complex) -> None:
    # FFT-like grid with k⊥ = 0, propagating, grazing (|k⊥| = n k0 exactly) and evanescent points
    k0 = 2 * np.pi / 0.5e-3
    k = np.arange(-8, 8) * k0 / 4
    kx, ky = np.meshgrid(k, k)
    for core_fn, ref_fn in (
        (_core.transmission_matrix, transmission_matrix),
        (_core.reflection_matrix, reflection_matrix),
    ):
        a = core_fn(kx, ky, n1, n2, k0)
        b = ref_fn(kx, ky, n1, n2, k0)
        assert a.shape == b.shape == (*kx.shape, 2, 2)
        assert np.all(np.isfinite(a))
        assert np.allclose(a, b, rtol=1e-13, atol=1e-13)


def test_facade_uses_core() -> None:
    assert fast.transmission_matrix is _core.transmission_matrix


def test_shape_mismatch_is_an_error() -> None:
    with pytest.raises(ValueError):
        _core.transmission_matrix(np.zeros(3), np.zeros(4), 1.0, 1.5, 1.0)
