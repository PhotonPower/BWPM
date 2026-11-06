"""B4: reference files for the C++/CUDA core (reference/*.result.json).

Every case of tools/reference_cases.py is recomputed and compared with its file (stale files fail
within the stored tolerance), independent values are checked within the physical tolerance, and
the encoding round-trips.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from bwpm.results import FIELD_REFERENCE, decode, encode, load_json, to_json

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from reference_cases import CASES  # noqa: E402


def _close(a: Any, b: Any, rel: float, abs_: float, where: str) -> None:
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), where
        for key in a:
            _close(a[key], b[key], rel, abs_, f"{where}/{key}")
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), where
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _close(x, y, rel, abs_, f"{where}/{i}")
    elif isinstance(a, np.ndarray):
        assert isinstance(b, np.ndarray) and a.dtype == b.dtype and a.shape == b.shape, where
        if a.dtype.kind in "fc":
            assert np.allclose(a, b, rtol=max(rel, 1e-10), atol=abs_, equal_nan=True), where
        else:
            assert np.array_equal(a, b), where
    elif isinstance(a, float | complex):
        assert abs(a - b) <= max(abs_, max(rel, 1e-10) * abs(b)) or (
            isinstance(a, float) and math.isinf(a) and a == b
        ), where
    else:
        assert a == b, where


@pytest.mark.parametrize("name", list(CASES))
def test_file_is_up_to_date(name: str) -> None:
    path = ROOT / "reference" / f"{name}.result.json"
    type_name, stored = load_json(path.read_text(encoding="utf-8"))
    assert type_name == FIELD_REFERENCE
    fresh = decode(encode(CASES[name]()))
    tol = fresh["tolerance"]
    _close(fresh["expected"], stored["expected"], tol["rel"], tol["abs"], f"{name}/expected")
    for key in ("case", "description", "source", "parameters", "inputs", "tolerance"):
        _close(fresh[key], stored[key], 1e-12, 0.0, f"{name}/{key}")


def test_independent_values_hold() -> None:
    _, data = load_json((ROOT / "reference" / "bwpm2d_slab_exact.result.json").read_text("utf-8"))
    r = complex(np.mean(data["expected"]["reflected"][:, 0] / data["inputs"]["field"][:, 0]))
    assert abs(r - data["reference"]["r_te_airy"]) < data["reference_tolerance"]["abs"]


@pytest.mark.parametrize("name", list(CASES))
def test_round_trip(name: str) -> None:
    path = ROOT / "reference" / f"{name}.result.json"
    type_name, stored = load_json(path.read_text(encoding="utf-8"))
    again = load_json(to_json(type_name, stored))[1]
    _close(stored, again, 0.0, 0.0, name)
