"""Facade that uses the compiled core `bwpm._core` when it is built (decision E7, D1).

Same functions and conventions as the Python reference; without the core the reference is used.
`CORE_AVAILABLE` tells which one is active.
"""

from __future__ import annotations

from bwpm import interface as _reference

try:
    from bwpm import _core  # type: ignore[attr-defined]

    CORE_AVAILABLE = True
except ImportError:  # pure Python installation (BWPM_BUILD_CORE=OFF)
    _core = None
    CORE_AVAILABLE = False

if _core is not None:
    transmission_matrix = _core.transmission_matrix
    reflection_matrix = _core.reflection_matrix
else:
    transmission_matrix = _reference.transmission_matrix
    reflection_matrix = _reference.reflection_matrix

__all__ = ["CORE_AVAILABLE", "reflection_matrix", "transmission_matrix"]
