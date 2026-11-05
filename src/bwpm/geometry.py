"""Material labels on the pixel grid from a geometry (decision E4, docs/spezifikation.md §7).

Pixel (j, n) covers z ∈ [z₀ + j·dz, z₀ + (j + 1)·dz] and y ∈ [y_n − dy/2, y_n + dy/2] (mm).

- `labels_at_centres`: standard rule, the material at the pixel centre (sharp staircase).
- `labels_area_fraction`: option for dielectric interfaces, the permittivity averaged over the
  pixel area (ε_eff = f·n₂² + (1 − f)·n₁², sub-sampled) and rounded to a few levels, so that the
  number of distinct indices M stays small.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from bwpm._types import FloatArray, IntArray
from bwpm.fpwpm import Grid2D

Inside = Callable[[FloatArray, FloatArray], np.ndarray[Any, Any]]


def labels_at_centres(grid: Grid2D, nz: int, z0: float, inside: Inside) -> IntArray:
    """Label 1 where the pixel centre lies inside the body, 0 elsewhere, shape (nz, ny).

    Args:
        grid: transverse grid and dz.
        nz: number of slabs.
        z0: position of the input plane (mm).
        inside: boolean function of (y, z) arrays in mm.
    """
    z = z0 + (np.arange(nz) + 0.5) * grid.dz
    zz, yy = np.meshgrid(z, grid.y, indexing="ij")
    return np.asarray(inside(yy, zz), dtype=np.int64)


def labels_area_fraction(
    grid: Grid2D,
    nz: int,
    z0: float,
    inside: Inside,
    index_out: complex,
    index_in: complex,
    levels: int = 8,
    subsamples: int = 8,
) -> tuple[list[complex], IntArray]:
    """Mixed indices from the area fraction f of each pixel inside the body.

    ε_eff = f·ñ_in² + (1 − f)·ñ_out², f rounded to `levels` equidistant values in (0, 1) plus the
    pure materials. Only meaningful for dielectric interfaces (decision E4).

    Returns:
        (indices, labels): index list (label 0 = outside, last = inside) and labels (nz, ny).
    """
    if levels < 1 or subsamples < 1:
        raise ValueError("levels and subsamples must be >= 1")
    offsets = (np.arange(subsamples) + 0.5) / subsamples - 0.5
    z = z0 + (np.arange(nz) + 0.5) * grid.dz
    fraction = np.zeros((nz, grid.ny))
    for oz in offsets:
        for oy in offsets:
            zz, yy = np.meshgrid(z + oz * grid.dz, grid.y + oy * grid.dy, indexing="ij")
            fraction += np.asarray(inside(yy, zz), dtype=np.float64)
    fraction /= subsamples * subsamples
    steps = np.rint(fraction * (levels + 1)).astype(np.int64)  # 0 … levels + 1
    f_values = np.arange(levels + 2) / (levels + 1)
    eps = f_values * complex(index_in) ** 2 + (1.0 - f_values) * complex(index_out) ** 2
    indices = [complex(np.sqrt(e)) for e in eps]
    return indices, np.asarray(steps)
