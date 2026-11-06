"""Facade that uses the compiled core `bwpm._core` when it is built (decision E7, D1).

Same functions, conventions and results as the Python reference (within rounding); without the
core the reference is used. `CORE_AVAILABLE` tells which one is active.

- `transmission_matrix`, `reflection_matrix`: interface matrices (B4b).
- `run_2d(propagator, field, iterations, keep)`: 2D FPWPM (iterations = 0) or BWPM, returns
  `bwpm.bidirectional.BidirectionalResult` like `Bidirectional2D.run` (B4c).
- `run_3d(solver, field, iterations, keep)`: like `Bidirectional3D.run`, returns `Result3D` (B4c).
"""

from __future__ import annotations

import os

import numpy as np

from bwpm import interface as _reference
from bwpm._types import ComplexArray
from bwpm.bidirectional import Bidirectional2D, BidirectionalResult
from bwpm.fpwpm import Propagator2D, power_flux
from bwpm.volume3d import Bidirectional3D, Result3D

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


def _threads(threads: int | None) -> int:
    return threads if threads and threads > 0 else (os.cpu_count() or 1)


def run_2d(
    propagator: Propagator2D,
    field: ComplexArray,
    iterations: int = 3,
    keep: bool = False,
    threads: int | None = None,
) -> BidirectionalResult:
    """2D run with the core (iterations = 0: forward FPWPM only), else with the reference."""
    g = propagator.grid
    if _core is None:
        if iterations == 0:
            planes = propagator.run(field, keep=keep)
            last = planes[-1] if keep else planes
            return BidirectionalResult(
                last, np.zeros_like(last), planes if keep else None, None, []
            )
        return Bidirectional2D(propagator).run(field, iterations, keep)
    labels = np.asarray(propagator.labels, dtype=np.int32)[:, :, None]
    nz = labels.shape[0]
    mask = None if propagator.mask is None else np.asarray(propagator.mask)[:, None]
    solver = _core.Solver(
        g.ny,
        1,
        g.dy,
        1.0,
        g.dz,
        g.k0,
        np.asarray(propagator.indices, dtype=np.complex128),
        labels,
        complex(propagator.index_before),
        mask,
        _threads(threads),
    )
    out = solver.run(np.asarray(field)[:, None, :], iterations, list(range(nz + 1)) if keep else [])
    forward = backward = None
    if keep:
        forward = np.stack([out["forward"][p][:, 0] for p in range(nz + 1)])
        if iterations > 0:
            backward = np.stack([out["backward"][p][:, 0] for p in range(nz + 1)])
    residual = []
    if iterations > 0:
        p_in = power_flux(np.asarray(field), g, propagator.index_before)
        last = complex(propagator.indices[int(np.asarray(propagator.labels)[-1][g.ny // 2])])
        for t, r in zip(out["transmitted_iter"], out["reflected_iter"], strict=True):
            p_r = power_flux(r[:, 0], g, propagator.index_before, direction=-1)
            p_t = power_flux(t[:, 0], g, last)
            residual.append(1.0 - (p_r + p_t) / p_in if p_in > 0 else float("nan"))
    return BidirectionalResult(
        transmitted=out["transmitted"][:, 0],
        reflected=out["reflected"][:, 0],
        forward=forward,
        backward=backward,
        energy_residual=residual,
    )


def run_3d(
    solver: Bidirectional3D,
    field: ComplexArray,
    iterations: int = 1,
    keep: tuple[int, ...] = (),
    threads: int | None = None,
) -> Result3D:
    """3D run with the core, else with the reference (`Bidirectional3D.run`)."""
    if _core is None:
        return solver.run(field, iterations, keep)
    g = solver.grid
    core = _core.Solver(
        g.ny,
        g.nx,
        g.dy,
        g.dx,
        g.dz,
        g.k0,
        np.asarray(solver.indices, dtype=np.complex128),
        np.asarray(solver.labels, dtype=np.int32),
        complex(solver.index_before),
        solver.mask,
        _threads(threads),
    )
    out = core.run(np.asarray(field), iterations, list(keep))
    zero = np.zeros_like(out["transmitted"])
    planes = {p: (out["forward"].get(p, zero), out["backward"].get(p, zero)) for p in keep}
    return Result3D(out["transmitted"], out["reflected"], planes, int(out["stored_values"]))


__all__ = ["CORE_AVAILABLE", "reflection_matrix", "run_2d", "run_3d", "transmission_matrix"]
