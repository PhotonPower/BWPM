"""Reference cases for the C++/CUDA core (stage B4), computed with the Python reference.

Every case stores inputs, the values of this implementation (`expected`, for the port: same
algorithm, tight tolerance) and, where available, an independent value (`reference`, exact,
with the physical tolerance). Small grids, so the files stay small. Written by
tools/export_reference.py to reference/<case>.result.json; checked by tests/test_reference_files.py.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D, absorbing_mask
from bwpm.interface import longitudinal_wavenumber, reflection_matrix, transmission_matrix
from bwpm.results import field_reference
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d

SILVER = 0.213 + 6.697j
TIGHT = {"rel": 1e-12, "abs": 1e-14}


def _tight(values: np.ndarray) -> dict[str, float]:  # type: ignore[type-arg]
    return {"rel": 1e-12, "abs": 1e-13 * float(np.max(np.abs(values)) or 1.0)}


def interface_matrices() -> dict[str, Any]:
    k0 = 2 * np.pi / 0.5e-3
    rng = np.random.default_rng(5)
    theta = rng.uniform(0, np.deg2rad(85), 24)
    phi = rng.uniform(0, 2 * np.pi, 24)
    pairs = np.array([[1.0, 1.5], [1.5, 1.0], [1.0, SILVER], [1.505, SILVER], [1.5, 1.511]])
    t_all, r_all = [], []
    kx = k0 * np.sin(theta) * np.cos(phi)
    ky = k0 * np.sin(theta) * np.sin(phi)
    kx = np.concatenate([kx, [0.0]])  # normal incidence (limit)
    ky = np.concatenate([ky, [0.0]])
    for n1, n2 in pairs:
        t_all.append(transmission_matrix(kx * n1.real, ky * n1.real, n1, n2, k0))
        r_all.append(reflection_matrix(kx * n1.real, ky * n1.real, n1, n2, k0))
    t, r = np.asarray(t_all), np.asarray(r_all)
    return field_reference(
        "interface_matrices",
        "Transmission and reflection matrices of a plane interface for (E_x, E_y), Convention A.",
        "tests/test_interface.py (K1); wende2023 eqs. (2), (3), (8), (9)",
        {"k0_per_mm": k0, "note": "k⊥ of pair j is (k_x, k_y)·Re(n_in); last sample k⊥ = 0"},
        {"k_x_per_mm": kx, "k_y_per_mm": ky, "index_pairs": pairs},
        {"T": t, "R": r},
        _tight(np.concatenate([t.ravel(), r.ravel()])),
    )


def _grid2d() -> Grid2D:
    return Grid2D(ny=64, dy=0.5e-3 / 16, dz=0.5e-3 / 16, wavelength_um=0.5)


def fpwpm2d_interface() -> dict[str, Any]:
    g = _grid2d()
    nz = 24
    labels = np.zeros((nz, g.ny), dtype=np.int64)
    labels[8:] = 1
    field = np.zeros((g.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((g.y / (6 * g.dy)) ** 2))
    field[:, 1] = 0.4j * np.exp(
        -((g.y / (6 * g.dy)) ** 2) + 1j * 3 * 2 * np.pi * g.y / (g.ny * g.dy)
    )
    out = Propagator2D(g, [1.0 + 0j, 1.5 + 0.01j], labels, 1.0 + 0j).run(field, keep=True)
    return field_reference(
        "fpwpm2d_interface",
        "Forward FPWPM in 2D through a plane interface air -> weakly absorbing glass.",
        "tests/test_fpwpm.py (K3); wende2023 eq. (1), table 1 i = 1",
        {"ny": g.ny, "dy_mm": g.dy, "dz_mm": g.dz, "wavelength_um": g.wavelength_um},
        {
            "indices": np.array([1.0, 1.5 + 0.01j]),
            "labels": labels,
            "index_before": 1.0 + 0j,
            "field": field,
        },
        {"planes": out},
        _tight(out),
    )


def bwpm2d_tilted_silver() -> dict[str, Any]:
    g = _grid2d()
    nz = 32
    z = (np.arange(nz) + 0.5) * g.dz
    zz, yy = np.meshgrid(z, g.y, indexing="ij")
    labels = np.where(zz > 12 * g.dz + 0.3 * yy, 1, 0).astype(np.int64)
    labels[zz > 24 * g.dz] = 2
    field = np.zeros((g.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(-((g.y / (8 * g.dy)) ** 2))
    field[:, 1] = 0.3 * np.exp(-((g.y / (8 * g.dy)) ** 2))
    mask = absorbing_mask(g, 8 * g.dy, 8)
    prop = Propagator2D(g, [1.0 + 0j, 1.505 + 0j, SILVER], labels, 1.0 + 0j, mask)
    res = Bidirectional2D(prop).run(field, iterations=2, keep=True)
    assert res.forward is not None and res.backward is not None
    return field_reference(
        "bwpm2d_tilted_silver",
        "BWPM in 2D, 2 iterations: tilted glass interface and silver, with absorbing mask.",
        "tests/test_bidirectional.py, tests/test_mirror.py; wende2023 eqs. (6), (7), (10), (11)",
        {
            "ny": g.ny,
            "dy_mm": g.dy,
            "dz_mm": g.dz,
            "wavelength_um": g.wavelength_um,
            "iterations": 2,
            "mask_width_mm": 8 * g.dy,
            "mask_order": 8,
        },
        {
            "indices": np.array([1.0, 1.505, SILVER]),
            "labels": labels,
            "index_before": 1.0 + 0j,
            "field": field,
            "mask": mask,
        },
        {
            "forward": res.forward,
            "backward": res.backward,
            "energy_residual": np.asarray(res.energy_residual),
        },
        _tight(np.concatenate([res.forward.ravel(), res.backward.ravel()])),
    )


def bwpm2d_slab_exact() -> dict[str, Any]:
    g = Grid2D(ny=64, dy=1e-3 / 20, dz=1e-3 / 20, wavelength_um=1.0)
    media = [1.0 + 0j] * 5 + [1.5 + 0j] * 50 + [1.0 + 0j] * 5
    labels = np.array([[0 if n == 1.0 else 1] * g.ny for n in media], dtype=np.int64)
    ky = 2 * np.pi * 2 / (g.ny * g.dy)
    field = np.zeros((g.ny, 2), dtype=np.complex128)
    field[:, 0] = np.exp(1j * ky * g.y)
    field[:, 1] = 0.5 * np.exp(1j * ky * g.y)
    res = Bidirectional2D(Propagator2D(g, [1.0 + 0j, 1.5 + 0j], labels, 1.0 + 0j)).run(field, 40)
    # independent: Airy for TE (component x); TM via tangential transfer matrix is in test K4
    k0 = g.k0
    k1 = longitudinal_wavenumber(1.0, k0, np.array(ky * ky))
    k2 = longitudinal_wavenumber(1.5, k0, np.array(ky * ky))
    r12 = (k1 - k2) / (k1 + k2)
    phase = np.exp(2j * k2 * 50 * g.dz)
    r_airy = (r12 - r12 * phase) / (1 - r12 * r12 * phase) * np.exp(2j * k1 * 5 * g.dz)
    return field_reference(
        "bwpm2d_slab_exact",
        "BWPM in 2D, 40 iterations: slab n = 1.5 (2.5 µm) at 39°, converged to the exact solution.",
        "tests/test_bidirectional.py (K4); wende2023 §2.1.3",
        {"ny": g.ny, "dy_mm": g.dy, "dz_mm": g.dz, "wavelength_um": 1.0, "iterations": 40},
        {
            "indices": np.array([1.0, 1.5]),
            "labels": labels,
            "index_before": 1.0 + 0j,
            "field": field,
        },
        {"transmitted": res.transmitted, "reflected": res.reflected},
        _tight(np.concatenate([res.transmitted.ravel(), res.reflected.ravel()])),
        {"r_te_airy": complex(np.asarray(r_airy).item())},
        {"abs": 1e-12},
    )


def bwpm3d_tilted() -> dict[str, Any]:
    g = Grid3D(nx=16, ny=16, dx=0.5e-3 / 8, dy=0.5e-3 / 8, dz=0.5e-3 / 8, wavelength_um=0.5)
    nz = 20
    z = (np.arange(nz) + 0.5) * g.dz
    zz, yy, xx = np.meshgrid(z, g.y, g.x, indexing="ij")
    labels = np.where(zz > 8 * g.dz + 0.2 * yy - 0.1 * xx, 1, 0).astype(np.int64)
    labels[zz > 15 * g.dz] = 2
    r2 = g.x[None, :] ** 2 + g.y[:, None] ** 2
    field = np.stack([np.exp(-r2 / (4 * g.dx) ** 2), 0.2j * np.exp(-r2 / (4 * g.dx) ** 2)], -1)
    mask = absorbing_mask_2d(g, 3 * g.dx)
    res = Bidirectional3D(g, [1.0 + 0j, 1.5 + 0j, SILVER], labels, 1.0 + 0j, mask).run(
        field, iterations=2, keep=(0, 10, nz)
    )
    fwd = np.stack([res.planes[p][0] for p in (0, 10, nz)])
    bwd = np.stack([res.planes[p][1] for p in (0, 10, nz)])
    return field_reference(
        "bwpm3d_tilted",
        "BWPM in 3D, 2 iterations: tilted glass and silver, absorbing mask, planes 0/10/nz.",
        "tests/test_volume3d.py (K10-K12)",
        {
            "nx": g.nx,
            "ny": g.ny,
            "dx_mm": g.dx,
            "dy_mm": g.dy,
            "dz_mm": g.dz,
            "wavelength_um": 0.5,
            "iterations": 2,
            "mask_width_mm": 3 * g.dx,
            "mask_order": 8,
        },
        {
            "indices": np.array([1.0, 1.5, SILVER]),
            "labels": labels,
            "index_before": 1.0 + 0j,
            "field": field,
            "mask": mask,
        },
        {
            "transmitted": res.transmitted,
            "reflected": res.reflected,
            "forward_planes": fwd,
            "backward_planes": bwd,
        },
        _tight(np.concatenate([res.transmitted.ravel(), res.reflected.ravel()])),
    )


CASES: dict[str, Callable[[], dict[str, Any]]] = {
    f.__name__: f
    for f in (
        interface_matrices,
        fpwpm2d_interface,
        bwpm2d_tilted_silver,
        bwpm2d_slab_exact,
        bwpm3d_tilted,
    )
}
