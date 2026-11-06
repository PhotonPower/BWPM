"""Fusion adapter of the BWPM (stage B4e).

Interface after the Fusion plan (Fieldtracing docs/rtt-fusion/plan.md, section 6.4) and
docs/b4-entwurf.md, section 4:

- `describe()`: method, conventions, ports, properties and validity limits as data.
- `prepare(structure, wavelengths)`: checks the structure, computes a content hash, returns a
  handle (no computation yet).
- `apply(handle, port_field, ...)`: field incident on the front port (+z) -> reflected field at
  the front port (−z) and transmitted field at the back port (+z).
- `estimate(handle, iterations, backend)`: run time and memory from the measurements of B3/B4d.
- `port_power(port_field)`: Poynting flux through a port (3D), the basis of the power
  normalisation against Raytatouille (Fieldtracing ADR 0004).

Ports (pilot, docs/pilot-imc-2d.md): plane z = const, direction +1/−1, medium at the port,
vacuum wavelength, transverse field (E_x, E_y) on the centred grid
x_m = (m − nx//2)·dx, y_n = (n − ny//2)·dy; spectra are referred to (0, 0).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from bwpm import __version__, fast
from bwpm._types import ComplexArray, IntArray
from bwpm.interface import longitudinal_wavenumber
from bwpm.volume3d import Bidirectional3D, Grid3D, absorbing_mask_2d, memory_estimate

# seconds per grid point and plane pass (one direction, one iteration has two passes), measured
# on the maintainer's machine (docs/testfaelle.md: B3 NumPy one core, B4d CPU core 32 cores, CUDA)
SECONDS_PER_POINT_PASS = {"python": 4.6e-6, "cpu": 1.5e-7, "cuda": 4.6e-8}


@dataclass(frozen=True)
class Port:
    """A port plane.

    Attributes:
        z: position of the plane (mm).
        direction: +1 for fields travelling in +z, −1 for −z.
        index: complex refractive index of the medium at the port.
        wavelength_um: vacuum wavelength (µm).
        nx, ny, dx, dy: centred transverse grid (mm).
    """

    z: float
    direction: int
    index: complex
    wavelength_um: float
    nx: int
    ny: int
    dx: float
    dy: float

    def __post_init__(self) -> None:
        if self.direction not in (1, -1):
            raise ValueError("direction must be +1 or -1")


@dataclass(frozen=True)
class PortField:
    """Transverse field (E_x, E_y), shape (ny, nx, 2), at a port."""

    port: Port
    field: ComplexArray

    def __post_init__(self) -> None:
        if np.asarray(self.field).shape != (self.port.ny, self.port.nx, 2):
            raise ValueError("field must have shape (ny, nx, 2) of the port grid")


@dataclass(frozen=True)
class Structure:
    """Index volume between the front plane z0 and the back plane z0 + nz·dz.

    Attributes:
        nx, ny, dx, dy, dz: grid (mm); slab j covers [z0 + j·dz, z0 + (j + 1)·dz].
        indices: complex index per material label.
        labels: label per pixel, shape (nz, ny, nx).
        index_before: medium in front of the volume (front port).
        z0: position of the front plane (mm).
        mask_width: width of the absorbing boundary zone (mm), 0 for none (rule K7: ≥ 20 λ).
    """

    nx: int
    ny: int
    dx: float
    dy: float
    dz: float
    indices: tuple[complex, ...]
    labels: IntArray
    index_before: complex
    z0: float = 0.0
    mask_width: float = 0.0


@dataclass
class Handle:
    """Prepared structure; `content_hash` identifies structure, wavelengths and BWPM version."""

    structure: Structure
    wavelengths_um: tuple[float, ...]
    content_hash: str
    _solvers: dict[float, Bidirectional3D] = field(default_factory=dict, repr=False)


def describe() -> dict[str, Any]:
    """Self-description for the Fusion planner (plan 6.4)."""
    return {
        "name": "BWPM",
        "version": __version__,
        "method": "fast bidirectional vector wave propagation (wende2023), FPWPM steps with "
        "Fresnel interface matrices, forward-backward iterations",
        "conventions": {
            "time": "exp(-i omega t), plane waves exp(i(k.r - omega t))",
            "fresnel": "Convention A (Byrnes arXiv:1603.02720v5 eq. 6)",
            "polarisation": "s = TE, p = TM; port fields are transverse (E_x, E_y)",
            "units": {"length": "mm", "wavelength": "um", "k": "1/mm"},
            "grid": "centred, x_m = (m - nx//2) dx; spectra referred to (0, 0)",
            "kz": "Im k_z >= 0; grazing rule |n^2 k0^2 - kperp^2| <= 1e-12 |n^2 k0^2|",
        },
        "ports": [
            {"name": "front", "plane": "z0", "in": "+z", "out": "-z (reflected)"},
            {"name": "back", "plane": "z0 + nz dz", "out": "+z (transmitted)"},
        ],
        "representation": "gridded transverse field (Fusion plan: G)",
        "linear": True,
        "reciprocal": True,
        "operator": False,
        "deterministic": {"python": "bitwise", "cpu": "bitwise", "cuda": "rel 1e-10"},
        "backends": _backends(),
        "validity": [
            "reflections quantitative only at z-normal interfaces; tilted surfaces: direction "
            "and phase only (E5, K8, K9)",
            "transmitted power at tilted surfaces only to a few percent (K8)",
            "waves trapped by total internal reflection do not converge; "
            "limit iterations (K4, K11)",
            "absorbing boundary zone >= 20 wavelengths, order 8 (K7)",
        ],
    }


def _backends() -> list[str]:
    out = ["python"]
    if fast.CORE_AVAILABLE:
        out.append("cpu")
    if fast.cuda_available():
        out.append("cuda")
    return out


def prepare(structure: Structure, wavelengths_um: list[float] | tuple[float, ...]) -> Handle:
    """Check the structure and return a handle with a content hash (nothing is computed)."""
    labels = np.asarray(structure.labels)
    if labels.ndim != 3 or labels.shape[1:] != (structure.ny, structure.nx):
        raise ValueError("labels must have shape (nz, ny, nx)")
    if labels.min() < 0 or labels.max() >= len(structure.indices):
        raise ValueError("labels must index into `indices`")
    if not wavelengths_um:
        raise ValueError("at least one wavelength is required")
    h = hashlib.sha256()
    h.update(f"bwpm {__version__}".encode())
    for value in (
        structure.nx,
        structure.ny,
        structure.dx,
        structure.dy,
        structure.dz,
        structure.index_before,
        structure.z0,
        structure.mask_width,
    ):
        h.update(repr(value).encode())
    h.update(np.asarray(structure.indices, dtype=np.complex128).tobytes())
    h.update(np.ascontiguousarray(labels, dtype=np.int32).tobytes())
    h.update(np.asarray(wavelengths_um, dtype=np.float64).tobytes())
    return Handle(structure, tuple(float(w) for w in wavelengths_um), h.hexdigest())


def _solver(handle: Handle, wavelength_um: float) -> Bidirectional3D:
    if wavelength_um not in handle.wavelengths_um:
        raise ValueError("wavelength was not prepared")
    if wavelength_um not in handle._solvers:
        s = handle.structure
        grid = Grid3D(s.nx, s.ny, s.dx, s.dy, s.dz, wavelength_um)
        mask = absorbing_mask_2d(grid, s.mask_width) if s.mask_width > 0 else None
        handle._solvers[wavelength_um] = Bidirectional3D(
            grid, list(s.indices), np.asarray(s.labels), s.index_before, mask
        )
    return handle._solvers[wavelength_um]


def front_port(handle: Handle, wavelength_um: float, direction: int = 1) -> Port:
    s = handle.structure
    return Port(s.z0, direction, s.index_before, wavelength_um, s.nx, s.ny, s.dx, s.dy)


def back_port(handle: Handle, wavelength_um: float) -> Port:
    s = handle.structure
    labels = np.asarray(s.labels)
    last = s.indices[int(labels[-1, s.ny // 2, s.nx // 2])]
    nz = labels.shape[0]
    return Port(s.z0 + nz * s.dz, 1, last, wavelength_um, s.nx, s.ny, s.dx, s.dy)


def apply(
    handle: Handle,
    incident: PortField,
    iterations: int = 3,
    backend: str = "auto",
) -> dict[str, PortField]:
    """Propagate a field incident on the front port.

    Args:
        handle: from `prepare`.
        incident: field at the front port with direction +1 and the prepared wavelength.
        iterations: forward-backward iterations (the source uses 1–10; 0 = forward only).
        backend: "python" (reference), "cpu", "cuda" or "auto" (fastest available).

    Returns:
        {"front": reflected field (direction −1), "back": transmitted field (direction +1)}.
    """
    p = incident.port
    s = handle.structure
    if p.direction != 1 or abs(p.z - s.z0) > 1e-15 or complex(p.index) != complex(s.index_before):
        raise ValueError("incident field must be on the front port, direction +1")
    if (p.nx, p.ny, p.dx, p.dy) != (s.nx, s.ny, s.dx, s.dy):
        raise ValueError("port grid differs from the structure grid")
    solver = _solver(handle, p.wavelength_um)
    if backend == "python" or (backend == "auto" and not fast.CORE_AVAILABLE):
        result = solver.run(np.asarray(incident.field), iterations)
    else:
        result = fast.run_3d(solver, np.asarray(incident.field), iterations, backend=backend)
    return {
        "front": PortField(front_port(handle, p.wavelength_um, -1), result.reflected),
        "back": PortField(back_port(handle, p.wavelength_um), result.transmitted),
    }


def estimate(handle: Handle, iterations: int = 3, backend: str = "auto") -> dict[str, float]:
    """Run time (s) and peak memory (bytes) from the measured models (B3, B4d)."""
    s = handle.structure
    if backend == "auto":
        backend = _backends()[-1]
    labels = np.asarray(s.labels)
    nz = labels.shape[0]
    passes = 2 * max(iterations, 1) if iterations > 0 else 1
    points = s.nx * s.ny * nz
    pairs = (
        len(
            {(int(a), int(b)) for a, b in zip(labels[:-1].ravel(), labels[1:].ravel(), strict=True)}
        )
        if nz > 1
        else 1
    )
    interfaces = int(np.sum(labels[:-1] != labels[1:])) if nz > 1 else 0
    return {
        "seconds": SECONDS_PER_POINT_PASS[backend] * points * passes,
        "memory_bytes": float(memory_estimate(s.nx, s.ny, interfaces, pairs, nz=nz, label_bytes=4)),
    }


def port_power(port_field: PortField) -> float:
    """Poynting flux through the port, per unit area normalisation of the grid (arbitrary units).

    Σ over plane waves (Parseval) of Re(k_z)·|E|², E_z = −direction·k⊥·E⊥/k_z, times dx·dy/N;
    evanescent waves carry no power. Valid for a homogeneous medium at the port.
    """
    p = port_field.port
    grid = Grid3D(p.nx, p.ny, p.dx, p.dy, 1.0, p.wavelength_um)
    kx, ky = grid.k_grid()
    spectrum = np.fft.fft2(np.asarray(port_field.field), axes=(0, 1))
    kz = longitudinal_wavenumber(p.index, grid.k0, kx * kx + ky * ky)
    with np.errstate(divide="ignore", invalid="ignore"):
        ez = np.where(
            kz != 0, -p.direction * (kx * spectrum[..., 0] + ky * spectrum[..., 1]) / kz, 0.0
        )
    magnitude = np.sum(np.abs(spectrum) ** 2, axis=-1) + np.abs(ez) ** 2
    n = p.nx * p.ny
    return float(np.sum(np.real(kz) * magnitude) * p.dx * p.dy / n)


def plane_wave(port: Port, mx: int, my: int, e_perp: tuple[complex, complex]) -> PortField:
    """Plane wave on the port grid with k⊥ = 2π(mx/(nx dx), my/(ny dy)), phase 1 at (0, 0)."""
    x = (np.arange(port.nx) - port.nx // 2) * port.dx
    y = (np.arange(port.ny) - port.ny // 2) * port.dy
    kx = 2 * np.pi * mx / (port.nx * port.dx)
    ky = 2 * np.pi * my / (port.ny * port.dy)
    phase = np.exp(1j * (kx * x[None, :] + ky * y[:, None]))
    return PortField(port, np.stack([e_perp[0] * phase, e_perp[1] * phase], axis=-1))


__all__ = [
    "Handle",
    "Port",
    "PortField",
    "Structure",
    "apply",
    "back_port",
    "describe",
    "estimate",
    "front_port",
    "plane_wave",
    "port_power",
    "prepare",
]
