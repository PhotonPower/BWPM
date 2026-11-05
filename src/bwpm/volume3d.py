"""BWPM in 3D (x, y, z) with reflections stored only at interface pixels (stage B3).

Specification: docs/spezifikation.md, sections 4 and 5; source wende2023, eqs. (1), (6), (7),
(10), (11) and table 1. Media and positions as in `bwpm.bidirectional` (derived in its module
docstring and checked in test K4): planes p = 0 … nz at z₀ + p·dz, plane p belongs to slab p − 1,
interface p between slabs p − 1 and p.

Memory: the source of the paper stores the reflected field in the whole volume (N_x N_y N_z
vectors). Reflections only arise where the index changes between two slabs (I_r,s with r ≠ s), so
here the reflection sources are kept as (pixel indices, values) of those pixels only. Fields of
intermediate planes are kept only for requested planes.

Fields have shape (ny, nx, 2) with components (E_x, E_y); grids are centred,
x_m = (m − nx//2)·dx, y_n = (n − ny//2)·dy (mm).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from bwpm._types import ComplexArray, FloatArray, IntArray
from bwpm.interface import longitudinal_wavenumber, reflection_matrix, transmission_matrix


@dataclass(frozen=True)
class Grid3D:
    """Transverse grid (ny × nx) and step size; lengths in mm, wavelength in µm."""

    nx: int
    ny: int
    dx: float
    dy: float
    dz: float
    wavelength_um: float

    def __post_init__(self) -> None:
        if min(self.nx, self.ny) < 2 or min(self.dx, self.dy, self.dz, self.wavelength_um) <= 0:
            raise ValueError("nx, ny >= 2 and positive steps and wavelength are required")

    @property
    def k0(self) -> float:
        """Vacuum wavenumber in 1/mm."""
        return float(2.0 * np.pi / (self.wavelength_um * 1e-3))

    @property
    def x(self) -> FloatArray:
        return np.asarray((np.arange(self.nx) - self.nx // 2) * self.dx, dtype=np.float64)

    @property
    def y(self) -> FloatArray:
        return np.asarray((np.arange(self.ny) - self.ny // 2) * self.dy, dtype=np.float64)

    def k_grid(self) -> tuple[FloatArray, FloatArray]:
        """(k_x, k_y) of the 2D FFT, shape (ny, nx) each, 1/mm (numpy order)."""
        kx = 2.0 * np.pi * np.fft.fftfreq(self.nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.dy)
        kxx, kyy = np.meshgrid(kx, ky)
        return np.asarray(kxx), np.asarray(kyy)


def absorbing_mask_2d(grid: Grid3D, width: float, order: int = 8, edge: float = 1e-6) -> FloatArray:
    """Product of 1D super-Gaussian masks in x and y (decision E2, rule of K7), shape (ny, nx)."""
    if width < 0 or order < 1 or not 0 < edge < 1:
        raise ValueError("width >= 0, order >= 1 and 0 < edge < 1 are required")

    def one(coordinates: FloatArray) -> FloatArray:
        a = np.abs(coordinates)
        if width == 0:
            return np.ones_like(a)
        u = np.clip((a - (float(np.max(a)) - width)) / width, 0.0, 1.0)
        return np.asarray(edge ** (u**order))

    return np.asarray(one(grid.y)[:, None] * one(grid.x)[None, :])


@dataclass(frozen=True)
class Result3D:
    """Fields after the last iteration.

    Attributes:
        transmitted: forward field at plane nz, (ny, nx, 2).
        reflected: backward field at plane 0, (ny, nx, 2).
        planes: requested planes {p: (forward, backward)}.
        stored_values: number of stored reflection source vectors (memory diagnostic).
    """

    transmitted: ComplexArray
    reflected: ComplexArray
    planes: dict[int, tuple[ComplexArray, ComplexArray]]
    stored_values: int


@dataclass
class Bidirectional3D:
    """Forward FPWPM and BWPM in 3D.

    Attributes:
        grid: transverse grid and step.
        indices: complex index per material label.
        labels: label per pixel, shape (nz, ny, nx).
        index_before: medium of the input field.
        mask: absorbing mask (ny, nx) applied after every transmission step, or None.
        sparse: keep reflection sources only at interface pixels (default) or as full planes.
    """

    grid: Grid3D
    indices: list[complex]
    labels: IntArray
    index_before: complex
    mask: FloatArray | None = None
    sparse: bool = True
    _cache: dict[tuple[str, complex, complex], ComplexArray] = field(
        default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        lab = np.asarray(self.labels)
        g = self.grid
        if lab.ndim != 3 or lab.shape[1:] != (g.ny, g.nx):
            raise ValueError("labels must have shape (nz, ny, nx)")
        if lab.min() < 0 or lab.max() >= len(self.indices):
            raise ValueError("labels must index into `indices`")
        self._kx, self._ky = g.k_grid()

    def _index(self, label: int) -> complex:
        return self.index_before if label < 0 else self.indices[label]

    def _plane_labels(self, p: int) -> IntArray:
        if p == 0:
            return np.full((self.grid.ny, self.grid.nx), -1, dtype=np.int64)
        return np.asarray(self.labels[p - 1])

    def _matrix(self, kind: str, r: complex, s: complex) -> ComplexArray:
        key = (kind, complex(r), complex(s))
        if key not in self._cache:
            g, kx, ky = self.grid, self._kx, self._ky
            kp2 = kx * kx + ky * ky
            if kind == "fwd_t":  # table 1, i = 1
                phase = np.exp(1j * longitudinal_wavenumber(s, g.k0, kp2) * g.dz)
                m = transmission_matrix(kx, ky, r, s, g.k0) * phase[..., None, None]
            elif kind == "fwd_r":  # i = 2
                m = reflection_matrix(kx, ky, r, s, g.k0)
            elif kind == "bwd_t":  # i = 3
                phase = np.exp(1j * longitudinal_wavenumber(r, g.k0, kp2) * g.dz)
                m = transmission_matrix(kx, ky, r, s, g.k0) * phase[..., None, None]
            else:  # "bwd_r", i = 4
                phase = np.exp(2j * longitudinal_wavenumber(r, g.k0, kp2) * g.dz)
                m = reflection_matrix(kx, ky, r, s, g.k0) * phase[..., None, None]
            self._cache[key] = np.asarray(m)
        return self._cache[key]

    def _apply(
        self, kind: str, field_in: ComplexArray, labels_r: IntArray, labels_s: IntArray
    ) -> ComplexArray:
        spectrum = np.fft.fft2(field_in, axes=(0, 1))
        out = np.zeros_like(spectrum)
        pairs = np.unique(np.stack([labels_r.ravel(), labels_s.ravel()], axis=1), axis=0)
        for r, s in pairs:
            ir, is_ = self._index(int(r)), self._index(int(s))
            if kind.endswith("_r") and ir == is_:
                continue
            weight = (labels_r == r) & (labels_s == s)
            m = self._matrix(kind, ir, is_)
            values = np.fft.ifft2(np.einsum("...ij,...j->...i", m, spectrum), axes=(0, 1))
            out[weight] = values[weight]
        if self.mask is not None and kind.endswith("_t"):
            out = out * self.mask[..., None]
        return np.asarray(out)

    def _interfaces(self, p_left: int, p_right: int) -> np.ndarray:  # type: ignore[type-arg]
        """Boolean (ny, nx): pixels where the index changes between the two planes' media."""
        a, b = self._plane_labels(p_left), self._plane_labels(p_right)
        ia = np.vectorize(self._index)(a) if p_left == 0 else np.asarray(self.indices)[a]
        ib = np.asarray(self.indices)[b]
        return np.asarray(ia != ib)

    def _store(self, values: ComplexArray, where: np.ndarray) -> tuple[np.ndarray, ComplexArray]:  # type: ignore[type-arg]
        if self.sparse:
            flat = np.flatnonzero(where)
            return flat, np.asarray(values.reshape(-1, 2)[flat])
        return np.arange(where.size), np.asarray(values.reshape(-1, 2).copy())

    def _restore(self, stored: tuple[np.ndarray, ComplexArray] | None) -> ComplexArray:  # type: ignore[type-arg]
        g = self.grid
        out = np.zeros((g.ny * g.nx, 2), dtype=np.complex128)
        if stored is not None:
            out[stored[0]] = stored[1]
        return out.reshape(g.ny, g.nx, 2)

    def run(
        self, field_in: ComplexArray, iterations: int = 1, keep: tuple[int, ...] = ()
    ) -> Result3D:
        """Forward-backward iterations; `iterations = 0` runs the forward FPWPM only."""
        g = self.grid
        nz = np.asarray(self.labels).shape[0]
        e_in = np.asarray(field_in, dtype=np.complex128)
        if e_in.shape != (g.ny, g.nx, 2):
            raise ValueError("field must have shape (ny, nx, 2)")
        interface = {p: self._interfaces(p, p + 1) for p in range(nz)}
        src_fwd: dict[int, tuple[np.ndarray, ComplexArray]] = {}  # type: ignore[type-arg]
        planes: dict[int, tuple[ComplexArray, ComplexArray]] = {}
        forward_keep: dict[int, ComplexArray] = {}
        backward_keep: dict[int, ComplexArray] = {}
        transmitted: ComplexArray = e_in
        reflected: ComplexArray = np.zeros_like(e_in)
        stored = 0
        for it in range(max(iterations, 1)):
            src_bwd: dict[int, tuple[np.ndarray, ComplexArray]] = {}  # type: ignore[type-arg]
            e: ComplexArray = e_in + self._restore(src_fwd.get(0))
            if 0 in keep:
                forward_keep[0] = e.copy()
            for p in range(nz):
                lr, ls = self._plane_labels(p), self._plane_labels(p + 1)
                if iterations > 0 and interface[p].any():
                    src_bwd[p] = self._store(self._apply("fwd_r", e, lr, ls), interface[p])
                e = self._apply("fwd_t", e, lr, ls)
                if p + 1 in src_fwd:
                    e = e + self._restore(src_fwd[p + 1])
                if p + 1 in keep:
                    forward_keep[p + 1] = e.copy()
            transmitted = e
            if iterations == 0:
                break
            src_fwd = {}
            e = self._restore(src_bwd.get(nz))
            if nz in keep:
                backward_keep[nz] = e.copy()
            for p in range(nz, 0, -1):
                lr, ls = self._plane_labels(p), self._plane_labels(p - 1)
                if interface[p - 1].any():
                    src_fwd[p] = self._store(self._apply("bwd_r", e, lr, ls), interface[p - 1])
                e = self._apply("bwd_t", e, lr, ls)
                if p - 1 in src_bwd:
                    e = e + self._restore(src_bwd[p - 1])
                if p - 1 in keep:
                    backward_keep[p - 1] = e.copy()
            reflected = e
            stored = max(
                stored,
                sum(v[0].size for v in src_bwd.values()) + sum(v[0].size for v in src_fwd.values()),
            )
            del it
        zero = np.zeros_like(e_in)
        for p in keep:
            planes[p] = (forward_keep.get(p, zero), backward_keep.get(p, zero))
        return Result3D(transmitted, reflected, planes, stored)


def memory_estimate(nx: int, ny: int, interface_pixels: int, pairs: int, kinds: int = 4) -> int:
    """Rough peak memory in bytes: cached operators, a few working planes, sparse sources.

    Operators: kinds × pairs × nx·ny 2×2 complex matrices (64 B per pixel each); working fields:
    about 8 planes of nx·ny 2-vectors (32 B per pixel); sources: 2 × interface_pixels × (32 B
    values + 8 B index).
    """
    plane = nx * ny
    return int(kinds * pairs * plane * 64 + 8 * plane * 32 + 2 * interface_pixels * 40)
