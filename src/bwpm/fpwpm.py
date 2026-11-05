"""Forward step of the FPWPM in 2D (y, z; the system is infinitely extended in x).

Specification: docs/spezifikation.md, sections 2 and 4; source wende2023, eqs. (1), (4), (5) and
table 1, case i = 1.

Discretisation: the volume consists of `nz` slabs of thickness dz along z, each sampled at `ny`
points in y (centred grid y_n = (n − ny//2)·dy, mm). Every pixel carries a material label into a
list of complex refractive indices. A field always sits at the right edge of its slab
(wende2023, table 1): the input field at z₀ lies at the right edge of a virtual slab −1 in the
medium `index_before`; after step j it lies at z₀ + (j + 1)·dz at the right edge of slab j.

Fields are the transverse components (E_x, E_y) on the grid, shape (ny, 2). In 2D with k_x = 0,
E_x is TE (s) and E_y is the transverse part of TM (p); E_z follows from eq. (5).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from bwpm._types import ComplexArray, FloatArray, IntArray
from bwpm.interface import longitudinal_wavenumber, transmission_matrix


@dataclass(frozen=True)
class Grid2D:
    """Transverse grid and step size.

    Attributes:
        ny: number of samples in y.
        dy: sample spacing in y (mm).
        dz: slab thickness (mm).
        wavelength_um: vacuum wavelength (µm).
    """

    ny: int
    dy: float
    dz: float
    wavelength_um: float

    def __post_init__(self) -> None:
        if self.ny < 2 or self.dy <= 0 or self.dz <= 0 or self.wavelength_um <= 0:
            raise ValueError("ny >= 2 and positive dy, dz, wavelength_um are required")

    @property
    def k0(self) -> float:
        """Vacuum wavenumber 2π/λ in 1/mm."""
        return float(2.0 * np.pi / (self.wavelength_um * 1e-3))

    @property
    def y(self) -> FloatArray:
        """Sample positions y_n = (n − ny//2)·dy in mm."""
        return np.asarray((np.arange(self.ny) - self.ny // 2) * self.dy, dtype=np.float64)

    @property
    def ky(self) -> FloatArray:
        """Transverse frequencies of the FFT in 1/mm (numpy order)."""
        return np.asarray(2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.dy), dtype=np.float64)


def absorbing_mask(grid: Grid2D, width: float, order: int = 4, edge: float = 1e-6) -> FloatArray:
    """Absorbing mask m(y), applied after every step (decision E2).

    m = 1 for |y| ≤ y_i and m = edge^(u^order) in the boundary zone, u = (|y| − y_i)/width ∈ [0, 1],
    y_i = y_max − width (super-Gaussian roll-off reaching `edge` at the window border).

    Args:
        grid: transverse grid.
        width: width of the boundary zone on each side (mm); 0 disables the mask.
        order: exponent of the roll-off (smoothness).
        edge: value of the mask at the outermost sample.
    """
    if width < 0 or order < 1 or not 0 < edge < 1:
        raise ValueError("width >= 0, order >= 1 and 0 < edge < 1 are required")
    y = np.abs(grid.y)
    if width == 0:
        return np.ones_like(y)
    inner = float(np.max(y)) - width
    u = np.clip((y - inner) / width, 0.0, 1.0)
    return np.asarray(edge ** (u**order), dtype=np.float64)


@dataclass
class Propagator2D:
    """Forward FPWPM through a labelled volume.

    Attributes:
        grid: transverse grid and step.
        indices: complex refractive index per material label.
        labels: material label per pixel, shape (nz, ny).
        index_before: index of the medium in which the input field lies.
        mask: absorbing mask applied after every step (shape (ny,)), or None.
    """

    grid: Grid2D
    indices: list[complex]
    labels: IntArray
    index_before: complex
    mask: FloatArray | None = None
    _cache: dict[tuple[complex, complex], ComplexArray] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        labels = np.asarray(self.labels)
        if labels.ndim != 2 or labels.shape[1] != self.grid.ny:
            raise ValueError("labels must have shape (nz, ny)")
        if labels.min() < 0 or labels.max() >= len(self.indices):
            raise ValueError("labels must index into `indices`")

    def _operator(self, index_r: complex, index_s: complex) -> ComplexArray:
        """T_r,s(k_y)·𝒫_s with 𝒫_s = exp(i k_z,s dz) (table 1, i = 1), shape (ny, 2, 2)."""
        key = (complex(index_r), complex(index_s))
        if key not in self._cache:
            g = self.grid
            ky = g.ky
            kx = np.zeros_like(ky)
            t = transmission_matrix(kx, ky, index_r, index_s, g.k0)
            phase = np.exp(1j * longitudinal_wavenumber(index_s, g.k0, ky * ky) * g.dz)
            self._cache[key] = np.asarray(t * phase[:, None, None])
        return self._cache[key]

    def step(
        self, field_in: ComplexArray, labels_before: IntArray, labels_after: IntArray
    ) -> ComplexArray:
        """One step z → z + dz, wende2023 eq. (1): Σ_r,s I_r,s·F⁻¹{T_r,s 𝒫_s F{E}}."""
        spectrum = np.fft.fft(field_in, axis=0)
        out = np.zeros_like(spectrum)
        pairs = np.unique(np.stack([labels_before, labels_after], axis=1), axis=0)
        for r, s in pairs:
            weight = (labels_before == r) & (labels_after == s)
            op = self._operator(self._index(int(r)), self.indices[int(s)])
            out[weight] = np.fft.ifft(np.einsum("nij,nj->ni", op, spectrum), axis=0)[weight]
        if self.mask is not None:
            out = out * self.mask[:, None]
        return np.asarray(out)

    def _index(self, label: int) -> complex:
        return self.index_before if label < 0 else self.indices[label]

    def run(self, field_in: ComplexArray, keep: bool = False) -> ComplexArray:
        """Propagate the input field (ny, 2) through all slabs.

        Returns:
            the field at the right edge of the last slab, shape (ny, 2); with `keep=True` all
            planes, shape (nz + 1, ny, 2), plane 0 being the input.
        """
        field_now: ComplexArray = np.asarray(field_in, dtype=np.complex128)
        if field_now.shape != (self.grid.ny, 2):
            raise ValueError("field must have shape (ny, 2)")
        before = np.full(self.grid.ny, -1)
        planes = [field_now] if keep else []
        for row in np.asarray(self.labels):
            field_now = self.step(field_now, before, row)
            before = row
            if keep:
                planes.append(field_now)
        return np.asarray(planes) if keep else field_now


def longitudinal_component(
    field_xy: ComplexArray, grid: Grid2D, index: complex, direction: int = 1
) -> ComplexArray:
    """E_z from the transverse components in a homogeneous medium, wende2023 eq. (5).

    E_z(k⊥) = −direction·k⊥·E⊥(k⊥)/k_z (direction +1 for +z, −1 for −z waves), shape (ny,).
    """
    spectrum = np.fft.fft(field_xy, axis=0)
    kz = longitudinal_wavenumber(index, grid.k0, grid.ky**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ez = np.where(kz != 0, -direction * grid.ky * spectrum[:, 1] / kz, 0.0)
    return np.asarray(np.fft.ifft(ez))


def power_flux(field_xy: ComplexArray, grid: Grid2D, index: complex, direction: int = 1) -> float:
    """Power through the plane per unit length in x, up to the constant 1/(2ωμ₀) (arbitrary units).

    Sum over plane waves (Parseval) of Re(k_z)·|E|² with E_z from eq. (5); evanescent waves carry
    no power. Defined for a homogeneous medium in the plane.
    """
    spectrum = np.fft.fft(field_xy, axis=0) / grid.ny
    kz = longitudinal_wavenumber(index, grid.k0, grid.ky**2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ez = np.where(kz != 0, -direction * grid.ky * spectrum[:, 1] / kz, 0.0)
    magnitude = np.sum(np.abs(spectrum) ** 2, axis=1) + np.abs(ez) ** 2
    return float(np.sum(np.real(kz) * magnitude) * grid.ny * grid.dy)
