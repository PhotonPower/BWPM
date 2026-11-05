"""Bidirectional wave propagation (BWPM) in 2D: reflections, backward pass, iterations.

Specification: docs/spezifikation.md, section 5; source wende2023, eqs. (6), (7), (10), (11) and
table 1. Assignment of the media and positions derived independently (decision E3):

Planes p = 0 … nz lie at z₀ + p·dz. Every field sits at the right edge of its slab, so plane p
belongs to slab p − 1 (plane 0 to the medium `index_before`). Interface j (at plane j) separates
slab j − 1 (left, medium r) and slab j (right). Behind the last slab its medium continues
(no reflection at plane nz).

- Forward step p → p + 1 (table 1, i = 1): crosses interface p from r = medium(plane p) into
  s = slab p, then runs dz in s: T_r,s·exp(i k_z,s dz).
- Forward reflection at interface p (i = 2): R_r,s, factor 1; the reflected field stays at plane p.
- Backward step p → p − 1 (i = 3): first runs dz through its own slab p − 1 (medium r), then
  crosses interface p − 1 into slab p − 2 (medium s): T_r,s·exp(i k_z,r dz).
- Backward reflection at interface p − 1 (i = 4): back through slab p − 1 to plane p:
  R_r,s·exp(2i k_z,r dz).

One iteration is a forward pass (input plus the sources from the last backward pass) and a
backward pass (sources from this forward pass). For planar layers and plane waves the iterations
converge to the exact multilayer solution of the discretised stack (test K4).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from bwpm._types import ComplexArray, IntArray
from bwpm.fpwpm import Propagator2D, power_flux
from bwpm.interface import longitudinal_wavenumber, reflection_matrix, transmission_matrix

Operator = Callable[[complex, complex], ComplexArray]


@dataclass(frozen=True)
class BidirectionalResult:
    """Fields after the last iteration.

    Attributes:
        transmitted: forward field at plane nz, shape (ny, 2), in the medium of the last slab.
        reflected: backward field at plane 0, shape (ny, 2), in the medium `index_before`.
        forward, backward: all planes (nz + 1, ny, 2) if requested, else None.
        energy_residual: ε = 1 − (P_R + P_T)/P_in after every iteration (wende2023 §2.1.3);
            meaningful for lossless volumes only.
    """

    transmitted: ComplexArray
    reflected: ComplexArray
    forward: ComplexArray | None
    backward: ComplexArray | None
    energy_residual: list[float]


class Bidirectional2D:
    """BWPM on top of the forward propagator (grid, indices, labels, mask)."""

    def __init__(self, propagator: Propagator2D) -> None:
        self.p = propagator
        self._cache: dict[tuple[str, complex, complex], ComplexArray] = {}
        labels = np.asarray(propagator.labels)
        before = np.full((1, propagator.grid.ny), -1, dtype=np.int64)
        # medium label of every plane: plane 0 -> before (-1), plane p -> slab p - 1
        self.plane_labels: IntArray = np.concatenate([before, labels], axis=0)

    # --- operators of table 1 ------------------------------------------------------------
    def _kz(self, index: complex) -> ComplexArray:
        g = self.p.grid
        return longitudinal_wavenumber(index, g.k0, g.ky**2)

    def _matrix(self, kind: str, r: complex, s: complex) -> ComplexArray:
        key = (kind, complex(r), complex(s))
        if key not in self._cache:
            g = self.p.grid
            kx = np.zeros_like(g.ky)
            if kind == "fwd_t":  # i = 1
                m = (
                    transmission_matrix(kx, g.ky, r, s, g.k0)
                    * np.exp(1j * self._kz(s) * g.dz)[:, None, None]
                )
            elif kind == "fwd_r":  # i = 2
                m = reflection_matrix(kx, g.ky, r, s, g.k0)
            elif kind == "bwd_t":  # i = 3
                m = (
                    transmission_matrix(kx, g.ky, r, s, g.k0)
                    * np.exp(1j * self._kz(r) * g.dz)[:, None, None]
                )
            else:  # "bwd_r", i = 4
                m = (
                    reflection_matrix(kx, g.ky, r, s, g.k0)
                    * np.exp(2j * self._kz(r) * g.dz)[:, None, None]
                )
            self._cache[key] = np.asarray(m)
        return self._cache[key]

    def _index(self, label: int) -> complex:
        return self.p.index_before if label < 0 else self.p.indices[label]

    def _apply(
        self, kind: str, field: ComplexArray, labels_r: IntArray, labels_s: IntArray
    ) -> ComplexArray:
        """Σ_r,s I_r,s·F⁻¹{M_r,s F{E}} for the operator kind; zero where r = s for reflections."""
        spectrum = np.fft.fft(field, axis=0)
        out = np.zeros_like(spectrum)
        pairs = np.unique(np.stack([labels_r, labels_s], axis=1), axis=0)
        for r, s in pairs:
            ir, is_ = self._index(int(r)), self._index(int(s))
            if kind.endswith("_r") and ir == is_:
                continue
            weight = (labels_r == r) & (labels_s == s)
            m = self._matrix(kind, ir, is_)
            out[weight] = np.fft.ifft(np.einsum("nij,nj->ni", m, spectrum), axis=0)[weight]
        if self.p.mask is not None and kind.endswith("_t"):
            out = out * self.p.mask[:, None]
        return np.asarray(out)

    # --- passes ----------------------------------------------------------------------------
    def run(
        self, field_in: ComplexArray, iterations: int = 3, keep: bool = False
    ) -> BidirectionalResult:
        """Run `iterations` forward-backward cycles for the input field (ny, 2) at plane 0."""
        if iterations < 1:
            raise ValueError("iterations must be >= 1")
        g = self.p.grid
        nz = self.plane_labels.shape[0] - 1
        lab = self.plane_labels
        field_in = np.asarray(field_in, dtype=np.complex128)
        if field_in.shape != (g.ny, 2):
            raise ValueError("field must have shape (ny, 2)")
        source_fwd = np.zeros((nz + 1, g.ny, 2), dtype=np.complex128)  # from backward reflections
        forward = np.zeros_like(source_fwd)
        backward = np.zeros_like(source_fwd)
        p_in = power_flux(field_in, g, self.p.index_before)
        residual: list[float] = []
        for _ in range(iterations):
            source_bwd = np.zeros_like(source_fwd)  # reflections of this forward pass
            e = field_in + source_fwd[0]
            forward[0] = e
            for p in range(nz):
                source_bwd[p] = self._apply("fwd_r", e, lab[p], lab[p + 1])
                e = self._apply("fwd_t", e, lab[p], lab[p + 1]) + source_fwd[p + 1]
                forward[p + 1] = e
            source_fwd = np.zeros_like(source_fwd)
            e = source_bwd[nz]
            backward[nz] = e
            for p in range(nz, 0, -1):
                source_fwd[p] = self._apply("bwd_r", e, lab[p], lab[p - 1])
                e = self._apply("bwd_t", e, lab[p], lab[p - 1]) + source_bwd[p - 1]
                backward[p - 1] = e
            last = self._index(int(lab[nz][g.ny // 2]))
            p_r = power_flux(backward[0], g, self.p.index_before, direction=-1)
            p_t = power_flux(forward[nz], g, last)
            residual.append(1.0 - (p_r + p_t) / p_in if p_in > 0 else float("nan"))
        return BidirectionalResult(
            transmitted=forward[nz].copy(),
            reflected=backward[0].copy(),
            forward=forward if keep else None,
            backward=backward if keep else None,
            energy_residual=residual,
        )
