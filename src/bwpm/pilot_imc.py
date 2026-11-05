"""Fusion pilot in 2D: imaging-fibre facet with a 3D-printed index-matching cap (IMC).

Setup after wende2023 §2.2.2: λ = 0.635 µm, TE, dry objective NA 0.4, fibre core n = 1.5 with
diameter 2 µm and NA_core = 0.39 (cladding n = sqrt(1.5² − 0.39²)), IMC of IP-Q (n = 1.511) and
50 µm thickness on the facet. Two configurations: bare fibre (facet in the focus of the
objective) and fibre with IMC (facet at the best focus behind the cap).

Chain (Raytatouille Fusion, Fieldtracing docs/rtt-fusion/plan.md, pilot K4):

1. **Field Tracing:** ideal aplanatic cylindrical objective (focal length F, pupil amplitude cos θ
   on the plane after the lens, spherical phase). Its pointwise Fourier transform
   (`fieldtracing.pft.pft_1d`) gives the plane-wave spectrum at the BWPM grid frequencies; the
   spectral phase is propagated analytically to the input plane of the BWPM window.
2. **BWPM:** forward and backward fields through air, IMC and fibre (`bwpm.bidirectional`).
3. **Field Tracing:** the reflected field is propagated with the exact angular spectrum
   (`fieldtracing.propagation.angular_spectrum`) to find its (virtual or real) focus.

Requires the optional package `fieldtracing` (extra `fusion`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from bwpm._types import ComplexArray, FloatArray
from bwpm.bidirectional import Bidirectional2D
from bwpm.fpwpm import Grid2D, Propagator2D, absorbing_mask, power_flux

WAVELENGTH_UM = 0.635
NA_OBJECTIVE = 0.4
N_CORE = 1.5
NA_CORE = 0.39
N_CLAD = float(np.sqrt(N_CORE**2 - NA_CORE**2))
N_IMC = 1.511
CORE_DIAMETER = 2e-3  # mm
IMC_THICKNESS = 50e-3  # mm
FOCAL_LENGTH = 10e-3  # mm, cylindrical objective; only the spectrum within NA matters


@dataclass(frozen=True)
class PilotResult:
    """Results of one configuration (lengths in mm, powers relative to the incident power).

    Attributes:
        reflected_power: power of the backward field at the input plane.
        backward_profile: backward power at every plane (nz + 1).
        z_planes: positions of the planes (mm), facet at z = 0.
        reflected_focus: position of the focus of the reflected light (mm, facet at 0).
        forward_spectrum_in_imc: (k_y, spectrum) of the forward field 1 µm before the facet
            (IMC case) or None.
    """

    reflected_power: float
    backward_profile: FloatArray
    z_planes: FloatArray
    reflected_focus: float
    forward_spectrum_in_imc: tuple[FloatArray, ComplexArray] | None


def make_grid(ny: int = 1152) -> Grid2D:
    """λ/(5 n_max) sampling as in wende2023 §2.2.2 (dy = dz ≈ 0.084 µm)."""
    d = WAVELENGTH_UM * 1e-3 / (5 * N_IMC)
    return Grid2D(ny=ny, dy=d, dz=d, wavelength_um=WAVELENGTH_UM)


def objective_field(grid: Grid2D, focus_distance: float) -> ComplexArray:
    """TE field (ny, 2) at the input plane of a beam that focuses `focus_distance` behind it in air.

    Field Tracing step: the exit-pupil field U·e^(iψ), ψ = −k(sqrt(F² + y²) − F), U = cos θ, is
    transformed with the PFT (`fieldtracing.pft.pft_1d`) at the stationary points
    y = −F k_y/k_z of the grid frequencies |k_y| ≤ k·NA. Its spectral phase k F − k_z F is replaced
    by the exact propagation to the input plane, −k_z·focus_distance (both analytic). U = cos θ
    gives the aplanatic plane-wave amplitude ∝ sqrt(cos θ) per angle.
    """
    from fieldtracing.pft import pft_1d

    k = grid.k0
    ky = grid.ky
    inside = np.abs(ky) <= k * NA_OBJECTIVE
    kz = np.sqrt(k * k - ky[inside] ** 2)
    f = FOCAL_LENGTH
    y = -f * ky[inside] / kz
    root = np.sqrt(f * f + y * y)
    phase = -k * (root - f)
    gradient = -k * y / root
    curvature = -k * f * f / root**3
    amplitude = (f / root).astype(np.complex128)
    _, residual, _ = pft_1d(y, amplitude, phase, gradient, curvature)
    spectrum = np.zeros(grid.ny, dtype=np.complex128)
    spectrum[inside] = residual * np.exp(-1j * kz * focus_distance)
    field = np.zeros((grid.ny, 2), dtype=np.complex128)
    # E(y_n) = Σ Ẽ(k_y) e^(i k_y y_n) on the centred grid y_n = (n − ny//2)·dy (Δk dropped)
    shift = np.exp(-1j * ky * (grid.ny // 2) * grid.dy)
    field[:, 0] = np.fft.ifft(spectrum * shift) * grid.ny
    return field


def _volume(
    grid: Grid2D, z_in: float, z_end: float, with_imc: bool
) -> tuple[list[complex], np.ndarray]:  # type: ignore[type-arg]
    """Labels (nz, ny): 0 air, 1 IMC, 2 core, 3 cladding; facet at z = 0."""
    nz = round((z_end - z_in) / grid.dz)
    z = z_in + (np.arange(nz) + 0.5) * grid.dz
    zz, yy = np.meshgrid(z, grid.y, indexing="ij")
    labels = np.zeros((nz, grid.ny), dtype=np.int64)
    if with_imc:
        labels[(zz > -IMC_THICKNESS) & (zz <= 0)] = 1
    fibre = zz > 0
    labels[fibre] = 3
    labels[fibre & (np.abs(yy) <= CORE_DIAMETER / 2)] = 2
    return [1.0 + 0j, N_IMC + 0j, N_CORE + 0j, N_CLAD + 0j], labels


def run_configuration(
    with_imc: bool, air_focus_depth: float, grid: Grid2D | None = None, iterations: int = 2
) -> PilotResult:
    """Run one configuration.

    Args:
        with_imc: cap on the facet or bare fibre.
        air_focus_depth: depth of the objective's focus in air behind the first surface (IMC
            front or bare facet), mm. The input plane lies 20 µm in front of that surface.
        grid: transverse grid (default `make_grid()`).
        iterations: BWPM iterations.
    """
    from fieldtracing.propagation import angular_spectrum

    g = grid or make_grid()
    front = -IMC_THICKNESS if with_imc else 0.0
    z_in = front - 20e-3
    z_end = 15e-3
    indices, labels = _volume(g, z_in, z_end, with_imc)
    field = objective_field(g, focus_distance=20e-3 + air_focus_depth)
    mask = absorbing_mask(g, 20 * WAVELENGTH_UM * 1e-3, 8)
    prop = Propagator2D(g, indices, labels, 1.0 + 0j, mask)
    res = Bidirectional2D(prop).run(field, iterations=iterations, keep=True)
    assert res.forward is not None and res.backward is not None
    p_in = power_flux(field, g, 1.0)
    z_planes = z_in + np.arange(labels.shape[0] + 1) * g.dz
    centre_label = np.concatenate([[0], labels[:, g.ny // 2]])
    backward = np.asarray(res.backward)
    profile = np.array(
        [
            power_flux(backward[i], g, indices[int(centre_label[i])], direction=-1) / p_in
            for i in range(backward.shape[0])
        ]
    )
    # focus of the reflected light: propagate the backward field from the input plane in its own
    # direction (−z, mirrored frame: positive steps) and search the peak intensity
    k = g.k0
    reflected = np.fft.fft(res.reflected[:, 0])
    reflected[np.abs(g.ky) >= k] = 0.0  # evanescent parts would grow for negative steps
    reflected = np.fft.ifft(reflected)
    steps = np.arange(-160e-3, 40e-3, 0.25e-3)

    def propagated(step: float) -> ComplexArray:
        # negative steps by conjugation: conj(AS(conj E, |s|)) = AS(E, s) for propagating waves
        if step >= 0:
            return np.asarray(angular_spectrum(reflected, g.dy, k, step))
        return np.asarray(np.conj(angular_spectrum(np.conj(reflected), g.dy, k, -step)))

    peaks = [float(np.max(np.abs(propagated(float(s))) ** 2)) for s in steps]
    j = int(np.argmax(peaks))
    a, b, c = peaks[max(j - 1, 0)], peaks[j], peaks[min(j + 1, len(peaks) - 1)]
    denominator = a - 2 * b + c
    refined = steps[j] + (0.5 * (a - c) / denominator * 0.25e-3 if denominator != 0 else 0.0)
    reflected_focus = z_in - refined  # a positive step moves towards −z
    spectrum = None
    if with_imc:
        plane = int(np.argmin(np.abs(z_planes - (-1e-3))))
        centred = np.fft.fft(res.forward[plane, :, 0]) * np.exp(1j * g.ky * (g.ny // 2) * g.dy)
        spectrum = (g.ky, centred)
    return PilotResult(
        reflected_power=float(profile[0]),
        backward_profile=profile,
        z_planes=z_planes,
        reflected_focus=float(reflected_focus),
        forward_spectrum_in_imc=spectrum,
    )


def strehl_from_spectrum(
    ky: FloatArray, spectrum: ComplexArray, k: float, offsets: FloatArray
) -> tuple[float, float]:
    """2D Strehl ratio max_z |Σ A e^(i k_z z)|² / (Σ |A|)² and the best focus offset (mm).

    The spectrum A(k_y) is referred to y = 0 (centred grid) and taken in a homogeneous medium with
    wavenumber k; offsets are distances from that plane along +z. The Strehl ratio is evaluated on
    the axis y = 0.
    """
    kz = np.sqrt(np.maximum(k * k - ky * ky, 0.0))
    norm = float(np.sum(np.abs(spectrum))) ** 2
    values = [float(np.abs(np.sum(spectrum * np.exp(1j * kz * dz))) ** 2) / norm for dz in offsets]
    j = int(np.argmax(values))
    return values[j], float(offsets[j])
