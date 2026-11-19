"""Vector against scalar Born for the Whittle–Matérn media of Scattuccino S9b-2 (#29, K17).

The scalar Born cross-section of a random index field is σ(θ) ∝ k⁴ Φ(2k sin(θ/2)) with the
Whittle–Matérn spectrum Φ(q) ∝ (1 + q² l_c²)^(−D/2) (constants cancel in the ratios below). The
vector Born cross-section for linearly polarised light carries the dipole factor
1 − sin²θ cos²φ, azimuthally averaged (1 + cos²θ)/2 (chew1995, dyadic Green's function).
For the forward FPWPM, K17 (docs/spezifikation.md §9) predicts the transverse field r(θ)·S with
r = cos θ (3 − cos θ)/2, so the measure of S9b-2 (|E⊥|² against scalar Born) is weighted by r².

Prints, per medium: 1 − g over the forward hemisphere (scalar Born, vector Born, K17 prediction
for the S9b-2 measure) and μs(vector)/μs(scalar) over the full sphere.

Run: python experiments/vector_born_wm.py
"""

from __future__ import annotations

import numpy as np

N0 = 1.38
WAVELENGTH_UM = 0.8
MEDIA = [(0.5, 4.0), (0.2, 4.0), (1.0, 4.0), (2.0, 4.0), (0.5, 3.5), (0.5, 5.0)]  # (l_c µm, D)


def one_minus_g(mu: np.ndarray, weight: np.ndarray) -> float:  # type: ignore[type-arg]
    fwd = mu > 0
    g = np.trapezoid(mu[fwd] * weight[fwd], mu[fwd]) / np.trapezoid(weight[fwd], mu[fwd])
    return float(1 - g)


def main() -> None:
    k = 2 * np.pi * N0 / WAVELENGTH_UM  # 1/µm
    mu = np.linspace(-1.0, 1.0, 400_001)
    q = 2 * k * np.sqrt((1 - mu) / 2)
    dipole = (1 + mu**2) / 2
    r2 = np.where(mu > 0, (mu * (3 - mu) / 2) ** 2, 0.0)
    print(
        "| l_c (µm) | D | 1 − g Born skalar | 1 − g Born vektoriell | 1 − g K17 (S9b-2-Messgröße) "
        "| μs vektoriell/skalar |"
    )
    print("| --- | --- | --- | --- | --- | --- |")
    for l_c, d in MEDIA:
        psd = (1 + (q * l_c) ** 2) ** (-d / 2)
        ratio = np.trapezoid(psd * dipole, mu) / np.trapezoid(psd, mu)
        print(
            f"| {l_c} | {d} | {one_minus_g(mu, psd):.4f} | {one_minus_g(mu, psd * dipole):.4f} "
            f"| {one_minus_g(mu, psd * r2):.4f} | {ratio:.3f} |"
        )


if __name__ == "__main__":
    main()
