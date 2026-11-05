# Konventionen

Verbindlich für allen Code. Sie stimmen mit Raytatouille (`docs/architecture.md`, Abschnitt
„Konventionen“, gelesen bei Commit `6078e55`) und Fieldtracing (`docs/konventionen.md`) überein,
damit die Kopplung keine Vorzeichen ändert.

## Einheiten und Koordinaten

| Größe | Einheit | Bemerkung |
| --- | --- | --- |
| Längen, Schrittweiten dx, dy, dz | mm | wie Raytatouille; Eingaben in µm nur an der Benutzeroberfläche mit Suffix `_um` |
| Vakuumwellenlänge | µm | k₀ = 2π/λ in 1/mm |
| Wellenzahl, Ortsfrequenz k⊥ = (k_x, k_y) | 1/mm | |
| Winkel | rad | Grad nur mit Suffix `_deg` |
| Brechzahl | absolut, komplex ñ = n + iκ mit κ ≥ 0 | Luft wie Raytatouille, wenn die Zahl von dort kommt |

Koordinatensystem rechtshändig, Hauptausbreitungsrichtung +z, Querschnitt (x, y). In 2D-Rechnungen
ist das System in x unendlich ausgedehnt (∂/∂x = 0), gerechnet wird in (y, z) wie in der Quelle.

## Felder

- **Zeitkonvention:** e^(−iωt), ebene Wellen ∝ exp(i(k·r − ωt)). Ein größeres Argument heißt
  Verzögerung. Die Quelle (Wende et al. 2023) schreibt die Ausbreitung als e^(i k_z dz), also
  dieselbe Konvention.
- **Longitudinale Wellenzahl:** k_z = √(ñ²k₀² − |k⊥|²) mit Im k_z ≥ 0 und, falls Im k_z = 0,
  Re k_z ≥ 0 (Raytatouille: q = ñ cos θ). Evaneszente Anteile klingen in Ausbreitungsrichtung ab.
- **Rücklaufende Felder** (−z) werden wie in der Quelle als vorwärtslaufend in einem lokal in z
  gespiegelten System behandelt; gespeichert werden die Querkomponenten (E_x, E_y) im globalen
  System. E_z folgt aus der Transversalität, für die vorlaufende Welle E_z = −k⊥·E⊥/k_z, für die
  rücklaufende E_z = +k⊥·E⊥/k_z.
- **Fourier-Transformation** wie Fieldtracing: F{f}(k) = (2π)^−d ∫ f(x) e^(−ik·x) dᵈx, zentrierte
  Gitter x_n = (n − N//2)·Δx, Δk = 2π/(N Δx). Für Propagationsschritte ist nur das Paar F⁻¹F
  maßgeblich; die Normierung kürzt sich.

## Fresnel und Polarisation

Convention A nach Byrnes (arXiv:1603.02720v5, Gl. (6)), wie Raytatouille, mit q_j = ñ_j cos θ_j
= k_z,j/k₀:

- r_s = (q₁ − q₂)/(q₁ + q₂), t_s = 2q₁/(q₁ + q₂),
- r_p = (ñ₂ cos θ₁ − ñ₁ cos θ₂)/(ñ₂ cos θ₁ + ñ₁ cos θ₂), t_p = 2ñ₁ cos θ₁/(ñ₂ cos θ₁ + ñ₁ cos θ₂),
- bei senkrechtem Einfall r_p = −r_s.
- Basis: N zeigt ins Medium 2, s = k_ein × N/|k_ein × N|, p = k̂ × s für jede Welle.
- **TE = s, TM = p.** Die Quelle verwendet t_TE, t_TM, r_TE, r_TM; sie sind genau die Koeffizienten
  der Convention A (Nachweis in `docs/spezifikation.md`, Abschnitt 3, und
  `tests/test_interface.py`).

## Ports und Leistung

Ports zu anderen Modulen sind Gitterfelder in Ebenen z = const, getrennt nach Richtung (+z, −z),
mit Querkomponenten (E_x, E_y) als Feldamplituden in einer für alle Ports gleichen, sonst
beliebigen Einheit (keine Leistungsnormierung im Feld). Die Leistung durch eine Ebene ist der Poynting-Fluss; die Umrechnung auf Raytatouilles
leistungsnormierte Strahlen geschieht an der Port-Grenze nach Fieldtracing ADR 0004. Einzelheiten
legt die Stufe B4 fest.

## Code

- Funktionen sind rein und deterministisch; Zufall nur mit explizitem Seed.
- Fehler in der Eingabe: `ValueError` mit klarer Meldung.
- Arrays als NumPy-Arrays je Größe, letzte Achse für Vektor- bzw. Matrixkomponenten.
