# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt

- Projektgerüst: Regeln für Agenten, Plan (`docs/roadmap.md`), Konventionen wie Raytatouille,
  Quellen, Testfälle, CI (ruff, mypy, pytest unter Linux und Windows).
- Spezifikation B0 (`docs/spezifikation.md`): Formeln von FPWPM und BWPM aus wende2023 in den
  Konventionen von Raytatouille, offene Punkte.
- `bwpm.interface`: Fresnel-Koeffizienten (Convention A), Transmissions- und Reflexionsmatrix nach
  wende2023 Gl. (2), (3), (8), (9); Nachweis, dass sie mit Convention A exakt der vektorielle
  Fresnel-Operator sind (K1).
- Entscheidungen E2–E4 (Randmaske, Prüfweg für den Rückwärtslauf, Brechzahlstufen) in Plan
  und Spezifikation; Testfälle K7, K8.
- B1: Vorwärtspropagation der FPWPM in 2D (`bwpm.fpwpm`: `Grid2D`, `Propagator2D`,
  `absorbing_mask`, `power_flux`, `longitudinal_component`) und Brechzahlstufen aus einer Geometrie
  (`bwpm.geometry`: Mittelpunkt, Flächenanteil). Testfälle K2, K3 (exakt), K7 (Regel für die
  Randzone), K8 (Gültigkeit an gekippten Flächen: Richtung bis 45°, Leistung nur auf einige Prozent).
