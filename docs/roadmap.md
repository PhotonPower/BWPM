# Plan und Stufen

> Maßgebliche Fassung des Plans. Änderungen per PR; Entscheidungen im Abschnitt „Entscheidungen“.

## Ziel

Eine schnelle bidirektionale vektorielle Wellenpropagation (BWPM) nach wende2023, die

1. als Python-Referenz gegen unabhängige Lösungen geprüft ist,
2. danach als C++/CUDA-Kern mit Python-Anbindung große 2D- und 3D-Volumina rechnet, und
3. über Gitterfeld-Ports mit Raytatouille und Field Tracing gekoppelt werden kann
   („Raytatouille Fusion“, Fieldtracing `docs/rtt-fusion/plan.md`, Abschnitt 6.6).

Leitlinien:

- Konventionen von Raytatouille von Anfang an (`docs/konventionen.md`).
- Jede Formel mit Herleitung oder Beleg, jeder Test mit unabhängigem Sollwert, bei Näherungen ein
  Konvergenztest.
- Erst 2D, dann 3D; erst Python, dann C++/CUDA, mit Sollwertdateien als Brücke.

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **B0 Spezifikation** | Formeln von FPWPM und BWPM aus wende2023 in den Konventionen von Raytatouille, Konventionsabgleich der Grenzflächenmatrizen, offene Punkte | Spezifikation freigegeben; K1 grün |
| **B1 FPWPM 2D** | Vorwärtspropagation in (y, z) mit diskreten Brechzahlen, Gewichtsfunktion, absorbierende Maske (E2), Brechzahlstufen (E4) | K2, K3; Konvergenz in dz und dy; Randzone am Gaußstrahl; Glaskeil 0–60° |
| **B2 BWPM 2D** | Reflexionen, Rückwärtspropagation, Iterationen, Ausbreitungsfaktoren nach Tabelle 1 (E3) | K4 gegen Transfer-Matrix des diskretisierten Stapels auf Maschinengenauigkeit und gegen Airy; K5 |
| **B3 3D** | (x, y, z) mit Speicherung der Reflexe nur an Grenzflächen; Speicher- und Laufzeitmodell | K10–K12; Messlauf auf dem Rechner des Maintainers (`docs/testfaelle.md`, Messlauf B3) ✅ |
| **B4 Kern und Ports** | C++/CUDA-Kern, Python-Anbindung, Adapter nach dem Fusionsplan (Abschnitt 6.4), Konformitätstests (6.5); Entwurf und Teilstufen B4a–B4e in [docs/b4-entwurf.md](b4-entwurf.md) | Sollwertdateien aus B1–B3 (`reference/`) bitgleich bzw. in Toleranz |
| **B5 Fusionspilot** | Faserfacette mit Index-Matching-Kappe zwischen Raytatouille-Objektiv und Lochblende (Fusionsplan K4) | K6 |

Reihenfolge streng: B0 → B1 → B2 → B3 → B4 → B5. Ausnahme (E6): Der Fusionspilot wurde in 2D
vor B3 gerechnet ([docs/pilot-imc-2d.md](pilot-imc-2d.md)); B5 bleibt der Pilot in 3D mit dem
Kern aus B4.

## Entscheidungen

| Nr. | Frage | Entscheidung |
| --- | --- | --- |
| E1 | Eigenes Repository | ja, vom Maintainer am 2026-11-05 angelegt (Fusionsplan D11) |
| E2 | Randbedingung quer (#1) | absorbierende Maske nach jedem Schritt; großes Fenster mit Nullrand als Testreferenz (2026-11-05) |
| E3 | Rückwärtslauf: Zuordnung r, s und Lage der Felder (#2) | unabhängig über die Transfer-Matrix des diskretisierten Stapels herleiten, Tabelle 1 daran prüfen (2026-11-05) |
| E4 | Brechzahlstufen aus einer Geometrie (#3) | Pixelmittelpunkt als Standard, Flächenanteil als Option für dielektrische Übergänge (2026-11-05) |
| E5 | Grenzen aus B1/B2 (#8) | Gültigkeitsgrenze festhalten: Reflexe quantitativ nur an z-senkrechten Flächen, an gekippten Flächen nur Richtung, Fokuslage und Phase; gefangene Wellen konvergieren nicht. Gegenvergleich mit hp-FEM später, Abhilfe nur bei Bedarf (2026-11-05) |
| E6 | Fusionspilot vor B3 | in 2D vorgezogen, mit `fieldtracing` als optionaler Abhängigkeit (Extra `fusion`, Commit festgelegt) (2026-11-05) |
| E7 | Aufbau des C++/CUDA-Kerns (#16) | D1–D8 aus `docs/b4-entwurf.md`: ein Paket mit optionalem `bwpm._core`, CMake-Presets und scikit-build-core, pybind11, pocketfft/cuFFT, complex128, T/R im GPU-Kernel, OpenMP, CUDA ≥ 12.6 sm_89 (2026-11-05) |

## Risiken

- **Speicher:** Die Quelle speichert die reflektierten Felder im ganzen Volumen (N_x N_y N_z
  Vektoren) und bleibt deshalb in 2D. B3 klärt, ob die Speicherung nur an Grenzflächen reicht.
- **Randbedingungen quer:** Die FFT ist periodisch; entschieden ist eine absorbierende Maske (E2).
  Restreflexion durch zu steile Masken misst B1.
- **Leistung an gekippten Flächen (gemessen in B1/B2, K8, K9):** Transmission auf einige Prozent,
  Reflexion an gekippten Flächen grob falsch und nicht konvergent; als Gültigkeitsgrenze festgehalten (E5).
- **Geführte Wellen (K4):** Die Iteration konvergiert nicht für durch Totalreflexion gefangene
  Wellen.
- **Wahl der Brechzahlstufen:** M wächst bei gekrümmten Flächen und Gradienten; der Aufwand geht mit
  M². B1 misst den Einfluss.
