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
| **B1 FPWPM 2D** | Vorwärtspropagation in (y, z) mit diskreten Brechzahlen, Gewichtsfunktion, Randbedingung quer | K2, K3; Konvergenz in dz und dy |
| **B2 BWPM 2D** | Reflexionen, Rückwärtspropagation, Iterationen, Ausbreitungsfaktoren nach Tabelle 1 | K4 (Platte gegen Airy, wie Abb. 2 der Quelle), K5 |
| **B3 3D** | (x, y, z) mit Speicherung der Reflexe nur an Grenzflächen; Speicher- und Laufzeitmodell | 3D-Platte, Linse; Messungen auf dem Rechner des Maintainers |
| **B4 Kern und Ports** | C++/CUDA-Kern, Python-Anbindung, Adapter nach dem Fusionsplan (Abschnitt 6.4), Konformitätstests (6.5) | Sollwertdateien aus B1–B3 bitgleich bzw. in Toleranz |
| **B5 Fusionspilot** | Faserfacette mit Index-Matching-Kappe zwischen Raytatouille-Objektiv und Lochblende (Fusionsplan K4) | K6 |

Reihenfolge streng: B0 → B1 → B2 → B3 → B4 → B5.

## Entscheidungen

| Nr. | Frage | Entscheidung |
| --- | --- | --- |
| E1 | Eigenes Repository | ja, vom Maintainer am 2026-11-05 angelegt (Fusionsplan D11) |

## Risiken

- **Speicher:** Die Quelle speichert die reflektierten Felder im ganzen Volumen (N_x N_y N_z
  Vektoren) und bleibt deshalb in 2D. B3 klärt, ob die Speicherung nur an Grenzflächen reicht.
- **Randbedingungen quer:** Die FFT ist periodisch; die Quelle nennt „perfectly transparent
  boundary conditions“ ohne Angabe des Verfahrens. B1 entscheidet (Frage im Issue).
- **Wahl der Brechzahlstufen:** M wächst bei gekrümmten Flächen und Gradienten; der Aufwand geht mit
  M². B1 misst den Einfluss.
