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
| **B4 Kern und Ports** | C++/CUDA-Kern, Python-Anbindung, Adapter nach dem Fusionsplan (Abschnitt 6.4), Konformitätstests (6.5); Entwurf und Teilstufen B4a–B4e in [docs/b4-entwurf.md](b4-entwurf.md) | Sollwertdateien aus B1–B3 (`reference/`) bitgleich bzw. in Toleranz ✅ (B4a–B4e, K13–K16) |
| **B5 Fusionspilot** | Faserfacette mit Index-Matching-Kappe zwischen Raytatouille-Objektiv und Lochblende (Fusionsplan K4) | K6 |
| **B6 Große Winkel (#29)** | Genauigkeit des Vorwärtsschritts bei Streuung in große Winkel an schwachen, feinen Brechzahländerungen; Sollwert erste Born-Näherung (E8) | K17 ✅ (Ursache erklärt, Grenze E9) |
| **B7 Schwache Fluktuationen (#37, E10)** | Zerlegung ñ = ñ_Struktur + δn: Struktur weiter mit WPM, δn je Scheibe mit dem vektoriellen Born-Inkrement im Wirtsmedium; B7a Python-Referenz, B7b Kern (C++/CUDA) | B7a ✅ (abgenommen 2026-11-19): K18 (K17-Sollwerte ohne Anpassung getroffen, Konvergenz in dz; Mehrfachstreuung gegen RCWA; Energie mit E_z; δn = 0 bitgleich, K1–K16 und `reference/` unberührt). B7b: neue Sollwertdateien aus B7a; CPU-Kern ✅ (K18 Punkt 7), CUDA ✅ auf RTX 2080 (K18 Punkt 8); B7b ✅ abgenommen 2026-11-19 |

Reihenfolge streng: B0 → B1 → B2 → B3 → B4 → B5. Ausnahme (E6): Der Fusionspilot wurde in 2D
vor B3 gerechnet ([docs/pilot-imc-2d.md](pilot-imc-2d.md)); B5 bleibt der Pilot in 3D mit dem
Kern aus B4. B5a (3D ohne Raytatouille) ist gemessen, siehe [docs/pilot-imc-3d.md](pilot-imc-3d.md);
B5b (Kopplung mit Raytatouille) wartet auf Raytatouille v1.0 und den Port. B6 (E8) hängt nicht
von B5b ab und lief vorher; abgeschlossen mit der Gültigkeitsgrenze E9. B7 (E10) hängt ebenfalls
nicht von B5b ab und läuft, solange B5b wartet; B7a ist abgenommen (Python-Referenz in 2D und
3D, K18). B7b ist abgenommen (2026-11-19): CPU- und CUDA-Kern, geprüft auf einer RTX 2080 mit
CUDA 13.4 und 13.5; ein Lauf auf sm_89 (RTX 6000 Ada) steht noch aus, ist aber keine Bedingung der
Abnahme. B7 ist damit abgeschlossen.

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
| E8 | Befund aus Scattuccino S9b-2 (#29) | als Stufe B6 vor B5b prüfen; Sollwert erste Born-Näherung einer dünnen Schicht, selbst gegen Airy und RCWA geprüft (2026-11-19) |
| E9 | Grenze bei Streuung in große Winkel (#29, K17) | Gültigkeitsgrenze festhalten, das Verfahren bleibt wie wende2023: Streuung an feinen Brechzahländerungen quantitativ (Amplitude < 1 %) nur bis θ ≈ 11° (s) bzw. 6,6° (p); darüber s-Anteil um r = cos θ (3 − cos θ)/2 zu klein, p-Anteil um r/cos²θ zu groß. Erweiterung für schwache Fluktuationen (Variante B3: Born-konsistentes Inkrement für δn) als Vorschlag in #37, nicht beschlossen (2026-11-19) |
| E10 | Erweiterung für schwache Fluktuationen (#37) | Variante B3 angenommen als Stufe B7: δn je Scheibe mit dem Born-Inkrement F⁻¹{e^(i k_z dz)·i dz/(2k_z)·P⊥·F{V_δn·E}} im homogenen Wirtsmedium, Struktur unverändert mit WPM; mit δn = 0 bitgleich zu heute. Erste Ausbaustufe: Wirtsmedium je Bereich homogen, nur Vorwärtsrichtung, streifende Wellen regularisiert oder auf propagierende beschränkt (in B7a festzulegen); Rückstreuung und δn in Linsen später. E9 gilt weiter für die reine WPM (2026-11-19) |
| E11 | Potenzial im Born-Inkrement (#46) | linearisiert, V = 2ñ₀k₀²δn statt k₀²((ñ₀ + δn)² − ñ₀²) nach chen2020 Gl. (1): mittlere Phase bei senkrechter Ausbreitung bis zur zweiten Ordnung exakt, Fehler gegen RCWA etwa ein Drittel; in erster Ordnung gleich (2026-11-19) |

## Risiken

- **Speicher:** Die Quelle speichert die reflektierten Felder im ganzen Volumen (N_x N_y N_z
  Vektoren) und bleibt deshalb in 2D. B3 klärt, ob die Speicherung nur an Grenzflächen reicht.
- **Randbedingungen quer:** Die FFT ist periodisch; entschieden ist eine absorbierende Maske (E2).
  Restreflexion durch zu steile Masken misst B1.
- **Leistung an gekippten Flächen (gemessen in B1/B2, K8, K9):** Transmission auf einige Prozent,
  Reflexion an gekippten Flächen grob falsch und nicht konvergent; als Gültigkeitsgrenze festgehalten (E5).
- **Streuung in große Winkel (K17, E9):** An feinen Brechzahländerungen streut der Schritt den
  s-Anteil zu schwach und den p-Anteil zu stark (Faktor 0,80 bzw. 1,68 in der Amplitude bei 46°);
  quantitativ nur bis etwa 6,6° (p) bzw. 11° (s). Betrifft Streumedien (Scattuccino), nicht die
  Brechung an glatten Flächen. Abhilfe für schwache Fluktuationen: B7 (E10).
- **Geführte Wellen (K4):** Die Iteration konvergiert nicht für durch Totalreflexion gefangene
  Wellen.
- **Wahl der Brechzahlstufen:** M wächst bei gekrümmten Flächen und Gradienten; der Aufwand geht mit
  M². B1 misst den Einfluss.
