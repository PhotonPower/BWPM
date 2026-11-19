# Hinweise für Scattuccino zu S9b-2 (aus BWPM K17, #29)

Vorschlag aus BWPM, nur lesend erarbeitet (AGENTS.md, Regel 7); Umsetzung und Entscheidung
liegen bei Scattuccino. Grundlage: `docs/bericht-s9b2.md`, `src/scattuccino/microstructure.py`,
`tools/calibrate_bwpm.py` und `tools/out/phase_bwpm_vs_born.csv` in Scattuccino (Stand nach
PR #25), in BWPM K17 (`docs/testfaelle.md`, `docs/spezifikation.md` Abschnitt 9).

## Befund in Kürze

Das in S9b-2 gemeldete Defizit bei großen Winkeln (1 − g −21 % bei l_c = 0,5 µm, Dichte 0,8 bei
40°) hat drei Teile:

1. **Die Referenz ist skalar.** `microstructure.born` nutzt σ(θ) = (2π k⁴/n₀²) Φ(2k sin(θ/2)) ohne
   Polarisationsfaktor. Die vektorielle erste Born-Näherung hat für linear polarisiertes Licht den
   Dipolfaktor 1 − sin²θ cos²φ, über den Azimut gemittelt (1 + cos²θ)/2 (0,79 bei 40°, 0,5 bei
   90°); für unpolarisiertes Licht gilt derselbe Mittelwert.
2. **Die Messgröße lässt E_z weg.** `angular_distribution` und `calibrate_with_bwpm` summieren
   |F{E_x}|² + |F{E_y}|². Die Leistung einer ebenen Welle braucht |E|² = |E⊥|² + |E_z|² mit
   E_z = −k⊥·E⊥/k_z (wende2023 Gl. (5), wie `bwpm.fpwpm.power_flux`); ohne E_z fehlt dem p-Anteil
   der Faktor 1/cos²θ.
3. **Die FPWPM weicht selbst ab (K17).** In erster Ordnung behandelt der Schritt E_x und E_y
   komponentenweise skalar und streut mit dem Schrägfaktor der einfallenden statt der gestreuten
   Welle: E⊥,FPWPM = r(θ)·S·ê mit r = cos θ (3 − cos θ)/2, S der skalaren Born-Amplitude und ê der
   Einfallspolarisation. In 3D an einer dünnen Zufallsschicht bestätigt (K17, auf 0,1–2 %).

Für die Messgröße von S9b-2 (|E⊥|² gegen skalares Born) folgt daraus das Verhältnis r²; gegen das
physikalisch richtige Paar (volles Feld gegen vektorielles Born) dagegen r²(1 + tan²θ cos²φ)/(1 −
sin²θ cos²φ), azimutal gemittelt r²/cos²θ – die FPWPM streut dann bei großen Winkeln **zu viel**
(1,25 bei 40°), weil der p-Anteil um r²/cos⁴θ zu groß ist und der s-Anteil um r² zu klein.

## Zahlen

`python experiments/vector_born_wm.py` (BWPM), gemessene Werte aus `docs/bericht-s9b2.md`:

| l_c (µm) | D | 1 − g Born skalar (Referenz S9b-2) | 1 − g Born vektoriell | 1 − g K17-Vorhersage für die Messgröße | 1 − g gemessen (S9b-2) | μs vektoriell/skalar |
| --- | --- | --- | --- | --- | --- | --- |
| 0,5 | 4 | 0,0538 | 0,0447 | 0,0386 | 0,043 | 0,951 |
| 0,2 | 4 | 0,1693 | 0,1384 | 0,1089 | 0,112 | 0,854 |
| 1,0 | 4 | 0,0191 | 0,0163 | 0,0147 | 0,0172 | 0,982 |
| 2,0 | 4 | 0,0062 | 0,0055 | 0,0051 | 0,0066 | 0,994 |
| 0,5 | 3,5 | 0,0783 | 0,0631 | 0,0515 | 0,0562 | 0,929 |
| 0,5 | 5 | 0,0275 | 0,0245 | 0,0229 | 0,0257 | 0,974 |

- Die skalaren Werte stimmen mit der Tabelle in S9b-2 überein (gleiche Referenz).
- Die K17-Vorhersage (einfache Streuung) liegt 3–23 % unter den Messwerten. Die Richtung passt
  zur Mehrfachstreuung bei der Auswertungstiefe 0,25/μs, die 1 − g vergrößert (S9b-2, Beobachtung
  1); quantitativ ist das nicht geprüft.
- Das μs-Defizit von 1–6 % gegen skalares Born hat die Größe des Faktors μs(vektoriell)/μs(skalar);
  wie die FPWPM die Extinktion des kohärenten Feldes wiedergibt, ist in K17 nicht untersucht.

Winkelverteilung (`phase_bwpm_vs_born.csv`, beide auf die Halbkugel normiert, Bereiche von 9°):
gemessen/skalar 1,06 / 1,06 / 1,04 / 0,91 / 0,75 / 0,55 bei 13,5° / 22,5° / 31,5° / 40,5° /
49,5° / 58,5°, K17-Vorhersage 1,04 / 0,98 / 0,90 / 0,78 / 0,63 / 0,45. Die Richtung stimmt, die
Kurve nicht genau (eine Realisierung, Mehrfachstreuung).

## Nachtrag: μs aus der Extinktion (Rückmeldung in #29, S9b-4)

Die Annahme „Kohärenzverlust = gestreute Querleistung, ∫ σ r² dΩ“ trifft nicht zu. Die FPWPM
erhält die **Norm** Σ|E⊥|² (in der Kalibrierung max |Σ|E⊥|² − 1| = 0,003), nicht den Fluss. Eine
Zelle im k⊥-Gitter mit |E⊥|² trägt zur Norm |E⊥|² bei, zum Fluss |E⊥|² cos θ (bei skalarer
Behandlung). Der Kohärenzverlust ist deshalb die gestreute Norm, im Ensemble
∫_vorwärts σ r²/cos θ dΩ. In BWPM je Realisierung bestätigt (K17, Punkt 8): Verlust = Σ r²|S|² auf
0,2–1,7 %, 8–9 % über dem vektoriellen Born-Fluss.

| l_c (µm) | D | μs vektoriell/skalar | ∫ σ r² dΩ (Annahme S9b-4) | μs FPWPM/skalar (∫ σ r²/cos θ dΩ) | gemessen S9b-2/S9b-4 |
| --- | --- | --- | --- | --- | --- |
| 0,5 | 4 | 0,951 | 0,932 | 0,978 | 0,96 ± 0,03 / 0,953 |
| 0,2 | 4 | 0,854 | 0,768 | 0,892 | 0,902 |
| 1,0 | 4 | 0,982 | 0,977 | 0,994 | 0,972 |
| 0,5 | 3,5 | 0,929 | 0,894 | 0,957 | 0,961 |
| 0,5 | 5 | 0,974 | 0,969 | 0,995 | 0,993 |

Die Einzelläufe passen bei l_c = 0,2 µm, D = 3,5 und D = 5 besser zur FPWPM-Vorhersage als zum
vektoriellen Born, bei l_c = 1,0 µm und im Mittel von l_c = 0,5 µm, D = 4 umgekehrt; die Streuung
der Realisierungen (±3 %) lässt im Ensemble keine sichere Aussage zu. Folgerung: Die Nähe zu
vektoriellem Born ist kein Beleg für vektorielle Streuung (wie in #29 vermutet). Das μs aus der
Extinktion der BWPM ist aber auch nicht das physikalische μs, sondern ≈ ∫ σ r²/cos θ dΩ, also nahe
am skalaren Born. Für eine Kalibrierung des Monte Carlo ist deshalb das vektorielle Born die
bessere Referenz als die gemessene Extinktion, solange #37 offen ist.

## Vorschläge

1. **Referenz vektoriell:** σ_vec(θ, φ) = σ(θ)·(1 − sin²θ cos²φ) für x-Polarisation, für die
   Phasenfunktion im Monte Carlo (unpolarisiert) σ(θ)·(1 + cos²θ)/2; μs und g daraus neu.
2. **Leistung mit E_z:** in `angular_distribution` und `calibrate_with_bwpm`
   |F{E_x}|² + |F{E_y}|² + |k_x F{E_x} + k_y F{E_y}|²/k_z² statt der Summe der Querkomponenten.
3. **Grenze der FPWPM beachten:** quantitativ (< 1 % Amplitude) nur bis θ ≈ 11° (s) bzw. 6,6° (p).
   Wird die Phasenfunktion aus der BWPM gemessen (Folgerung D9), trägt sie den Fehler aus K17; mit
   den Vorschlägen 1 und 2 ist er für 1 − g bei l_c ≥ 0,5 µm klein gegen den Fehler der skalaren
   Referenz, bei l_c = 0,2 µm nicht. Ob BWPM das Verfahren korrigiert, ist offen (#29).
