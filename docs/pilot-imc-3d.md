# Fusionspilot in 3D: Faserfacette mit Index-Matching-Kappe (B5a)

Teil B5a von B5: der Pilot aus [pilot-imc-2d.md](pilot-imc-2d.md) in 3D mit dem CUDA-Kern aus B4,
ohne Raytatouille. Die Kopplung mit Raytatouille (B5b) ist nicht Teil dieses Berichts.
Code: `bwpm.pilot_imc3d` (Auswertung), `experiments/b3_scaling.py` (Messlauf); Tests der
Auswertung: `tests/test_pilot_imc3d.py`; Messungen: Abschnitte „Messlauf B5a“ und „Diagnose B5a“ in
[testfaelle.md](testfaelle.md) (Testfall K6 in 3D), Aufträge:
[agent-server-b5a.md](agent-server-b5a.md), [agent-server-b5a-diagnose.md](agent-server-b5a-diagnose.md).

**Stand:** gemessen, Abnahmekriterien noch nicht festgelegt (Maintainer). Die Auswertung ist durch
Tests gegen analytische Werte abgesichert, die Messwerte selbst entstehen nur auf dem Server und
sind nicht Teil der CI.

## Aufbau

Geometrie und Abtastung wie in 2D (wende2023 §2.2.2): λ = 0,635 µm, trockenes Objektiv NA 0,4,
Faserkern n = 1,5 (Radius 1 µm), Mantel n = √(1,5² − 0,39²) = 1,4485, Kappe aus IP-Q, n = 1,511,
50 µm dick, Facette bei z = 0. Eingangsebene 20 µm vor der Kappe. Abtastung
dx = dy = dz = λ/(f·n_max); f = 5 entspricht der Quelle. Fenster 96 oder 128 µm mit Randmaske
(20 λ), zwei Iterationen, CUDA-Kern.

Das einfallende Feld ist das Richards-Wolf-Spektrum eines aplanatischen Objektivs (x-polarisierte
Pupille, Amplitude ∝ cos^(−1/2) θ, mit den Vektorfaktoren der Querkomponenten), analytisch bis zur
Eingangsebene geführt. Der Luftfokus liegt 32,42 µm hinter der Kappenoberfläche, übernommen aus
dem 2D-Pilot; der beste Fokus in 3D wird gemessen und nicht vorgegeben.

**Unterschied zu 2D:** Die Kette läuft in 3D vollständig im Skript. Das Objektivfeld ist analytisch
(kein Field Tracing), und der Fokus des reflektierten Feldes wird mit `pilot_imc3d.axial_focus`
aus dem Winkelspektrum bestimmt. Port-Kopplung mit Field Tracing oder Raytatouille wurde in 3D
nicht gerechnet.

## Auswertung

Alle Größen aus `bwpm.pilot_imc3d`, Einheiten wie in `docs/konventionen.md`:

- **Reflexionsgrad:** Poynting-Fluss des reflektierten Feldes am Eingang geteilt durch den des
  einfallenden Feldes. Referenz: flussgewichtete Fresnel-Reflexion Luft → Kappe, jede Mode mit
  R_s bzw. R_p ihres eigenen Winkels (`reference_reflectance`).
- **Reflexfokus:** Maximum der Intensität auf der Achse beim Rückwärtsführen des reflektierten Feldes
  im Winkelspektrum (`axial_focus`). Referenz: Spiegelbild des gebauten Luftfokus,
  −(50 + 32,42) µm = −82,42 µm.
- **Rücklaufender Fluss an der Facette:** rücklaufender durch vorlaufenden Fluss 1 µm vor der
  Facette, in der Kappe. Referenz (grobe Näherung): Fresnel bei senkrechtem Einfall Kappe → Kern
  (1,3·10⁻⁵) bzw. → Mantel (4,5·10⁻⁴), gewichtet mit dem gemessenen Kernanteil des Vorwärtsflusses
  (r ≤ 1 µm).
- **Strehl-Verhältnis** und bester Fokus auf der Facette wie in B3/K15.

Die Hilfsfunktionen sind gegen analytische Werte getestet: cos-θ-Gesetz des Flusses für s und p,
Fresnel bei senkrechtem Einfall und am Brewster-Winkel, bekannter Defokus.
Der Fokusfinder reproduziert in jedem Lauf den gebauten Fokus des einfallenden Feldes auf
8·10⁻¹⁴ µm (`incident_focus_check_um`).

## Ergebnisse

Mit Kappe, f = 5 (Quellauflösung, 1142² × 1011 bei 96 µm, 1524² × 1011 bei 128 µm), 2 Iterationen,
RTX 6000 Ada:

| Größe | Fenster 96 µm | Fenster 128 µm | Referenz | 2D (K6) | wende2023 |
| --- | --- | --- | --- | --- | --- |
| Reflexionsgrad | 4,143 % | 4,153 % | 4,158 % (Fresnel, spektrumgewichtet, s und p) | 4,45 % (TE) | ≈ 4 % |
| Fokus des Reflexes (Facette bei 0) | −82,42 µm | −82,42 µm | −82,42 µm | −82,40 µm | Δz ≈ 83 µm |
| rücklaufender Fluss an der Facette | 8,35·10⁻⁵ | 8,42·10⁻⁵ | 8,9·10⁻⁵ / 9,0·10⁻⁵ (grob) | 5,6·10⁻⁵ | < 0,1 % |
| Kernanteil des Vorwärtsflusses | 0,826 | 0,824 | – | nicht gemessen | – |
| Strehl-Verhältnis | 0,9969 | 0,9972 | – (in 3D keine analytische Referenz) | 0,9972 | 0,97 |
| Laufzeit, Speicher | 245 s, 11,3 GB | 410 s, 20,0 GB | | | |

Beste Fokuslage auf der Facette in den Läufen mit 96 µm: +0,12 bis +0,26 µm.

### Konvergenz

- **Iterationen** (f = 2,5): Von 1 über 2 auf 4 Iterationen ändern sich Reflexionsgrad
  (< 10⁻⁵ relativ), Reflexfokus (nicht), rücklaufender Fluss (0,03 % zwischen 1 und 2, danach
  nicht mehr) und Strehl (6·10⁻⁶) praktisch nicht.
- **Abtastung** (f = 2,5 / 3,5 / 5): Der rücklaufende Fluss ändert sich um 0,2 %, der Kernanteil
  um höchstens 0,004, der Strehl um 4·10⁻⁵, der Reflexfokus bleibt bei −82,40 / −82,32 /
  −82,42 µm (±0,1 µm, in der Größe von dz = 0,17 / 0,12 / 0,084 µm). Der Reflexionsgrad sinkt mit f
  (96 µm: 4,1578 → 4,1525 → 4,1426 %; 128 µm: 4,1673 → 4,1664 → 4,1525 %), in beiden Fenstern um
  etwa 0,015 Prozentpunkte von f = 2,5 auf 5. Die Schritte werden nicht kleiner; ein Plateau ist
  bis f = 5 nicht erreicht.
- **Fenster** (96 gegen 128 µm, gleiches f): Der Reflexionsgrad ist bei 128 µm um 0,23–0,34 %
  (relativ) höher, der Strehl um etwa 3·10⁻⁴. Das 96-µm-Fenster hat damit einen messbaren
  Abschneidefehler in dieser Größe.

Der Abstand des Reflexionsgrads zur Referenz bleibt unter 0,4 %. Sein Vorzeichen hängt vom Fenster
ab (96 µm: durchweg negativ, 128 µm: +0,22 / +0,19 / −0,14 %); der Trend mit f ist davon
unabhängig.

## Bewertung

- **Reflexfokus:** stimmt mit dem Spiegelbild des Luftfokus auf 0,1 µm überein, gleich dem 2D-Ergebnis
  (−82,40 µm). Die Kappe verschiebt den Rückreflex wie in der Quelle (Δz ≈ 83 µm) aus der
  Zwischenbildebene.
- **Reflexionsgrad:** im Rahmen von 0,4 % an der unabhängigen Referenz. Der Unterschied zu den
  4,45 % in 2D liegt in dem, was verglichen wird: 2D rechnet TE, die Fresnel-Reflexion von TE steigt
  mit dem Winkel; in 3D mischen s und p, und p fällt zum Brewster-Winkel ab.
- **Rücklaufender Fluss:** 3D-Wert 8,4·10⁻⁵, unter der Quellgrenze von 0,1 %. Er ist in Abtastung,
  Iterationen und Fenster stabil (0,2 %, 0,03 %, 0,9 %) und liegt 6–7 % unter der groben Referenz.
  Der Faktor 1,5 gegenüber 2D (5,6·10⁻⁵) passt zum Kernanteil: In 3D trifft 82–83 % des Flusses den
  Kern. Mit diesem Anteil liefert die Fresnel-Rechnung 8,8–9,0·10⁻⁵; für 5,6·10⁻⁵ wären etwa 90 %
  nötig. Dass ein runder Fokus mehr Mantelfläche trifft als ein Linienfokus, ist plausibel, aber
  nicht gemessen, weil der Kernanteil in 2D nicht ausgewertet wurde.
- **Strehl:** 0,9969–0,9972 in 3D wie in 2D (0,9972). Die Frage der Definition gegenüber der
  Quelle (0,97) bleibt offen, siehe 2D-Bericht.
- **Gültigkeit (E5):** Alle reflektierenden Flächen stehen senkrecht zu z. Der Fall liegt im
  quantitativen Bereich der BWPM.

## Offene Punkte

1. **Abwärtstrend des Reflexionsgrads mit f.** Etwa −0,37 % relativ von f = 2,5 auf 5, in beiden
   Fenstern gleich, ohne erkennbares Plateau. Ursache nicht untersucht. Ein Lauf bei f = 7 und 96 µm
   (etwa 31 GB) würde zeigen, ob der Trend weiterläuft; er liegt über der Quellauflösung.
2. **Fenster.** Der Fehler durch das 96-µm-Fenster (0,2–0,3 % im Reflexionsgrad, 3·10⁻⁴ im Strehl)
   ist gemessen, aber ohne Konvergenz im Fenster (nur 96 und 128 µm).
3. **Kernanteil in 2D** nicht gemessen; die Erklärung des Faktors 1,5 ist damit nur in 3D
   belegt, nicht im Vergleich.
4. **Keine unabhängigen Referenzen** für Strehl und für den rücklaufenden Fluss in 3D (nur grobe
   Fresnel-Näherung bei senkrechtem Einfall, Abweichung 6–7 %, nicht untersucht).
5. **Kopplung (B5b):** Port-Übergabe an Raytatouille und Leistungsnormierung gegenüber Raytatouille
   hängen von Raytatouille v1.0 und dem Port ab (Fusionsplan).

## Reproduktion

```bash
cd experiments
python b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2 --core --backend cuda
python b3_scaling.py --fractions 2.5 3.5 5 --window-um 128 --iterations 2 --core --backend cuda
```

Die JSON-Zeilen enthalten `reflectance`, `reflectance_ref`, `reflex_focus_um`,
`reflex_focus_ref_um`, `backward_over_forward_at_facet`, `backward_over_forward_ref`,
`core_flux_fraction`, `incident_focus_check_um`, `strehl`, `best_focus_um`. Ein Lauf bei f = 1,2 mit
dem CPU-Kern (`--backend cpu`, 20 s) prüft die Kette ohne GPU.
