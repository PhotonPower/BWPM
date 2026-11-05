# Fusionspilot in 2D: Faserfacette mit Index-Matching-Kappe

Vorgezogener Teil von B5 (Fusionsplan, Pilot K4; hier Testfall K6), in 2D und in Python.
Code: `bwpm.pilot_imc`, Tests: `tests/test_pilot_imc.py` (brauchen das Extra `fusion`, also das
Paket `fieldtracing`, festgelegt auf Commit `6da9983`).

## Aufbau

Nach wende2023 §2.2.2: λ = 0,635 µm, TE, trockenes Objektiv NA 0,4; Faserkern n = 1,5, Durchmesser
2 µm, NA_Kern = 0,39 (Mantel n = √(1,5² − 0,39²) = 1,4485); Kappe aus IP-Q, n = 1,511, 50 µm dick.
Abtastung dy = dz = λ/(5·n_max) ≈ 0,084 µm wie in der Quelle, Fenster 1152 Punkte (97 µm) mit
Randmaske (20 λ, Ordnung 8), zwei Iterationen.

- **Ohne Kappe:** Facette im Fokus des Objektivs.
- **Mit Kappe:** Der Luftfokus des Objektivs liegt x = 32,42 µm hinter der Kappenoberfläche; so
  liegt der beste Fokus in der Kappe genau auf der Facette (bestimmt aus dem Strehl-Verhältnis,
  paraxial wären es 50 µm/1,511 = 33,09 µm).

## Kette

1. **Field Tracing** (`fieldtracing.pft.pft_1d`): ideales aplanatisches Zylinderobjektiv,
   Pupillenfeld U·e^(iψ) mit ψ = −k(√(F² + y²) − F) und U = cos θ (ergibt die aplanatische
   Amplitude ∝ √cos θ je Winkel). Die PFT liefert das Spektrum an den Gitterfrequenzen der BWPM;
   die Spektralphase wird analytisch bis zur Eingangsebene geführt, 20 µm vor der ersten Fläche.
2. **BWPM** (`bwpm.bidirectional`): vor- und rücklaufende Felder durch Luft, Kappe und Faser.
3. **Field Tracing** (`fieldtracing.propagation.angular_spectrum`): Das reflektierte Feld wird mit
   dem exakten Winkelspektrum weitergeführt, um seinen Fokus zu finden.

Die Übergabe zwischen 1 und 2 und zwischen 2 und 3 ist jeweils ein Gitterfeld in einer Ebene
z = const mit Querkomponenten in den gemeinsamen Konventionen – genau der Port, den der
Fusionsplan für die BWPM vorsieht.

## Ergebnisse

| Größe | ohne Kappe | mit Kappe | Referenz | wende2023 (3D) |
| --- | --- | --- | --- | --- |
| reflektierte Leistung am Eingang | 4,23 % | 4,45 % (an der Kappenoberfläche) | Fresnel, über das Spektrum gemittelt: Kern 4,3 %, Mantel 3,7 %; Luft → IP-Q 4,45 % | ≈ 4 % |
| Fokus des reflektierten Lichts (Facette bei 0) | 0,02 µm | −82,40 µm | Spiegelbild des Luftfokus: −(50 + 32,42) = −82,42 µm | Δz ≈ 83 µm |
| rücklaufende Leistung zwischen Kappe und Facette | (≈ 4 %) | 5,6·10⁻⁵ | Fresnel IP-Q → Kern 1,3·10⁻⁵, → Mantel 4,5·10⁻⁴ | < 0,1 % |
| Strehl-Verhältnis auf der Facette | – | 0,9972 | analytisch (Plattenphase, t_s): 0,9976 | 0,97 |
| bester Fokus | Facette | Facette (± 0,02 µm) | – | – |

Laufzeit je Konfiguration etwa 2–4 s.

**Bewertung:**

- Die Kette funktioniert: Ein Feld aus Field Tracing geht als Gitterfeld in die BWPM, das
  reflektierte Feld als Gitterfeld zurück in Field Tracing, und alle Größen stimmen mit
  unabhängigen Referenzen überein (Fokus des Reflexes auf 0,02 µm, Strehl auf 4·10⁻⁴).
- Die Kappe verschiebt den Rückreflex um 82,4 µm aus der Zwischenbildebene und senkt die
  Rückreflexion an der Facette von 4,2 % auf 5,6·10⁻⁵ – das Ergebnis der Quelle wird in 2D
  bestätigt.
- Das Strehl-Verhältnis ist höher als in der Quelle (0,997 statt 0,97). **Korrektur:** Auch die
  Quelle rechnet in 2D (wende2023 §2.1.2), der Unterschied ist also kein 2D/3D-Effekt; ein erster
  kleiner 3D-Lauf (B3, grobes Gitter) ergibt 0,993. Wahrscheinlich unterscheiden sich die
  Definitionen: Die Quelle vergleicht normierte Spitzenintensitäten der Fokusse mit und ohne
  Kappe (Abb. 5(e)); darin stecken zusätzlich etwa 4 % Fresnel-Verlust an der Kappenoberfläche
  und das Objektivmodell. Hier ist es das Verhältnis zum aberrationsfreien Spektrum gleicher
  Amplitude. Offen, nicht als Abweichung gemeldet.
- Der Fall liegt im quantitativen Bereich der BWPM (E5): alle reflektierenden Flächen stehen
  senkrecht zu z. Die geführten Moden im Kern laufen nur vorwärts; ihre Rückreflexe innerhalb der
  Faser wären der nicht konvergierende Fall aus K4 und kommen hier nicht vor.

## Was der Pilot für B4 und die Fusion festlegt

1. **Port:** Querfeld (E_x, E_y) auf einem zentrierten Gitter y_n = (n − N//2)·dy in einer Ebene,
   Richtung (+z oder −z), Medium am Port, Wellenlänge. Das Spektrum ist auf y = 0 zu beziehen
   (Phasenfaktor e^(i k_y y_0) gegenüber der FFT); ohne ihn lag der Fokus im ersten Versuch am
   Fensterrand.
2. **Übergabe Field Tracing → BWPM:** Spektrum an den Gitterfrequenzen der BWPM aus der PFT,
   Phase analytisch bis zur Port-Ebene (keine Interpolation der Phase).
3. **Übergabe BWPM → Field Tracing:** Gitterfeld; Rückwärtsausbreitung über das Winkelspektrum
   nur für propagierende Anteile (evaneszente wachsen bei negativen Schritten).
4. **Leistung:** Poynting-Fluss je Port (`power_flux`); die Leistungsnormierung gegenüber
   Raytatouille folgt in B4.
