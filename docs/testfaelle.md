# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ offen.

| Nr. | Prüfung | Referenz | Abnahme | Stufe | Status | Test |
| --- | --- | --- | --- | --- | --- | --- |
| K1 | Transmissions- und Reflexionsmatrix der Quelle (Gl. (2), (3), (8), (9)) an einer ebenen Grenzfläche | vektorieller Fresnel-Operator in Convention A (unabhängig aufgebaut) | Abweichung < 1e−13 relativ, auch absorbierend, Totalreflexion, senkrechter Einfall | B0 | ✅ | `test_interface.py` |
| K2 | homogenes Medium | ebene Welle (Dispersionsrelation), einmalige Winkelspektrum-Propagation | < 1e−12 | B1 | ✅ | `test_fpwpm.py` |
| K3 | ebene Grenzfläche, schräger Einfall, FPWPM | Fresnel mit Winkeln (TE: t_s, TM: t_p cos θ₂/cos θ₁); Leistung eines Gaußstrahls aus dem Winkelspektrum | < 1e−12, auch absorbierend und Silber | B1 | ✅ | `test_fpwpm.py` |
| K4 | Platte und Mehrschichtstapel schräg beleuchtet, BWPM | Transfer-Matrix des treppenförmig diskretisierten Stapels (exakt); Airy-Formel | Abweichung nach n Iterationen → Maschinengenauigkeit, TE und TM, auch Silber; Konvergenz über Iterationen wie wende2023 Abb. 2 | B2 | ○ | |
| K7 | Randzone der absorbierenden Maske | gleiche Rechnung im vierfach breiten Fenster ohne Maske | Fehler im Inneren < 1e−5 bei Randzone ≥ 20 λ, Ordnung 8 | B1 | ✅ | `test_boundary_and_wedge.py` |
| K8 | gekippter Glaskeil (Flächenwinkel 0–60°) | Snellius und Fresnel an der gekippten Fläche | Richtung < 0,02° bis 45°; Leistung weicht bis 7 % ab (Grenze des Verfahrens); R folgt in B2 | B1 | ✅ | `test_boundary_and_wedge.py` |
| K5 | Parabolspiegel (Silber) | Brennweite f′ | Fokuslage, wende2023 Abschnitt 2.1.4 als Plausibilität | B2/B3 | ○ | |
| K6 | Faserfacette mit Index-Matching-Kappe | Fresnel an der Facette; wende2023 Abschnitt 2.2.2 als Plausibilität | Rückreflex, Fokusversatz | B5 | ○ | |

## K7: Randzone der Maske (B1)

TE-Gaußstrahl w₀ = 2 µm (4 λ), um 30° gekippt, 40 µm durch Luft, λ = 0,5 µm, dy = dz = λ/10,
Fenster 1024 Punkte (51,2 µm); Referenz im Fenster mit 4096 Punkten ohne Maske. Größter Betrag der
Abweichung im Inneren (außerhalb der Randzone), über Ebenen alle 2,5 µm, Feldamplitude anfangs 1:

| Randzone | ohne Maske | 10 λ | 15 λ | 20 λ | 30 λ |
| --- | --- | --- | --- | --- | --- |
| Ordnung 8 | 0,48 | 3,9·10⁻⁴ | 1,1·10⁻⁵ | 1,7·10⁻⁶ | 5,9·10⁻⁸ |

Ordnung 4 ist bei 20 λ etwa fünfmal, Ordnung 16 etwa achtmal schlechter als Ordnung 8.
**Regel:** Randzone ≥ 20 λ, Ordnung 8. **Grenze:** Bei 60° Kippung und w₀ = 4 λ bleibt ein Fehler
von etwa 10⁻² unabhängig von der Randzone; dort hat das Spektrum merkliche Anteile nahe streifendem
Einfall, die die Maske nicht sauber abfängt. Für Ausbreitungswinkel bis 30° (NA 0,5) gilt die Regel.

## K8: gekippter Glaskeil (B1)

Fläche z = 10 µm + y·tan α, n = 1,5, Gaußstrahl w₀ = 5 µm senkrecht einfallend, λ = 0,5 µm,
Fenster 102 µm mit Maske (20 λ, Ordnung 8). Richtung des durchgelassenen Strahls (Schwerpunkt des
Spektrums) gegen Snellius, Leistung gegen die Fresnel-Transmission an der gekippten Fläche:

| α | Regel, d | Richtungsfehler | Leistung (Fresnel) TE |
| --- | --- | --- | --- |
| 0° | Mittelpunkt, λ/10 | 0 | 0,96000 (0,96000; Abweichung 1,4·10⁻⁵ aus der Winkelbreite) |
| 10° | Mittelpunkt, λ/10 / λ/20 | −0,008° / −0,002° | 0,969 / 0,997 (0,958) |
| 20° | Mittelpunkt, λ/10 / λ/20 | +0,012° / +0,001° | 1,009 / 1,012 (0,953) |
| 30° | Mittelpunkt, λ/10 / λ/20 | 0,000° / 0,000° | 1,009 / 1,000 (0,942) |
| 45° | Mittelpunkt, λ/10 / λ/20 | −0,002° / −0,002° | 0,955 / 0,945 (0,908) |
| 60° | Mittelpunkt, λ/10 / λ/20 | +1,98° / +1,85° | 0,846 / 0,838 (0,823) |
| 30° | Flächenanteil, λ/10 / λ/20 | −0,001° / −0,001° | 1,010 / 0,999 (0,942) |

TM bei λ/10: Richtungsfehler ≤ 0,004° bis 45°, Leistung 0,972/0,959/0,969 gegen 0,962/0,975/0,992
(10°/30°/45°); bei 60° wieder etwa 2°.

**Befund:**

1. Die **Richtung** stimmt bis 45° Flächenwinkel auf ≤ 0,01°; bei 60° weicht sie um etwa 2° ab.
2. Die **Leistung** weicht um bis zu 7 % von Fresnel an der gekippten Fläche ab, für TE sogar über
   1 (keine Energieerhaltung), und **konvergiert nicht in dz**. Ursache ist die Näherung des
   Verfahrens: Jede Stufe der Treppe ist eine z-senkrechte Grenzfläche, Fresnel wird zum Winkel gegen
   z angewandt, und die Maskierung mit I_r,s ist nicht unitär (Spezifikation, Abschnitt 4).
3. Die Regel Flächenanteil (E4) ändert daran nichts Wesentliches; sie bleibt eine Option.
4. Folge für die Fusion: Die BWPM ist für **Richtungen und Phasen** gekippter Flächen bis 45°
   verlässlich, für **Leistungen** an Flächen über etwa 10° nur auf einige Prozent. Der Planer
   braucht diese Grenze als Kennzahl.

## K1: Grenzflächenmatrizen (B0)

Messwerte und Fälle: `docs/spezifikation.md`, Abschnitt 3.
