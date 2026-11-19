# Entwurf B8: Rückstreuung an schwachen Fluktuationen δn

> Entwurf zur Freigabe (2026-11-19). Grundlage: Stufe B7 (Born-Inkrement vorwärts,
> `docs/spezifikation.md` §10, K18), chen2020 mit Zusatzmaterial (Volltext geprüft), BWPM-Ablauf
> nach wende2023 (§5). Anlass: Scattuccino braucht das rückgestreute Feld (OCT, diffuse Reflexion);
> heute streut nur der Vorwärtslauf an δn (E10).

## 1. Was fehlt

In B7 erzeugt δn nur ein Vorwärtsinkrement. In der BWPM fehlen damit drei Beiträge:

| Beitrag | Herkunft | chen2020 Zusatzmaterial |
| --- | --- | --- |
| **R1** Rückstreuung des Vorwärtsfeldes an δn | Quelle für den Rückwärtslauf | U_bs, Gl. (S3), (S4) |
| **R2** Vorwärtsstreuung des Rückwärtsfeldes an δn (Inkrement in −z) | Rückwärtslauf | U_bi,fs, Gl. (S2) |
| **R3** Rückstreuung des Rückwärtsfeldes an δn (wieder nach +z) | Quelle für die nächste Iteration | in chen2020 nicht enthalten (dort nur ein Rückweg) |

R1 ist das OCT-Signal in einfacher Rückstreuung. R2 ist die Mehrfachstreuung auf dem Rückweg,
die bei dicken Proben die Bildqualität bestimmt. R3 ergibt sich in der BWPM ohne Zusatzaufwand,
wenn R1 umgesetzt ist (gleiche Formel, Richtung gespiegelt), und ist nötig, damit die Iteration
der BWPM wie bei Grenzflächen gegen die exakte Lösung konvergiert (K4).

## 2. Formel

Wie §10, mit der gestreuten Welle in −z-Richtung, k_s = (k_x, k_y, −k_z):

ΔE⊥,rück(z_j) = F⁻¹{ i/(2k_z) · [P(k_s)·∫_Scheibe e^(i k_z (ζ − z_j)) F{V·E(ζ)} dζ]⊥ }.

- Vorwärts ist die Phasendifferenz im Integral (k − k_z)dz klein; die Mittelpunktregel (sinc ≈ 1)
  genügt (B7).
- Rückwärts ist sie (k_z + k_z′)dz ≈ 2k·dz, also groß. chen2020 (Zusatzmaterial, nach Gl. (S4))
  sagt ausdrücklich, dass sinc ≈ 1 dann nicht gilt, und entwickelt sinc als Taylor-Reihe bis zur
  Ordnung t (Gl. (S5)–(S7), t = 8 in ihren Rechnungen, dort mit dz = 2,44 µm).

Fehler des Scheibenintegrals für λ = 0,8 µm, ñ₀ = 1,38 (Faktor gegen das exakte Integral):

| dz | Richtung, Winkel | Mittelpunkt | Gauß, 2 Punkte | Taylor t = 2 |
| --- | --- | --- | --- | --- |
| 0,04 µm | rück, 0° | 3,2·10⁻² | 1,3·10⁻⁴ | 3,0·10⁻⁴ |
| 0,04 µm | rück, 60° | 1,8·10⁻² | 4,2·10⁻⁵ | 9,5·10⁻⁵ |
| 0,02 µm | rück, 0° | 7,9·10⁻³ | 8,2·10⁻⁶ | 1,9·10⁻⁵ |
| 0,04 µm | vor, 60° | 2,0·10⁻³ | 5,1·10⁻⁷ | 1,2·10⁻⁶ |

Kosten je Scheibe (FFTs über das Querfeld, je Komponente): Gauß mit M Punkten braucht M-mal
„Feld an ζ_m, mal V, zurück transformieren“; Taylor bis t braucht t + 1 solche Paare. Die
Quellen F{V·E(ζ_m)} können Vorwärts- und Rückwärtsinkrement **gemeinsam** nutzen.

## 3. Entscheidungen

| Nr. | Frage | Vorschlag | Begründung |
| --- | --- | --- | --- |
| F1 | Quadratur des Scheibenintegrals | **Gauß-Legendre mit 2 Punkten** für vorwärts und rückwärts, gemeinsame Quellen | Fehler 10⁻⁴ bei dz = 0,04 µm, weniger FFTs als Taylor t = 2; weicht von chen2020 (Taylor) ab → Issue mit Label `paper`. Das Vorwärtsinkrement wird dadurch genauer als in B7; K18 bleibt grün, die Sollwertdateien von B7b werden neu erzeugt |
| F2 | Umfang | **R1, R2 und R3** | R1 allein (einfache Rückstreuung) unterschätzt den Rückweg durch dickes Gewebe; R3 kostet nichts zusätzlich |
| F3 | Speicher der Rückstreuquellen | **volle Ebenen je Scheibe mit δn**, complex128; complex64 als spätere Option | wie die Reflexquellen der BWPM, aber nicht dünn besetzt. Beispiel: 256² × 250 Scheiben ≈ 0,5 GB; 512² × 1000 ≈ 8,4 GB (passt auf die RTX 6000 Ada mit 48 GB, nicht auf die RTX 2080 mit 8 GB) |
| F4 | Iterationen | wie heute (`iterations`), δn-Quellen wie Grenzflächenquellen behandelt | ein Ablauf für beides, Konvergenz wie K4 prüfbar |
| F5 | Wirtsmedium | wie B7: je Scheibe ein Wirt | unverändert |

Alternative zu F3, falls der Speicher nicht reicht: nur einfache Rückstreuung (R1) und die
Beiträge aller Scheiben im Spektrum einer Ebene aufsummieren (Speicher eine Ebene); gilt nur,
wenn zwischen Streuer und Eintrittsebene nur der Wirt liegt, und lässt R2 und R3 weg. Nicht
empfohlen für OCT in dickem Gewebe.

## 4. Prüfweg (K19)

Sollwerte unabhängig vom Code:

1. **Homogene Schicht ñ₀ + δn** (2D, senkrecht und schräg, s und p): Reflexion gegen die
   Airy-Formel. In erster Ordnung in δn (einfache Rückstreuung, eine Iteration) und mit
   Iterationen gegen die volle Airy-Lösung, Konvergenz wie K4.
2. **Bragg-Gitter in z** (δn periodisch in z, Periode λ/(2ñ₀)): Reflexionsmaximum gegen die
   Transfer-Matrix des diskretisierten Stapels. Mehrfachstreuung und Phasenanpassung.
3. **Gitter quer (Lamellen) in Reflexion**: Beugungsordnungen der Reflexion gegen die RCWA
   (`tests/rcwa_reference.py` um Reflexion erweitern) und gegen eine Born-Referenz der Rückstreuung
   (`bwpm.born` erweitern, wie K17).
4. **Energie**: R + T = 1 für verlustfreies δn mit Iterationen, Leistung mit E_z.
5. **Konvergenz** in dz (Gauß zweiter Ordnung) und in den Iterationen.
6. **δn = 0 bitgleich**, K1–K18 und Sollwertdateien unverändert (bis auf die Neuerzeugung nach F1).

## 5. Ablauf

- **B8a** Python-Referenz in 2D und 3D, K19.
- **B8b** CPU- und CUDA-Kern, neue Sollwertdateien.

Wie B7 unabhängig von B5b.
