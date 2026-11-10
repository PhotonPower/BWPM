# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ offen.

| Nr. | Prüfung | Referenz | Abnahme | Stufe | Status | Test |
| --- | --- | --- | --- | --- | --- | --- |
| K1 | Transmissions- und Reflexionsmatrix der Quelle (Gl. (2), (3), (8), (9)) an einer ebenen Grenzfläche | vektorieller Fresnel-Operator in Convention A (unabhängig aufgebaut) | Abweichung < 1e−13 relativ, auch absorbierend, Totalreflexion, senkrechter Einfall | B0 | ✅ | `test_interface.py` |
| K2 | homogenes Medium | ebene Welle (Dispersionsrelation), einmalige Winkelspektrum-Propagation | < 1e−12 | B1 | ✅ | `test_fpwpm.py` |
| K3 | ebene Grenzfläche, schräger Einfall, FPWPM | Fresnel mit Winkeln (TE: t_s, TM: t_p cos θ₂/cos θ₁); Leistung eines Gaußstrahls aus dem Winkelspektrum | < 1e−12, auch absorbierend und Silber | B1 | ✅ | `test_fpwpm.py` |
| K4 | Platte und Mehrschichtstapel schräg beleuchtet, BWPM | Transfer-Matrix des treppenförmig diskretisierten Stapels (tangentiale Felder, exakt); Airy-Formel | < 1e−12 nach 60 Iterationen, TE und TM, 0–70°, auch Silber und evaneszent; ε fällt geometrisch | B2 | ✅ | `test_bidirectional.py` |
| K10 | 3D homogen, ebene Welle mit k_x, k_y ≠ 0 | Dispersionsrelation | < 1e−12 | B3 | ✅ | `test_volume3d.py` |
| K11 | 3D ebener Stapel (Luft/1,5/2,1/Silber/1,5), beliebiger Azimut und Polarisation | Transfer-Matrix für s und p getrennt (tangentiale Felder) | < 1e−12 nach 30 Iterationen; prüft die Kopplung von E_x und E_y | B3 | ✅ | `test_volume3d.py` |
| K12 | 3D in x homogen = 2D; dünn gespeicherte Reflexe = voll gespeicherte | 2D-BWPM; volle Speicherung | < 1e−12 bzw. bitgleich, Speicher < 1/5 | B3 | ✅ | `test_volume3d.py` |
| K13 | C++-Kern: Grenzflächenmatrizen (`bwpm._core`) | Sollwertdatei `interface_matrices`; Python-Referenz auf einem FFT-Gitter mit k⊥ = 0, streifend, evaneszent | Dateitoleranz bzw. < 1e−13 | B4b | ✅ | `test_core_interface.py` |
| K14 | C++-Kern: FPWPM/BWPM in 2D und 3D (`bwpm._core.Solver`, `bwpm.fast.run_2d/run_3d`) | alle Sollwertdateien unter `reference/`; Python-Referenz auf zufälligen 3D-Volumina (vier Materialien inkl. Silber, mit und ohne Maske, dx ≠ dy) | Dateitoleranz bzw. < 1e−12 des Eingangs; gespeicherte Reflexwerte gleich; Threads ändern nichts | B4c | ✅ | `test_core_solver.py` |
| K15 | CUDA-Kern (`bwpm._core.CudaSolver`, `backend="cuda"`) | alle Sollwertdateien; CPU-Kern auf zufälligen 3D-Volumina | rel 1e−10, abs 1e−11 des Maximums; CUDA gegen CPU ≈ 1e−16 | B4d | ✅ (Server, 2026-11-06) | `test_core_cuda.py` |
| K16 | Fusion-Adapter, Konformität (`bwpm.adapter`) | Transfer-Matrix (s, p getrennt) und Airy; Energieerhaltung; Reziprozität der Transmission; Determinismus; Verschiebung am zentrierten Port | < 1e−12 (CUDA 1e−10 beim Determinismus), für jedes verfügbare Backend | B4e | ✅ | `test_adapter_conformance.py` |
| K7 | Randzone der absorbierenden Maske | gleiche Rechnung im vierfach breiten Fenster ohne Maske | Fehler im Inneren < 1e−5 bei Randzone ≥ 20 λ, Ordnung 8 | B1 | ✅ | `test_boundary_and_wedge.py` |
| K8 | gekippter Glaskeil (Flächenwinkel 0–60°) | Snellius und Fresnel an der gekippten Fläche | Richtung < 0,02° bis 45°; Leistung weicht bis 7 % ab (Grenze des Verfahrens); R folgt in B2 | B1 | ✅ | `test_boundary_and_wedge.py` |
| K5 | Parabolspiegel (Silber) | Brennweite f′ (Parabel ohne Aberration) | Fokus 49,87 µm vom Scheitel bei f′ = 50 µm (wende2023: 49,85 µm), stabil bei Verfeinerung | B2 | ✅ (slow) | `test_mirror.py` |
| K9 | gekippter ebener Silberspiegel | Reflexionsgesetz, Fresnel | Richtung → 2α bei Verfeinerung; Leistung weit unter Fresnel und nicht konvergent (Grenze des Verfahrens) | B2 | ✅ (slow) | `test_mirror.py` |
| K6 | Faserfacette mit Index-Matching-Kappe, Kette Field Tracing → BWPM → Field Tracing (2D) | Fresnel, Spiegelbild des Luftfokus, analytische Plattenphase | Fokus des Reflexes −82,40 µm (Soll −82,42), Rückreflex an der Facette 4,2 % → 5,6·10⁻⁵, Strehl 0,9972 (Soll 0,9976) | B5 (2D vorgezogen) | ✅ | `test_pilot_imc.py` ([Bericht](pilot-imc-2d.md)) |
| K6 (3D) | Faserfacette mit Index-Matching-Kappe in 3D, Kette Richards-Wolf-Spektrum → BWPM → Winkelspektrum | Fresnel je Winkel und Polarisation (flussgewichtet), Spiegelbild des Luftfokus, Fresnel Kappe → Kern/Mantel am Facettenfeld | Fokus des Reflexes −82,42 µm, Reflexionsgrad = Referenz, rücklaufender Fluss an der Facette in der Größe der Referenz, Strehl wie 2D | B5a | Hilfsfunktionen ✅, Messlauf und Diagnose B5a auf dem Server gemessen (Abnahme offen) | `test_pilot_imc3d.py`, `experiments/b3_scaling.py` ([Bericht](pilot-imc-3d.md)) |

## K6 in 3D (B5a)

Erweiterung von `experiments/b3_scaling.py` (mit Kappe): Reflexionsgrad aus dem Poynting-Fluss
(`bwpm.pilot_imc3d.plane_power`, statt Summe der |E|²) gegen die flussgewichtete Fresnel-Referenz
des einfallenden Spektrums (`reference_reflectance`, je Mode R_s bzw. R_p seines Winkels),
Fokus des reflektierten Feldes (`axial_focus`, Sollwert: Spiegelbild des gebauten Luftfokus,
−(50 + 32,42) µm), rücklaufender Fluss 1 µm vor der Facette gegen eine grobe Referenz (Fresnel bei
senkrechtem Einfall Kappe → Kern bzw. Mantel, gewichtet mit dem Vorwärtsfluss über Kern und Mantel).
Die Hilfsfunktionen sind in `tests/test_pilot_imc3d.py` gegen analytische Werte geprüft (cos-θ-Gesetz,
Fresnel bei senkrechtem Einfall und am Brewster-Winkel, bekannter Defokus).

**Funktionsprüfung im Container** (CPU-Kern, ein Kern, 274² × 243, λ/(1,2 n) – zu grob für
Aussagen, 20 s): Reflexionsgrad 4,153 % gegen Referenz 4,158 %; Reflexfokus −82,48 µm gegen −82,42 µm
(Schritt dz = 0,35 µm); rücklaufend 8,3·10⁻⁵ gegen grobe Referenz 8,9·10⁻⁵ (2D, K6: 5,6·10⁻⁵);
der Fokus des einfallenden Feldes stimmt mit dem gebauten auf 10⁻¹³ µm überein; Strehl 0,9969.
Die Messung in der Auflösung der Quelle steht im folgenden Abschnitt ([Auftrag](agent-server-b5a.md)).

### Messlauf B5a

Umgebung: NVIDIA RTX 6000 Ada (49 GB, Treiber 590.48.01), nvcc 12.9.41, GCC 14.3.1, Python 3.13.5,
CUDA-Kern aus B4d (`fast.CORE_AVAILABLE` und `fast.cuda_available()` beide `True`),
`tests/test_pilot_imc3d.py` und `tests/test_core_cuda.py` 15 Tests grün. Alle Läufe mit Kappe,
Fenster 96 µm, Branch `b5a-pilot-3d-analysis`, im Ordner `experiments`:

```bash
python b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2 --core --backend cuda
python b3_scaling.py --fractions 2.5 --window-um 96 --iterations 1 --core --backend cuda
python b3_scaling.py --fractions 2.5 --window-um 96 --iterations 4 --core --backend cuda
```

**Abtastung** (2 Iterationen; f = Bruchteil der Quellauflösung, n = Gitterpunkte je Seite, nz = Schichten):

| Feld | f = 2,5 (572², 506) | f = 3,5 (800², 708) | f = 5 (1142², 1011) |
| --- | --- | --- | --- |
| `reflectance` | 4,1578 % | 4,1525 % | 4,1426 % |
| `reflectance_ref` | 4,1584 % | 4,1584 % | 4,1585 % |
| `reflectance` / `reflectance_ref` − 1 | −0,016 % | −0,14 % | −0,38 % |
| `reflex_focus_um` | −82,40 | −82,32 | −82,42 |
| `reflex_focus_ref_um` | −82,42 | −82,42 | −82,42 |
| `backward_over_forward_at_facet` | 8,356·10⁻⁵ | 8,338·10⁻⁵ | 8,348·10⁻⁵ |
| `backward_over_forward_ref` | 9,00·10⁻⁵ | 8,82·10⁻⁵ | 8,90·10⁻⁵ |
| `incident_focus_check_um` | 7,8·10⁻¹⁴ | 7,8·10⁻¹⁴ | 7,8·10⁻¹⁴ |
| `strehl` | 0,99694 | 0,99695 | 0,99691 |
| `best_focus_um` | 0,26 | 0,12 | 0,18 |
| `seconds` | 24,5 | 63,7 | 273,2 |
| `peak_rss_mb` | 1571 | 4074 | 11279 |

**Iterationen** (f = 2,5):

| Feld | 1 Iteration | 2 Iterationen | 4 Iterationen |
| --- | --- | --- | --- |
| `reflectance` | 4,157776 % | 4,157781 % | 4,157781 % |
| `reflex_focus_um` | −82,40 | −82,40 | −82,40 |
| `backward_over_forward_at_facet` | 8,3583·10⁻⁵ | 8,3555·10⁻⁵ | 8,3555·10⁻⁵ |
| `backward_over_forward_ref` | 9,0001·10⁻⁵ | 8,9977·10⁻⁵ | 8,9977·10⁻⁵ |
| `strehl` | 0,996929 | 0,996935 | 0,996935 |
| `seconds` | 26,8 | 24,5 | 43,0 |

Befund (nur Zahlen, keine Erklärung, nichts angepasst):

1. Zwischen f = 2,5, 3,5 und 5 ändern sich der Reflexionsgrad um 0,015 Prozentpunkte (relativ 0,36 %)
   und der rücklaufende Fluss um 0,2 %. Der Reflexionsgrad liegt für alle f unter der Referenz, und
   die Abweichung wächst mit f (−0,016 %, −0,14 %, −0,38 %); sie bleibt unter 1 %. Der Reflexfokus
   liegt bei f = 2,5 und 5 bis 0,02 µm an der Referenz, bei f = 3,5 bei −82,32 µm (0,10 µm daneben).
   Der Strehl ist auf 4·10⁻⁵ konstant, `best_focus_um` schwankt zwischen 0,12 und 0,26 µm.
2. Zwischen 1, 2 und 4 Iterationen ändern sich Reflexionsgrad (< 10⁻⁵ relativ) und Reflexfokus
   (nicht) praktisch nicht; der rücklaufende Fluss ändert sich zwischen 1 und 2 Iterationen um 0,03 %
   und danach nicht mehr; der Strehl um 6·10⁻⁶ und danach nicht mehr.
3. `backward_over_forward_at_facet` liegt in 3D bei 8,3·10⁻⁵ gegen 5,6·10⁻⁵ in 2D (K6), also um den
   Faktor 1,5 höher. Die grobe Referenz (8,8–9,0·10⁻⁵) liegt 6–8 % darüber.
4. Der Test der Auswertung (`incident_focus_check_um` < 10⁻⁶ µm) ist mit 7,8·10⁻¹⁴ µm erfüllt.

Die Erwartungen im Auftrag sind keine Abnahmekriterien; die Abnahme legt der Maintainer fest.
Rohdaten (`experiments/out/*.json`) sind nicht eingecheckt.

### Diagnose B5a

Umgebung wie im Messlauf B5a (RTX 6000 Ada, nvcc 12.9.41, GCC 14.3.1, Python 3.13.5, CUDA-Kern aus B4d; Branch `b5a-pilot-3d-analysis` ab `2121cf2`). Alle Läufe mit Kappe, 2 Iterationen, im Ordner `experiments`:

```bash
for w in 96 128; do
  python b3_scaling.py --fractions 2.5 3.5 --window-um $w --iterations 2 --core --backend cuda
done
python b3_scaling.py --fractions 5 --window-um 96 --iterations 2 --core --backend cuda
python b3_scaling.py --fractions 5 --window-um 128 --iterations 2 --core --backend cuda   # optional
```

**Fenster 96 µm**

| Feld | f = 2,5 | f = 3,5 | f = 5,0 |
| --- | --- | --- | --- |
| `reflectance` | 4,1578 % | 4,1525 % | 4,1426 % |
| `reflectance_ref` | 4,1584 % | 4,1584 % | 4,1585 % |
| Verhältnis − 1 | −0,016 % | −0,142 % | −0,382 % |
| `reflex_focus_um` | −82,40 | −82,32 | −82,42 |
| `backward_over_forward_at_facet` | 8,356·10⁻⁵ | 8,338·10⁻⁵ | 8,348·10⁻⁵ |
| `backward_over_forward_ref` | 8,998·10⁻⁵ | 8,816·10⁻⁵ | 8,905·10⁻⁵ |
| `core_flux_fraction` | 0,8234 | 0,8276 | 0,8255 |
| Kernanteil aus gemessenem rücklaufenden Fluss | 0,8377 | 0,8381 | 0,8379 |
| `n`, `nz` | 572², 506 | 800², 708 | 1142², 1011 |
| `seconds` | 24 | 61 | 245 |
| `peak_rss_mb` | 1574 | 4077 | 11276 |

**Fenster 128 µm**

| Feld | f = 2,5 | f = 3,5 | f = 5,0 |
| --- | --- | --- | --- |
| `reflectance` | 4,1673 % | 4,1664 % | 4,1525 % |
| `reflectance_ref` | 4,1584 % | 4,1584 % | 4,1584 % |
| Verhältnis − 1 | +0,215 % | +0,193 % | −0,141 % |
| `reflex_focus_um` | −82,40 | −82,32 | −82,42 |
| `backward_over_forward_at_facet` | 8,428·10⁻⁵ | 8,412·10⁻⁵ | 8,421·10⁻⁵ |
| `backward_over_forward_ref` | 9,067·10⁻⁵ | 8,887·10⁻⁵ | 8,972·10⁻⁵ |
| `core_flux_fraction` | 0,8218 | 0,8259 | 0,8240 |
| Kernanteil aus gemessenem rücklaufenden Fluss | 0,8360 | 0,8364 | 0,8362 |
| `n`, `nz` | 762², 506 | 1066², 708 | 1524², 1011 |
| `seconds` | 50 | 134 | 410 |
| `peak_rss_mb` | 2663 | 7105 | 19950 |

Fresnel bei senkrechtem Einfall: r_core = 1,335·10⁻⁵, r_mantel = 4,460·10⁻4.

Befund (nur Zahlen, keine Deutung, nichts angepasst):

1. **Drift und Fenster.** Bei gleichem f unterscheiden sich 96 und 128 µm im Abstand Reflexionsgrad zu
   Referenz um 0,23 / 0,34 / 0,24 Prozentpunkte (f = 2,5 / 3,5 / 5). Bei 96 µm ist der Abstand für
   alle f negativ (−0,016 / −0,142 / −0,382 %), bei 128 µm positiv bei f = 2,5 und 3,5
   (+0,215 / +0,193 %) und negativ bei f = 5 (−0,141 %). In beiden Fenstern sinkt der Reflexionsgrad
   mit f (96 µm: 4,1578 → 4,1426 %; 128 µm: 4,1673 → 4,1525 %), um jeweils 0,015 Prozentpunkte.
   Der Strehl liegt bei 128 µm um 2,5–3·10⁻⁴ höher (0,9972 gegen 0,9969).
2. **Kernanteil.** `core_flux_fraction` liegt bei 96 µm bei 0,8234 / 0,8276 / 0,8255 und bei 128 µm
   bei 0,8218 / 0,8259 / 0,8240 (f = 2,5 / 3,5 / 5). Mit f ändert er sich um höchstens 0,004
   (nicht monoton, bei f = 3,5 am größten), mit dem Fenster um 0,0016–0,0017 (bei 128 µm kleiner).
3. **Kernanteil und rücklaufender Fluss.** f_Kern · r_core + (1 − f_Kern) · r_mantel mit dem gemessenen
   Kernanteil ergibt 8,975 / 8,794 / 8,882·10⁻⁵ (96 µm) und stimmt mit `backward_over_forward_ref`
   auf 0,3 % überein (Referenz 8,998 / 8,816 / 8,905·10⁻⁵). Der gemessene rücklaufende Fluss
   (8,356 / 8,338 / 8,348·10⁻⁵) liegt bei 0,929 / 0,946 / 0,938 dieses Werts (128 µm: 0,930 / 0,947 /
   0,939). Der Kernanteil, den der gemessene Fluss ergäbe, ist 0,838 (96 µm) bzw. 0,836 (128 µm),
   gemessen sind 0,82–0,83; die Differenz ist 0,011–0,015.
4. Der Kernanteil in 2D wurde nicht gemessen; der Auftrag nennt aus den Fresnel-Werten rückwärts
   gerechnet 0,90 (2D) und 0,84 (3D). Der gemessene Kernanteil in 3D (0,82–0,83) liegt in der Nähe
   des 3D-Werts, nicht bei 0,90.

Rohdaten (`experiments/out/*.json`) sind nicht eingecheckt.

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

## K4: ebene Schichtstapel (B2)

λ = 1 µm, dy = dz = λ/20, Fenster 64 Punkte; ebene Wellen e^(i k_y y) mit sin θ = 0; 0,31; 0,63;
0,94 in Vakuum (0°, 18°, 39°, 70°) und ein evaneszenter Fall. Stapel: Platte n = 1,5 der Dicke
2,5 µm in Vakuum (wie wende2023 §2.1.3), Mehrschichtstapel 1,45/2,1/1,45/1,52, IP-S mit 0,4 µm
Silber. Nach 60 Iterationen stimmen r und t der tangentialen Felder in allen 30 Fällen auf
< 10⁻¹² mit der Transfer-Matrix des diskretisierten Stapels überein, die Platte zusätzlich mit der
Airy-Formel. Damit ist die unabhängig hergeleitete Zuordnung der Medien und Ausbreitungsfaktoren
bestätigt, und sie **stimmt mit Tabelle 1 von wende2023 überein** (E3).

Energieresiduum |ε| = |1 − (P_R + P_T)/P_in| der Platte über die Iterationen 1 … 8:

| Fall | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 18°, TE | 4·10⁻² | 2·10⁻³ | 2·10⁻⁴ | 4·10⁻⁶ | 2·10⁻⁷ | 2·10⁻⁸ | 4·10⁻¹⁰ | 2·10⁻¹¹ |
| 39°, TM | 1·10⁻² | 3·10⁻⁴ | 8·10⁻⁶ | 2·10⁻⁸ | 2·10⁻⁹ | 3·10⁻¹¹ | 6·10⁻¹⁴ | 2·10⁻¹⁵ |
| 70°, TE | 3·10⁻¹ | 4·10⁻³ | 3·10⁻² | 1·10⁻² | 1·10⁻³ | 5·10⁻⁴ | 3·10⁻⁴ | 7·10⁻⁵ |

Wie in wende2023 Abb. 2(e, f) fällt der Fehler etwa geometrisch, bei großen Winkeln langsamer;
drei Iterationen reichen für < 1 % außer bei streifendem Einfall.

**Grenze:** Eine Welle, die außen evaneszent und in einer hochbrechenden Schicht geführt ist
(sin θ = 1,56 im Stapel mit n = 2,1 zwischen 1,45 und 1,52), konvergiert nicht; der Fehler bleibt
nach 60 und 200 Iterationen bei etwa 10⁻². Die Iteration ist eine Neumann-Reihe der
Mehrfachreflexionen, deren Spektralradius bei Totalreflexion an beiden Seiten nicht unter 1 liegt.
Für geführte Moden (Wellenleiter, Faserkern) ist die BWPM-Iteration daher kein geeignetes Verfahren
für die Rückreflexe innerhalb der Führung.

## K5: Parabolspiegel (B2)

wende2023 §2.1.4 nachgestellt: λ = 0,78 µm, dy = dz = λ/(10·1,505), f′ = 50 µm, Silber 0,4 µm auf
IP-S, kollimierter TE-Gaußstrahl w₀ = 30 µm, eine Iteration. Maximum der Intensität des
rücklaufenden Felds auf der Achse: 49,87 µm vor dem Scheitel (λ/20: ebenfalls 49,87 µm),
wende2023: 49,85 µm. Die Abweichung von 0,13 µm gegen f′ ist nicht durch den Phasensprung des
Silbers (0,02 µm) erklärt und bleibt als Beobachtung offen. Die reflektierte Leistung beträgt nur
0,82 (λ/10) bzw. 0,77 (λ/20) der eingestrahlten, obwohl |r|² = 0,98 – siehe K9.

## K9: gekippter Silberspiegel (B2)

λ = 0,5 µm, Spiegelebene z = 15 µm + y·tan α, TE-Gaußstrahl w₀ = 8 µm senkrecht, eine Iteration:

| α | Richtung λ/10 / λ/20 / λ/40 (Soll 2α) | Leistung λ/10 / λ/20 / λ/40 | Fresnel |r|² |
| --- | --- | --- | --- |
| 5° | 6,63° / 9,99° / 10,00° | 0,857 / 0,835 / 0,778 | 0,982 |
| 10° | 19,86° / 19,96° / 20,00° | 0,732 / 0,678 / 0,586 | 0,982 |
| 20° | 39,47° / 39,99° / 39,99° | 0,485 / 0,362 / 0,279 | 0,983 |
| 30° | 56,63° (λ/10) | 0,203 (λ/10) | 0,984 |

Für Glas (n = 1,5) bei λ/10: R = 0,035 / 0,028 / 0,015 / 0,006 gegen 0,040 / 0,042 / 0,047 / 0,058
(5°/10°/20°/30°).

**Befund:** Die Richtung der Reflexion konvergiert, die **reflektierte Leistung an gekippten
Flächen nicht** – sie fällt mit feinerem Gitter weiter. Ursache ist dieselbe Näherung wie in K8,
hier stärker: Jede Treppenstufe reflektiert als z-senkrechte Fläche, und der Reflex läuft beim
Rückweg pixelweise in die schon metallischen Nachbarpixel der Treppe, wo er absorbiert wird;
mit feinerem Gitter gibt es mehr Stufen. wende2023 zeigt den Parabolspiegel nur qualitativ.
Folge: Die BWPM ist für Reflexe an **z-senkrechten** Flächen (Facetten, Schichten) quantitativ,
an gekippten Flächen nur für Richtung und Phase brauchbar (Issue zur Abhilfe).

## K10–K12: 3D (B3)

Gitter 32 × 32, dx = dy = dz = λ/16, λ = 1 µm. K11 mit 30° längs x, 45° unter den Azimuten 45° und
135° und senkrechtem Einfall, komplexe Polarisation: Transmission und Reflexion auf < 10⁻¹² nach 30
Iterationen.

**Beobachtung (Stabilität):** Mit deutlich mehr Iterationen wächst im 3D-Stapel mit n = 2,1 und
Silber der Fehler wieder (45°: 2·10⁻¹⁷ nach 30, 4·10⁻¹⁴ nach 60, 2·10⁻⁴ nach 150, divergent nach 300
Iterationen; etwa Faktor 1,3 je Iteration). Rundungsrauschen in Wellen, die in der Schicht mit
n = 2,1 gefangen sind (außen evaneszent), wird von der Iteration verstärkt – die Erweiterung der
Grenze aus K4: Dort konvergierten solche Wellen nicht, hier wachsen sie. In 2D trat das bei den
geprüften Stapeln bis 300 Iterationen nicht auf. **Folge:** Iterationen begrenzen (die Quelle
nutzt 1–10) und am Energieresiduum abbrechen; für viele Iterationen wäre ein Filter der gefangenen
Anteile nötig.

## Messlauf B3 (Rechner des Maintainers, 2026-11-05)

`experiments/b3_scaling.py`, Faserfacette mit Index-Matching-Kappe in 3D (Richards-Wolf-Spektrum,
NA 0,4, x-polarisiert, Kern rund mit 2 µm), Fenster 96 µm, 2 Iterationen, NumPy-FFT auf einem Kern
(vor der Option `--workers`). Speichermodell nachträglich mit den Materialkennungen
(1 B je Gitterpunkt):

| Abtastung | Gitter | Laufzeit | Speicher gemessen | Modell (mit Kennungen) | gespeicherte Reflexwerte | Rückreflex | Strehl auf der Facette | bester Fokus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| λ/(2,5 n) | 572² × 506 | 51 min | 0,94 GB | 0,86 GB | 1,31·10⁶ (2 Ebenen) | 4,185 % | 0,99694 | +0,26 µm |
| λ/(3,5 n) | 800² × 708 | 2 h 26 min | 2,04 GB | 1,80 GB | 2,56·10⁶ (2 Ebenen) | 4,180 % | 0,99695 | +0,12 µm |

**Befund:**

1. Rückreflex und Strehl-Verhältnis sind konvergiert (Änderung in der 5. Stelle). Der 3D-Strehl
   0,997 gleicht dem 2D-Wert (K6, 0,9972); der Unterschied zu 0,97 in wende2023 ist kein
   Dimensionseffekt (vermutlich die Definition, `docs/pilot-imc-2d.md`).
2. Die Laufzeit wächst etwa linear mit der Zahl der Gitterpunkte (Faktor 2,9 bei 2,7-mal so vielen
   Punkten); λ/(5 n) wie in der Quelle (≈ 1142² × 1012) bräuchte auf einem Kern ≈ 7 h. Der Lauf
   wurde deshalb abgebrochen. Etwa 3 s je Ebene und Durchlauf bei 800² – die FFTs je Brechzahlpaar
   dominieren; B4 bringt sie nach C++/CUDA, `--workers -1` nutzt bis dahin alle Kerne.
3. Die Reflexe liegen nur an Kappenoberfläche und Facette: zwei volle Ebenen statt N_z, also
   < 1 % des Speichers der Quelle für die Reflexe.
4. Das Modell liegt mit den Kennungen bei 88–91 % des gemessenen Spitzenspeichers; der Rest sind
   temporäre Arrays und der Python-Prozess.

**Vergleich ohne Kappe** (`--bare`, λ/(2,5 n), 2 Iterationen, `--workers -1`, 2026-11-06):
572² × 208, 19,6 min, 0,63 GB (Modell 0,57 GB), 6,5·10⁵ gespeicherte Reflexwerte (nur die
Facette, je Richtung eine Ebene), Rückreflex 3,93 % – in 3D etwas weniger als in 2D (4,23 %), weil
der runde Fokus mehr Fläche des Mantels (n = 1,4485, Fresnel 3,7 %) trifft als der Linienfokus.
Mit Kappe kommen 4,18 % von der Kappenoberfläche, aber 82 µm vor der Facette (K6).

**Threads in NumPy:** 5,6 s je Ebene (beide Iterationen) gegen 6,1 s im einkernigen Lauf mit
Kappe – die Option `--workers` bringt kaum etwas, weil nicht die FFTs, sondern die
elementweisen Schritte (Matrix × Spektrum, Maskierung) in NumPy dominieren. Der Gewinn kommt
aus dem C++-Kern (K14: einkernig ×3,7) und der GPU (B4d).

## K14: CPU-Kern (B4c)

Alle fünf Sollwertdateien in ihrer Toleranz (relativ 10⁻¹², absolut 10⁻¹³ des Maximums), dazu die
Python-Referenz auf zufälligen 3D-Volumina auf < 10⁻¹² des Eingangs. Gefunden dabei: Liegen
Gitterpunkte bis auf Rundung auf |k⊥| = n·k₀ (hier mit dx ≠ dy), löscht k_z = √(n²k₀² − k⊥²) aus;
rundete der Kern k als 2π·m/(N·d) statt wie NumPy als 2π·(m·(1/(N·d))), wichen die Ergebnisse um
2·10⁻¹⁰ ab. Der Kern rundet jetzt wie `numpy.fft.fftfreq`.

Laufzeit (Container, ein Kern, Faserfacette mit Kappe in 3D, 214² × 303, eine Iteration):
NumPy 16,1 s, C++-Kern 4,4 s (Faktor 3,7), gleiche Ergebnisse (Strehl 0,99310, Rückreflex 4,076 %).
Mit OpenMP und pocketfft-Threads skaliert der Kern zusätzlich mit der Kernzahl; Messung auf dem
Server mit `experiments/b3_scaling.py --core`.

## K15: CUDA-Kern (B4d, Rechner des Maintainers, 2026-11-06)

Umgebung: RTX 6000 Ada (48 GB), nvcc 12.9 (HPC SDK), GCC 14.3 bzw. gcc-toolset-15 als
Host-Compiler, Python 3.13 (conda). Build mit RPATH auf cuFFT; wegen der älteren libstdc++ von
conda mit `BWPM_STATIC_LIBSTDCXX=ON`.

Tests: 6 von 7 grün; `bwpm3d_tilted` wich um 10⁻⁹ ab – ebenso im CPU-Kern auf dem Server, CUDA
gegen CPU ≈ 10⁻¹⁶. Ursache war das Rundungsrauschen an streifenden Gitterpunkten (#25), jetzt
durch eine einheitliche Regel in Python, C++ und CUDA behoben (Spezifikation §3, Punkt 5).

**Bestätigungslauf** (Commit `17fce34`): mit CUDA und gcc-toolset-15 alle Tests grün (124 gesamt;
25 in den Kern- und CUDA-Tests, darunter die vier Tests aus #25 und der neue Test zur Regel);
CPU-Kern ohne CUDA mit gcc-toolset-15 18/18. Mit dem System-GCC 14 lädt das Modul in der
conda-Umgebung nicht (libstdc++ zu alt), dort ist die Regel nicht direkt geprüft; sie ist aber
durch die Schwelle 10⁻¹² unabhängig von der Rundungsreihenfolge, und der Test
`test_grazing_rule_is_independent_of_rounding` prüft genau das auf jeder Plattform.

Messung (Faserfacette, Fenster 96 µm, 2 Iterationen):

| Abtastung | Gitter | NumPy (1 Kern) | CPU-Kern (32 Kerne) | CUDA | Host / GPU-Speicher (CUDA) | Strehl | Rückreflex |
| --- | --- | --- | --- | --- | --- | --- | --- |
| λ/(2,5 n) | 572² × 506 | 51 min | 99 s | 23 s | 1,6 GB / ≤ 1,6 GB | 0,99694 | 4,185 % |
| λ/(3,5 n) | 800² × 708 | 2 h 26 min | – | 59 s | 4,1 GB | 0,99695 | 4,180 % |
| λ/(5 n) | 1142² × 1011 | (≈ 7 h, nicht gerechnet) | – | 240 s | 11,3 GB | 0,99691 | 4,170 % |
| ohne Kappe, λ/(2,5 n) | 572² × 208 | 19,6 min | – | 9,9 s | 0,7 GB / 1,3 GB | – | 3,925 % |

CUDA ist 130–150-mal schneller als NumPy und 4,3-mal schneller als der CPU-Kern auf 32 Kernen;
bei λ/(5 n) etwa 0,06 s je Ebene und Durchlauf – das Ziel aus `docs/b4-entwurf.md` (< 0,1 s) ist
ohne Bündelung der Paare erreicht. Strehl und Rückreflex wie die Referenz; bei λ/(5 n) Strehl
0,99691 (Konvergenz in der 5. Stelle). Der Host-Speicher wächst mit dem Kennungsfeld (int32 im
Kern, 4 B je Gitterpunkt), die GPU braucht wenig, weil T und R im Kernel berechnet werden.

## K16: Konformität des Adapters (B4e)

Für jedes verfügbare Backend (Python-Referenz, CPU-Kern, CUDA) über die Adapter-Schnittstelle
(`prepare`/`apply`, Ports in den Konventionen von Raytatouille):

- Schichtstapel Luft/1,5/2,1 + 0,02i/1,45/Luft, ebene Welle 0°, 30° und 45° unter den Azimuten
  45° und 135°, komplexe Polarisation, 30 Iterationen: reflektiertes und durchgelassenes
  Querfeld gegen die Transfer-Matrix auf < 10⁻¹²; Platte n = 1,5 gegen Airy auf < 10⁻¹².
- Energieerhaltung (verlustfreier Stapel, bandbegrenztes Zufallsfeld, 60 Iterationen):
  (P_R + P_T)/P_in − 1 < 10⁻¹², mit dem Poynting-Fluss der Ports (`port_power`).
- Reziprozität: asymmetrischer, absorbierender Stapel mit Silber, Vakuum auf beiden Seiten;
  Transmission von vorn und von hinten (gespiegelter Stapel) gleich auf < 10⁻¹², Reflexion
  verschieden (> 10⁻³).
- Determinismus: gleiche Struktur → gleicher Inhalts-Hash, gleiche Felder (Python und CPU bitgleich,
  CUDA in 10⁻¹⁰); andere Struktur → anderer Hash.
- Port-Konvention: Verschiebung des Eingangs um ganze Pixel verschiebt den Ausgang gleich
  (< 10⁻¹³) – kein Bezug auf den Fensterrand.

Offen als Querprobe des Fusionsplans: derselbe Schichtstapel mit FourierGrating, sobald dessen
Adapter existiert.

## K1: Grenzflächenmatrizen (B0)

Messwerte und Fälle: `docs/spezifikation.md`, Abschnitt 3.
