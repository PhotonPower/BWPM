# Spezifikation (Stufe B0)

> Grundlage: wende2023 (Volltext, `docs/quellen.md`). Gleichungsnummern „(n)“ beziehen sich auf
> diese Quelle. Konventionen nach `docs/konventionen.md`. Status: Entwurf zur Freigabe.

## 1. Aufgabe

Gegeben ist ein Volumen mit stückweise konstanter, komplexer Brechzahl ñ(x, y, z) auf einem
Gitter (dx, dy, dz) und ein einfallendes Feld in der Ebene z = z₀. Gesucht sind die vor- und
rücklaufenden Felder E⁺, E⁻ in allen Ebenen, insbesondere an den Rändern z₀ und z₁ des Volumens
(Ports). Das Verfahren ist die FPWPM (wende2022, zusammengefasst in Abschnitt 2.1.1 der Quelle),
erweitert um Reflexion und Rückwärtsausbreitung (Abschnitt 2.1.2).

## 2. Diskretisierung und Bezeichnungen

- Querschnitt (x, y) mit N_x × N_y Punkten, Schritt in z: dz. In 2D (Quelle, Abschnitt 2.1.2
  und 2.1.4) ist x unendlich ausgedehnt, N_x = 1.
- In jedem Schritt z → z + dz kommen M verschiedene Brechzahlen ñ₁ … ñ_M vor (M ist die Zahl der
  verschiedenen Werte im ganzen Volumen). Index r gehört zu ñ(r⊥, z), Index s zu ñ(r⊥, z + dz).
- **Gewichtsfunktion** (4): I_r,s(r⊥) = 1, wenn ñ(r⊥, z) = ñ_r und ñ(r⊥, z + dz) = ñ_s, sonst 0.
  Die Paare (r, s) zerlegen den Querschnitt; r ≠ s sind Grenzflächen in diesem Schritt.
- k_r = ñ_r k₀, k_z,r = √(k_r² − k⊥²) mit Im k_z,r ≥ 0 (Konventionen), k⊥² = k_x² + k_y².
- Das Feld wird nur über seine Querkomponenten E⁽²⁾ = (E_x, E_y) geführt; E_z folgt aus (5):
  E_z(k⊥) = −k⊥·E⁽²⁾(k⊥)/k_z,s für die vorlaufende Welle (für die rücklaufende mit + und k_z
  des jeweiligen Mediums, Konventionen).

## 3. Grenzflächenmatrizen und Konventionsabgleich

Die Quelle gibt für eine ebene Grenzfläche zwischen den Medien r (bei z) und s (bei z + dz) die
Transmissionsmatrix (2), (3)

T_r,s(k⊥) = 𝒫_s/(k⊥² k_r k_s) · [[t_TM k_x² Γ_r,s + t_TE k_y² k_r k_s, (t_TM Γ_r,s − t_TE k_r k_s) k_x k_y],
[(t_TM Γ_r,s − t_TE k_r k_s) k_x k_y, t_TM k_y² Γ_r,s + t_TE k_x² k_r k_s]],
Γ_r,s = k_z,r k_z,s + k_z,s k⊥²/k_z,r,

und die Reflexionsmatrix (8), (9)

R_r,s(k⊥) = 𝒫⁽ⁱ⁾/(k⊥² k_r²) · [[−r_TM k_x² Γ_r + r_TE k_y² k_r², −(r_TM Γ_r + r_TE k_r²) k_x k_y],
[−(r_TM Γ_r + r_TE k_r²) k_x k_y, −r_TM k_y² Γ_r + r_TE k_x² k_r²]],
Γ_r = k_z,r² + k⊥² = k_r².

Die Quelle nennt für t und r nur „Fresnel’s equations“ (Born & Wolf). Die Vorzeichenwahl von r_TM
ist zwischen Lehrbüchern verschieden; die Reflexionsmatrix enthält ein ausdrückliches −r_TM.
**Ergebnis des Abgleichs:** Mit den Koeffizienten der Convention A (t_TE = t_s, t_TM = t_p,
r_TE = r_s, r_TM = r_p, Byrnes Gl. (6)) sind beide Matrizen exakt der vektorielle
Fresnel-Operator auf den Querkomponenten, im globalen System, für die durchgelassene (+z) und die
reflektierte (−z) Welle. Anschaulich: Bei senkrechtem Einfall gilt in Convention A r_p = −r_s,
und der p-Vektor der reflektierten Welle hat die umgekehrte x-Komponente; das −r_TM der Quelle
macht daraus wieder r_s für beide Querkomponenten.

Geprüft in `tests/test_interface.py` (K1) gegen einen unabhängig aufgebauten Fresnel-Operator
(Ergänzung zum Volltransversalfeld, Zerlegung in s und p, Koeffizienten mit Winkeln nach Byrnes),
je 300 zufällige Richtungen bis 85°, beliebige komplexe Polarisation (Seed fest), Fehler relativ
zum einfallenden Feld:

| Grenzfläche | T | R |
| --- | --- | --- |
| Luft → Glas (1,5) | 7,9·10⁻¹⁶ | 3,1·10⁻¹⁶ |
| Glas → Luft, mit Totalreflexion | 1,3·10⁻¹⁵ | 7,5·10⁻¹⁶ |
| Luft → Silber (0,213 + 6,697i, Quelle Abschnitt 2.1.4) | 1,2·10⁻¹⁵ | 7,4·10⁻¹⁶ |
| IP-S (1,505) → Silber | 9,8·10⁻¹⁶ | 9,3·10⁻¹⁶ |
| Index-Matching 1,5 → 1,511 (Quelle Abschnitt 2.2.2) | 1,1·10⁻¹⁵ | 2,2·10⁻¹⁶ |
| gleiches Medium | 6,0·10⁻¹⁶ | 0 |

Dazu: Grenzwert k⊥ → 0 (T = t·I, R = r_s·I, stetig auf 10⁻⁹), Brewster-Winkel (r_p = 0) und die
Leistungsbilanz R + T = 1 auf einem 41 × 41-Gitter in k⊥ (Abweichung < 10⁻¹³).

**Folgerungen für die Umsetzung:**

1. Fresnel-Koeffizienten in Convention A, nicht aus einer anderen Quelle; Code
   `bwpm.interface.fresnel_coefficients`.
2. Bei k⊥ = 0 sind beide Formeln 0/0; dort gilt der Grenzwert (t_s = t_p, r_p = −r_s).
3. R vereinfacht sich wegen Γ_r = k_r² zu (1/k⊥²)·[[−r_p k_x² + r_s k_y², −(r_p + r_s) k_x k_y], …];
   der Code hält die Form der Quelle bei, damit sie nachprüfbar bleibt.
4. Totalreflexion und absorbierende Medien brauchen keine Sonderbehandlung, wenn k_z nach der
   Konvention gewählt wird.
5. **Streifender Einfall (#25):** Liegt ein Gitterpunkt bis auf Rundung auf |k⊥| = ñk₀, ist
   ñ²k₀² − |k⊥|² reines Rundungsrauschen, dessen Wert von der Reihenfolge der Rechenschritte
   abhängt (gemessen: GCC 14 gegen NumPy T₀₀ = 2 + 1,2·10⁻⁸ i gegen 1,99999999999982; Folge in 3D
   10⁻⁹). Regel in Python, C++ und CUDA gleich: Ist |ñ²k₀² − |k⊥|²| ≤ 10⁻¹²·|ñ²k₀²|, gilt die Welle
   als exakt streifend, k_z = i·10⁻¹²·k₀. Betroffen sind nur Wellen, die bis auf 10⁻⁶ streifend
   sind und keine Leistung tragen.

## 4. Vorwärtsschritt (FPWPM)

(1), (6): E⁺(r⊥, z + dz) = Σ_r Σ_s I_r,s(r⊥) · F⁻¹{ T_r,s(k⊥) F{E⁺(r⊥, z)} }, mit
𝒫_s = e^(i k_z,s dz) (Tabelle 1, Fall i = 1).

Bemerkungen zur Umsetzung:

- F{E⁺} wird einmal je Schritt gerechnet; für jedes vorkommende Paar (r, s) ein F⁻¹, maskiert mit
  I_r,s. Der Aufwand ist O(P·N log N) mit P ≤ M² vorkommenden Paaren (Quelle: O(M²[…])).
- **Näherung des Verfahrens:** Jedes Medium r wird für das ganze Querfeld als unendlich
  ausgedehnt behandelt und erst danach örtlich ausgewählt (wie WPM). Fehlerquellen sind daher
  die Treppe gekrümmter Flächen (dy, dz), die Wahl der Brechzahlstufen und Felder nahe Kanten
  von Querschnittsänderungen. Die Quelle wählt dy = dz = λ/(10|ñ_max|) (Abschnitt 2.1.4) bzw.
  λ/(5 n_max) (Abschnitt 2.2.2).
- Mehrfachreflexionen innerhalb eines Schritts werden nicht erfasst; sie entstehen über die
  Iterationen (Abschnitt 5).
- **Gemessene Grenze (K8, B1):** An gekippten Flächen stimmt die Richtung bis 45° Flächenwinkel auf
  ≤ 0,01°, die durchgelassene Leistung aber nur auf einige Prozent (bis 7 % bei 20–30°, nicht
  konvergent in dz, auch über 1). Bei 60° weicht auch die Richtung um etwa 2° ab.
- **Randmaske (K7, B1):** Randzone ≥ 20 λ, Ordnung 8 (Fehler im Inneren < 10⁻⁵ bis 30°
  Ausbreitungswinkel). Code `bwpm.fpwpm.absorbing_mask`.

## 5. Bidirektionales Verfahren (BWPM)

- **Reflexion im Vorwärtslauf** (7): E⁻(r⊥, z) = Σ_r Σ_s I_r,s · F⁻¹{ R_r,s F{E⁺(r⊥, z)} },
  gespeichert je Ebene (Tabelle 1, Fall i = 2: 𝒫 = 1).
- **Rückwärtslauf** (10), (11): Am Ende des Volumens wird die Richtung umgekehrt; rücklaufende
  Felder werden als vorlaufend in einem in z gespiegelten System behandelt. Der rücklaufende
  Eingang jedes Schritts ist die Summe aus dem weitergeführten E⁻ und dem im Vorwärtslauf
  gespeicherten Reflex an dieser Ebene; die Reflexion im Rückwärtslauf erzeugt Beiträge zu E⁺, die
  im nächsten Vorwärtslauf eingehen. Ein Hin- und Rücklauf ist eine **Iteration**. Der
  rücklaufende Anteil startet mit null.
- **Ausbreitungsfaktoren** (Tabelle 1): i = 1 durchgelassen vorwärts e^(i k_z,s dz); i = 2
  Speichern des Reflexes vorwärts 1; i = 3 durchgelassen rückwärts e^(i k_z,r dz); i = 4
  Reflex rückwärts e^(2i k_z,r dz) (hin und zurück durch ñ(z)). Sie sind so gewählt, dass jedes Feld
  am rechten Rand seines Pixels liegt (Abb. 1(b) der Quelle). Die genaue Zuordnung von r und s im
  Rückwärtslauf legt B2 fest und prüft sie an der Platte (K4).
- **Gemessene Grenzen (B2):** (a) Für Wellen, die durch Totalreflexion in einer Schicht gefangen
  sind, konvergiert die Iteration nicht (Neumann-Reihe mit Spektralradius ≥ 1, K4). (b) An
  gekippten reflektierenden Flächen konvergiert die Richtung des Reflexes, die reflektierte
  Leistung nicht; sie liegt weit unter Fresnel (Silber 10°: 0,59–0,73 statt 0,98, K9).
- **Stabilität (B3, K11):** Bei vielen Iterationen kann Rundungsrauschen in gefangenen Wellen
  wachsen (3D-Stapel mit n = 2,1 und Silber: Faktor ≈ 1,3 je Iteration). Iterationen begrenzen.
- **Streifende Wellen:** k_z = 0 auf dem Gitter wird durch i·10⁻¹²·k₀ ersetzt (sonst 0/0 in
  Fresnel; tritt in 3D bei |k⊥| = n k₀ genau auf).
- **Abbruch:** Die Quelle misst die Konvergenz an ε = 1 − |R + T| (Abschnitt 2.1.3) und empfiehlt
  drei Iterationen; bei stark reflektierenden Spiegeln reicht eine (Abschnitt 2.1.4). Für die
  Umsetzung: ε als Diagnose je Iteration, Abbruch bei einer Schwelle oder einer Höchstzahl.
- **Aufwand:** Ein Schritt etwa doppelt so teuer wie bei der FPWPM, eine Iteration etwa viermal
  (Quelle, Abschnitt 2.1.2). **Speicher:** gespeicherte Reflexe N_x N_y N_z Vektoren in der
  allgemeinen Form; die Quelle rechnet deshalb in 2D. Reflexe entstehen nur an Grenzflächen
  (I_r,s mit r ≠ s); B3 speichert nur diese.

## 6. Ports

- Eingang: Querfeld E⁺ in z₀ (vorlaufend) und optional E⁻ in z₁ (rücklaufend), Gitter und
  Wellenlänge, Medium am Port.
- Ausgang: E⁺ in z₁ (durchgelassen), E⁻ in z₀ (reflektiert), optional Felder in Zwischenebenen.
- Leistung durch eine Ebene über den Poynting-Fluss im Frequenzraum (vorlaufend und rücklaufend
  getrennt, mit E_z nach Abschnitt 2); Grundlage für ε und für die Kopplung (B4).

## 7. Entscheidungen und offene Punkte

Entschieden vom Maintainer am 2026-11-05 (Issues #1–#3, `docs/roadmap.md`, E2–E4):

1. **Randbedingung quer (E2, #1):** absorbierende Maske nach jedem Schritt (glatte Funktion,
   Mitte 1, Rand 0, etwa super-gaußförmig) als Standard; großes Fenster mit Nullrand als Referenz
   in den Tests. Eine Absorberschicht über die Brechzahl scheidet aus, weil die BWPM an jeder Stufe
   der Rampe einen Reflex speichern und zurücktragen würde und M wüchse. B1 misst die nötige Breite
   der Randzone an einem gekippten Gaußstrahl, der das Fenster verlässt.
2. **Zuordnung r, s und Lage der Felder im Rückwärtslauf (E3, #2):** unabhängig herleiten und
   Tabelle 1 daran prüfen. Für ebene Wellen und ebene Schichten ist die BWPM nach unendlich vielen
   Iterationen exakt die Transfer-Matrix-Lösung des treppenförmig diskretisierten Stapels; K4
   prüft deshalb auf Maschinengenauigkeit (Platte und Mehrschichtstapel, schräg, TE und TM, auch
   Silber), zusätzlich gegen die Airy-Formel. **Ergebnis (B2):** Die Herleitung (Kopf von
   `bwpm.bidirectional`) stimmt mit Tabelle 1 überein; K4 < 10⁻¹².
3. **Brechzahlstufen (E4, #3):** Pixelmittelpunkt als Standard; Flächenanteil als Option nur für
   dielektrische Übergänge, auf wenige Stufen gerundet. B1 misst an einem gekippten Glaskeil
   (Flächenwinkel 0–60°) Richtung und Amplitude gegen Fresnel an der gekippten Fläche; daraus folgt
   die Gültigkeitsgrenze im Flächenwinkel, die auch der Fusionsplaner braucht.

Offen, zu klären in den genannten Stufen:

4. **Speicher in 3D (B3):** umgesetzt in `bwpm.volume3d`: Reflexquellen nur an Pixeln, an denen
   sich die Brechzahl zwischen zwei Scheiben ändert, als (Pixelindex, Wert); Zwischenebenen nur auf
   Anforderung. Modell `memory_estimate`: Operatoren (Arten × Paare × N_x N_y × 64 B), etwa acht
   Arbeitsebenen (× 32 B) und die Quellen (× 40 B). Messungen auf dem Rechner des Maintainers mit
   `experiments/b3_scaling.py`.
5. **Leistungsnormierung der Ports (B4):** Kopplung an Raytatouille nach Fieldtracing ADR 0004.

## 8. Abnahme B0

- K1 grün (`tests/test_interface.py`), Werte in Abschnitt 3.
- Spezifikation vom Maintainer freigegeben; Entscheidungen E2–E4 eingetragen.

## 9. Erste Born-Näherung als Prüfweg (B6, K17, #29)

Zweck: Genauigkeit des Vorwärtsschritts bei Streuung in große Winkel an schwachen, feinen
Brechzahländerungen prüfen, unabhängig vom Verfahren. Code `bwpm.born.born_forward_2d`.

- **Aufbau (2D, k_x = 0):** Pixel (j, n) ist die Zelle z ∈ [j·dz, (j + 1)·dz], |y − y_n| ≤ dy/2
  mit konstantem ñ; Hintergrund ñ₀, periodisch in y. Ebene Welle e^(i k z), k = ñ₀k₀, senkrecht;
  TE: E = E_x, TM: E = E_y.
- **Formel (bornwolf1999, chew1995; nur bibliografisch, daher ohne Gleichungsnummern):** mit
  V = k₀²(ñ² − ñ₀²) und k_z = sqrt(k² − k_y²), Im k_z ≥ 0, gilt hinter allen Streuern

  E_s,m(z) = i/(2k_z) · Σ_j V̂_j,m · ∫_Scheibe j e^(i k_z (z − z′)) e^(i k z′) dz′,

  V̂_j,m = DFT der Pixelwerte · sinc(k_y dy/2) (exakter Fourier-Koeffizient der Treppe), das
  z-Integral je Scheibe geschlossen. TM: Die dyadische Green-Funktion liefert für die
  Querkomponente E_y den Faktor 1 − k_y²/k² = k_z²/k² (= cos²θ).
- **Prüfung der Referenz:** gegen die Airy-Formel der homogenen Platte und gegen eine RCWA des
  Lamellengitters (moharam1995, Faktorisierung nach li1996; `tests/rcwa_reference.py`); die
  Abweichung ist von erster Ordnung in δn (≈ 10⁻⁴ bei δn = 10⁻⁵).
- **Vergleichsebene:** hinter der Schicht im Hintergrundmedium. In der Schicht fehlt dem
  vorlaufenden Feld der Reflex; es weicht dort um den Fresnel-Faktor 1 − δn/(2ñ₀) vom
  Gesamtfeld ab, was bei der Born-Amplitude (∝ δn) ein Fehler erster Ordnung ist.
- **Gemessene Grenze (K17):** Paraxial stimmt die FPWPM mit Born überein (< 10⁻³ bei 0,8°). Bei
  großen Winkeln weicht sie ab, konvergiert in dy = dz und unabhängig von der Schichtdicke:
  TE-Amplitude 0,93 / 0,80 / 0,58 bei 29° / 46° / 63°; für TM ist die Querkomponente E_y gleich
  der TE-Amplitude, der Vektorfaktor cos²θ fehlt (E_y zu groß um 1/cos²θ). Die Phase konvergiert
  linear in der Schrittweite gegen Born. Werte in `docs/testfaelle.md`, K17.
- **Ursache (Analyse erster Ordnung in δn, Lamellengitter, senkrechter Einfall):** Mit
  k_z,s ≈ k_z + k₀²ñ₀δn_s/k_z wird der Schritt (wende2023 Gl. (1), r = s im Inneren) zu
  E(z + dz) = 𝒫E + i k₀²ñ₀ dz · δn(y) · F⁻¹{e^(i k_z dz) Ê(q)/k_z(q)}. Die Quelle der Streuung
  trägt also den Schrägfaktor 1/k_z(q) der **einfallenden** Welle q; in Born steht 1/k_z(k_y) der
  **gestreuten** Welle. Für senkrechten Einfall fehlt damit im Volumen der Faktor k/k_z = 1/cos θ.
  Dazu kommen die Fresnel-Faktoren an Ein- und Austritt (bei k_y = 0: 1 ∓ δn/(2ñ₀), pixelweise
  maskiert), die zur Ordnung +1 den Beitrag cos θ (1 − cos θ)/2 relativ zu Born liefern. Zusammen:

  TE: FPWPM/Born = cos θ · (3 − cos θ)/2,   TM: FPWPM/Born = (3 − cos θ)/(2 cos θ).

  Für TM verhält sich E_y wie E_x, weil im Inneren T = 1 ist und an Ein- und Austritt bei
  senkrechtem Einfall TE und TM gleich sind; die Kopplung über die Flächen y = const (stetiges
  D_y) kommt im Schritt nicht vor. Die Formeln treffen die Rechnung nach Extrapolation im Schritt
  auf ≤ 5·10⁻⁴ bis 63° (Test `test_fpwpm_deviation_matches_first_order_analysis`). Die Abweichung
  bleibt unter 1 % bis θ ≈ 11° (TE) bzw. 6,6° (TM).
- **3D (K17):** `bwpm.born.born_forward_3d` setzt S mit dem Querblock der dyadischen
  Green-Funktion, E⊥ = S·(I − k⊥k⊥ᵀ/k²)·ê, um (Eigenwerte k_z²/k² längs k⊥, 1 quer; für ein in x
  oder y homogenes Gitter gleich dem 2D-Ergebnis). Die FPWPM in 3D ist in erster Ordnung
  komponentenweise skalar: E⊥ = r(θ)·S·ê mit r = cos θ (3 − cos θ)/2, gestreutes E_y bei
  x-Polarisation ≈ 0. Gegen das volle Feld (mit E_z nach Gl. (5)) und vektorielles Born folgt
  azimutal gemittelt r²/cos²θ (s-Anteil r² zu klein, p-Anteil r²/cos⁴θ zu groß); gemessen an
  einer dünnen Zufallsschicht auf 0,1–2 %.
- **Extinktion (K17):** Die FPWPM erhält Σ|E⊥|², nicht den Fluss. Der Verlust des kohärenten
  Feldes ist deshalb die gestreute Quernorm, im Ensemble μs,FPWPM = ∫_vorwärts σ r²/cos θ dΩ
  (σ skalares Born), nicht das vektorielle μs.
- **Entscheidung E9 (2026-11-19):** Gültigkeitsgrenze festhalten, das Verfahren bleibt wie
  wende2023. Streuung an feinen Brechzahländerungen ist quantitativ (Amplitude < 1 %) nur bis
  θ ≈ 11° (s) bzw. 6,6° (p). Eine Erweiterung für schwache Fluktuationen δn (Born-konsistentes
  Inkrement im Wirtsmedium, WPM für die Struktur) ist als Vorschlag in #37 beschrieben.
- **Entscheidung E10 (2026-11-19):** Variante B3 aus #37 angenommen als Stufe B7. Formeln,
  Formeln, Quelle und Regeln für streifende Wellen in Abschnitt 10.

## 10. Schwache Fluktuationen mit Born-Inkrement (B7, E10, #37)

Zweck: Streuung an schwachen, feinen Brechzahlschwankungen δn bei großen Winkeln richtig
wiedergeben (Abhilfe zu E9), ohne den WPM-Schritt für die Struktur zu ändern. Prüfung K18.

- **Zerlegung:** ñ(r) = ñ_S(r) + δn(r). Die Struktur ñ_S (Brechzahlstufen, Grenzflächen) läuft
  unverändert durch den WPM-Schritt (Abschnitt 4). δn ist reell oder komplex, je Voxel gegeben und
  darf in der ersten Ausbaustufe nur in Pixeln ≠ 0 sein, deren Struktur-Brechzahl gleich dem
  Wirtsmedium ñ₀ der Scheibe ist (je Scheibe genau ein Wirtsmedium; sonst Fehler).
- **Schritt:** mit k = ñ₀k₀, k_z = sqrt(k² − |k⊥|²), V = k₀²((ñ₀ + δn)² − ñ₀²) und dem vollen
  Feld in der Scheibenmitte, im Wirtsmedium um dz/2 weiterpropagiert,
  E_m = F⁻¹{e^(i k_z dz/2) (Ê⊥, Ê_z)}, Ê_z = −k⊥·Ê⊥/k_z (wende2023 Gl. (5)):

  E⊥(z + dz) = WPM-Schritt(E⊥(z)) + F⁻¹{ e^(i k_z dz/2) · i dz/(2k_z) · [P(k)·F{V·E_m}]⊥ },

  P(k) = I − k kᵀ/k², k = (k_x, k_y, k_z) (Projektor der dyadischen Green-Funktion, chew1995;
  wie `bwpm.born.born_forward_3d`, dort ohne E_z des einfallenden Feldes). Ausgeschrieben für die
  Querkomponenten: [P·Ŝ]⊥ = Ŝ⊥ − k⊥(k⊥·Ŝ⊥ + k_z Ŝ_z)/k², Ŝ = F{V·E}. In 2D (k_x = 0): TE nur
  E_x mit Faktor 1; TM E_y mit (1 − k_y²/k²) Ŝ_y − k_y k_z Ŝ_z/k².
- **Mittelpunktregel:** chen2020 legt die Felder in die Scheibenmitten (Algorithmus 1); das ist
  die Mittelpunktregel für das z-Integral in Gl. (4). Die Felder der WPM liegen am rechten
  Scheibenrand (Abschnitt 2), daher die halben Schritte. Mit dem Feld am Scheibenanfang und
  e^(i k_z dz) wäre der Fehler von erster Ordnung (Phase (k − k_z)dz/2 je Scheibe, gemessen 6 %
  bei 63° und dz = 0,02 µm), mit Mittelpunkt von zweiter (0,2 %).
- **Quelle (chen2020, Volltext geprüft):** skalares Multi-Layer-Born, Gl. (1) (Streupotenzial
  V = k₀²(n_b² − n²), Vorzeichen umgekehrt zu hier), Gl. (3) (Rekursion
  U^(n+1) = U^n + U_s^n), Gl. (5) (Ausbreitung), Gl. (6) (Spektrum der Green-Funktion
  G̃(u, z) = −i P(u, |z|)/(4π η(u)), η = sqrt((n_b/λ)² − |u|²)) und Gl. (9) mit sinc ≈ 1 für
  kleines dz wie in Algorithmus 1. Umrechnung: u in 1/Länge ohne 2π, also k_z = 2πη und
  1/(4πη) = 1/(2k_z); G = e^(i k₀ n_b r)/(4πr) gehört zu e^(−iωt) wie hier. Die vektorielle
  Erweiterung (Projektor P) ist nicht aus chen2020, sondern aus der dyadischen Green-Funktion
  (chew1995, bibliografisch) und durch K17 (`born_forward_3d`) geprüft.
- **Streifende und evaneszente Wellen:** Das Inkrement gilt nur für propagierende Zellen
  (|k⊥| < k); evaneszente Zellen und |k⊥| = k bekommen kein Inkrement (sie erreichen die
  Vergleichsebene nicht, 1/k_z wäre dort unbeschränkt bzw. singulär). Ob das bei feinen δn
  nahe dem Grenzwinkel die Leistung verfälscht, misst K18.
- **Ordnung:** Gegenüber der ersten Born-Näherung der ganzen Schicht (K17-Sollwert) bleiben der
  Fehler der Mittelpunktregel (sinc in chen2020 Gl. (9)) und die Abtastung des Potenzials (ohne
  den Faktor sinc(k_y dy/2) der Pixelzelle), beide zweiter Ordnung in dz bzw. dy. Je Scheibe
  wächst die Norm um O((k₀δn dz)²) (Fluss erster Ordnung in dz). Mehrfachstreuung zwischen
  Scheiben ist enthalten (chen2020), Rückstreuung nicht (E10): Gegen eine exakte Rechnung bleibt
  ein Fehler der Ordnung δn/ñ₀ (fehlende Fresnel-Reflexe an Ein- und Austritt) und im Mittelwert
  der Phase k₀ δn²/(2ñ₀)·L (der Term −V²/(8k³) der Ausbreitungskonstante fehlt in der
  Born-Reihe ohne Rückwärtsanteil).
- **Umsetzung 2D:** `bwpm.fpwpm.Propagator2D(..., fluctuation=δn)`, δn der Form (nz, ny); das
  Inkrement wird vor der Maske addiert. Gemessen in K18.
- **δn = 0:** Ohne δn wird das Inkrement nicht berechnet; das Ergebnis ist bitgleich zu heute
  (K1–K16, `reference/`).
