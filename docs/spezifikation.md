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

## 7. Offene Punkte

Als Issues mit Label `frage` vor bzw. in den genannten Stufen:

1. **Randbedingung quer (B1):** Die FFT ist periodisch. Die Quelle spricht von „perfectly
   transparent boundary conditions“ ohne Verfahren. Optionen: Nullauffüllung mit genügend Rand,
   absorbierende Randschicht (komplexe Brechzahl mit Profil), PML in der Ausbreitungsrichtung
   nicht nötig. Entscheidung mit Messung des Randfehlers.
2. **Zuordnung r, s und Lage der Felder im Rückwärtslauf (B2):** aus Tabelle 1 und Abb. 1(b);
   Prüfung an der Platte gegen die Airy-Formel (K4).
3. **Brechzahlstufen (B1):** Wie werden gekrümmte Flächen auf das Gitter gebracht (Mittelpunkt,
   Flächenanteil)? Die Quelle erzeugt die Matrix aus einem Optikdesign (Abschnitt 2.1.4) ohne
   Angabe der Regel.
4. **Speicher in 3D (B3):** nur Grenzflächen speichern; Aufwand und Speicher messen.
5. **Leistungsnormierung der Ports (B4):** Kopplung an Raytatouille nach Fieldtracing ADR 0004.

## 8. Abnahme B0

- K1 grün (`tests/test_interface.py`), Werte in Abschnitt 3.
- Spezifikation vom Maintainer freigegeben; offene Punkte als Issues angelegt.
