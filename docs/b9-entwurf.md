# Entwurf B9: Stabiler, flusserhaltender Born-Schritt

> Freigegeben (2026-11-19): G1–G5 wie vorgeschlagen angenommen. Anlass: Instabilität des expliziten Born-Inkrements bei
> σ_n ≈ 0,06 (Scattuccino S9b-8, #29, Entscheidung E13). Grundlage: §10, §11, K18, K19.

## 1. Befund

Der Born-Schritt aus B7/B8 ist explizit: E_aus = P·E_ein + B·E_ein. Für reelles δn ist der
Streuoperator schiefhermitesch bezüglich des Flusses Σ k_z|E|², mit E dem vollen Feld einschließlich
E_z. Deshalb erhält der Schritt den Fluss nur bis zur ersten Ordnung. Je Scheibe wächst er um etwa
1 + Γ², mit Γ ≈ k₀ñ₀δn·dz/k_z. Nahe dem Grenzwinkel ist Γ groß.

Messung am Medium aus #29 (96², l_c = 0,5 µm, D = 4, σ_n = 0,0572, 6 µm, μs·z ≈ 1,1). Der Fluss
ist Σ k_z(|E⊥|² + |E_z|²) relativ zum Eingang:

| Schritt | dz = 0,1 µm | dz = 0,05 µm | Anteil \|E⊥\|² bei cos θ < 0,3 (0,1 µm) |
| --- | --- | --- | --- |
| explizit (B7/B8) | 1,187 | 1,071 | 0,50 |
| Cayley, 3 Fixpunktschritte | 0,999993 | 0,999999 | 0,042 |
| Cayley, GMRES bis 10⁻¹² | 1,000000 | 1,000000 | 0,042 |

Der Flussverlust von 3,5 %, den ich in #29 für den ersten Versuch genannt hatte, war ein Messfehler:
Das Flussmaß ließ E_z weg.

## 2. Formel

Strang-Aufteilung mit der Cayley-Form in der Scheibenmitte:

E_aus = P_½ · (1 − G/2)⁻¹(1 + G/2) · P_½ · E_ein,

- P_½ = F⁻¹ e^(i k_z dz/2) F im Wirt ñ₀;
- G = F⁻¹{ i dz/(2k_z) · [P(k)·F{V·e}]⊥ }, mit V = 2ñ₀k₀²δn (E11) und e = (E⊥, E_z),
  E_z = −k⊥·Ê⊥/k_z (wende2023 Gl. (5));
- nur propagierende, nicht streifende Wellen (§10).

**Warum der Fluss erhalten bleibt:** Im vollen Vektor ist G = i·W·P(k)VP(k) mit W = dz/(2k_z)
diagonal und P(k)VP(k) hermitesch. G ist also schiefhermitesch bezüglich K = diag(k_z).
Die Cayley-Form ist dann K-unitär, und P_½ ist unitär und vertauscht mit K. Der Fluss bleibt für
reelles δn exakt erhalten, nach Lösung des Gleichungssystems.

**Bezug zu B7/B8:** In erster Ordnung in G ist das der Mittelpunkt-Schritt aus B7. Die Werte von
K18 bleiben deshalb bis auf die zweite Ordnung gleich. Die Gauß-Punkte aus F1 sind vorwärts nicht
nötig, weil dort die Phasenfehlanpassung klein ist; für die Rückstreuung (R1, R3) bleiben sie.

## 3. Entscheidungen

| Nr. | Frage | Vorschlag | Begründung |
| --- | --- | --- | --- |
| G1 | Form des Vorwärtsschritts | **Strang-Cayley wie oben**, ersetzt das explizite Inkrement (vorwärts und R2 rückwärts) | Fluss exakt erhalten, stabil für alle Γ; erste Ordnung unverändert |
| G2 | Lösung von (1 − G/2)y = (1 + G/2)x | **GMRES** (matrixfrei, Neustart 30) bis zur relativen Toleranz 10⁻¹⁰, Startwert x + G·x | Fixpunktiteration: 18 Schritte bei σ_n = 0,057, 32 bei 0,1, keine volle Konvergenz bei 0,2; GMRES: 8–10 Schritte bei 0,057. Im Kern GMRES von Hand (klein, ohne Abhängigkeit) |
| G3 | Rückstreuung R1, R3 | bleiben explizit als Quellen mit Gauß-Punkten (§11) | Kopplung zwischen den Richtungen ist von erster Ordnung klein; die Energiebilanz der BWPM bleibt bis auf O(δn²) |
| G4 | Komplexes δn (Absorption) | gleiche Formel; Fluss nimmt dann ab, wie physikalisch erwartet | G ist dann nicht schiefhermitesch; die Cayley-Form bleibt kontraktiv für Im δn ≥ 0 |
| G5 | Sollwertdateien | `fpwpm2d_fluctuation`, `bwpm3d_fluctuation` neu erzeugen | Änderung in zweiter Ordnung |

**Kosten:** je Scheibe mit δn etwa 1 + m Anwendungen von G, mit m der Zahl der GMRES-Schritte
(8–10 bei σ_n = 0,06). Eine Anwendung von G kostet etwa so viel wie das heutige Inkrement mit einem
Gauß-Punkt. Gegenüber B7 wird eine Scheibe mit δn also etwa 5- bis 10-mal teurer. Gegenüber
Brechzahlstufen bleibt der Ansatz bei breitem Spektrum günstiger. Das misst B9b.

## 4. Prüfweg (K20)

1. **Flusserhaltung** für reelles δn: Medium aus #29 (σ_n = 0,057 und 0,1, dz = 0,1 und 0,05 µm,
   6 µm). Fluss auf ≤ 10⁻⁸, nur von der GMRES-Toleranz begrenzt.
2. **Stabilität:** Bei σ_n = 0,057 und dz = 0,1 µm bleibt der Anteil nahe dem Grenzwinkel
   beschränkt und hängt nicht von dz ab (0,1 und 0,05 µm innerhalb der Streuung).
3. **Absorption:** komplexes δn mit Im δn > 0. Der Fluss nimmt ab, und zwar um die erwartete
   Dämpfung in erster Ordnung.
4. **K17–K19 unverändert grün.** Die Messwerte von K18 und K19 werden neu notiert.
5. **Ensemble gegen vektorielles Born:** Winkelverteilung des Flusses in einer dünnen Schicht
   (einfache Streuung) wie K18 Punkt 5, jetzt über den Fluss.
6. **Kern** gegen Python auf 10⁻¹⁰, mit der GMRES-Toleranz als Grenze. CUDA gegen CPU.

## 5. Ablauf

- **B9a:** Python-Referenz und K20.
- **B9b:** CPU- und CUDA-Kern, neue Sollwertdateien, Laufzeit.
