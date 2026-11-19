# Auftrag für einen Agenten: Drift des Reflexionsgrads und Kernanteil (B5a, Diagnose)

> Für denselben Agenten und Rechner wie `docs/agent-server-b5a.md`. Es gelten `AGENTS.md` und
> `docs/konventionen.md`. Zusammenhang: Abschnitt „Messlauf B5a“ in `docs/testfaelle.md`, Issue #30.

## Ziel

Zwei Befunde aus dem Messlauf B5a einordnen, ohne etwas anzupassen:

1. **Drift:** Der Reflexionsgrad `reflectance` liegt bei f = 2,5 / 3,5 / 5 bei 4,1578 / 4,1525 /
   4,1426 %, die Referenz konstant bei 4,1584 %. Mit feinerem Gitter wächst also der Abstand
   (−0,016 / −0,14 / −0,38 %). Gefragt ist, ob er von der Fensterbreite abhängt.
2. **Faktor 1,5 zu 2D:** `backward_over_forward_at_facet` ist in 3D 8,35·10⁻⁵, in 2D (K6)
   5,6·10⁻⁵. Hypothese, bisher nicht geprüft: In 3D trifft weniger Fluss den Kern (Fresnel Kappe →
   Kern 1,3·10⁻⁵, → Mantel 4,5·10⁻⁴). Das Skript gibt dafür jetzt `core_flux_fraction` aus
   (Anteil des Vorwärtsflusses innerhalb r ≤ 1 µm, 1 µm vor der Facette). Aus den Fresnel-Werten
   folgt rückwärts gerechnet etwa 0,84 (3D) bzw. 0,90 (2D). Ein grober Container-Lauf (λ/(1,2 n))
   ergab 0,825.

Der Code steht auf dem Branch `b5a-pilot-3d-analysis` (nur die Ausgabe `core_flux_fraction` ist neu).

## Regeln

Wie in `docs/agent-server-b5a.md`: auf dem Branch arbeiten, nicht mergen, nichts an Tests,
Sollwerten, Kern oder `src/bwpm/*.py` ändern, keine Pfade, Hostnamen oder Tokens, bei Unklarheit
anhalten und berichten. Es werden nur Zahlen berichtet; Deutungen sind nicht Teil dieses Auftrags.

## Schritt 1: Stand prüfen

```bash
cd ~/BWPM && git fetch && git checkout b5a-pilot-3d-analysis && git pull
python -c "from bwpm import fast; print(fast.CORE_AVAILABLE, fast.cuda_available())"   # True True
```

## Schritt 2: Läufe

Alle mit Kappe, 2 Iterationen, CUDA-Kern, im Ordner `experiments`. Die JSON-Datei in
`experiments/out/` wird je Aufruf überschrieben: nach jedem Aufruf umbenennen, nicht einchecken.

```bash
for w in 96 128; do
  python b3_scaling.py --fractions 2.5 3.5 --window-um $w --iterations 2 --core --backend cuda
done
python b3_scaling.py --fractions 5 --window-um 96 --iterations 2 --core --backend cuda
```

Der letzte Lauf ist nur für den Kernanteil bei der feinsten Abtastung (etwa 4,5 min, 11 GB).
Optional, falls Zeit und Speicher reichen: `--fractions 5 --window-um 128` (etwa 20 GB GPU-Speicher,
mehrere Minuten). Bei Speichermangel überspringen und das vermerken.

## Schritt 3: Berichten

Eine Tabelle je Fensterbreite, Spalten f = 2,5, 3,5 (und 5 bei 96 µm), mit:

| Feld | Hinweis |
| --- | --- |
| `reflectance`, `reflectance_ref`, Verhältnis − 1 | wie im Messlauf B5a |
| `reflex_focus_um` | wie im Messlauf B5a |
| `backward_over_forward_at_facet`, `backward_over_forward_ref` | wie im Messlauf B5a |
| `core_flux_fraction` | neu |
| `n`, `nz`, `seconds`, `peak_rss_mb` | Größe und Aufwand |

Zusätzlich berichten (nur Zahlen):

1. Ändert sich der Abstand `reflectance` zu `reflectance_ref` zwischen 96 und 128 µm bei gleichem f?
   Ändert sich die Richtung des Drifts mit f?
2. `core_flux_fraction` bei f = 2,5, 3,5 und 5. Ändert er sich mit f oder mit dem Fenster?
3. Der Wert, den `core_flux_fraction` zusammen mit den Fresnel-Werten für den rücklaufenden Fluss
   ergibt, gegen den gemessenen Wert. Rechnung: r_core = ((1,5 − 1,511)/(1,5 + 1,511))²,
   r_mantel = ((1,4485 − 1,511)/(1,4485 + 1,511))², Erwartung f·r_core + (1 − f)·r_mantel. Das
   Skript gibt dasselbe als `backward_over_forward_ref` aus; nur prüfen, ob es mit dem Kernanteil
   übereinstimmt.
4. Umgebung wie im Auftrag B5a (GPU, nvcc, GCC, Python).

Die Tabellen als Kommentar am zugehörigen Issue und als Abschnitt „Diagnose B5a“ in
`docs/testfaelle.md` (unter „Messlauf B5a“) auf dem Branch einchecken. Sonst nichts ändern.

## Nicht Teil dieses Auftrags

Eine Ursache des Drifts suchen oder beheben, Änderungen an der Auswertung oder am Kern, die
Kopplung mit Raytatouille (B5b).
