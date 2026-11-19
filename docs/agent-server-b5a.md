# Auftrag für einen Agenten: 3D-Pilot (K6 in 3D, B5a) auf dem Server messen

> Für einen Programmieragenten auf dem Rechner des Maintainers (RTX 6000 Ada, viel Hauptspeicher,
> CUDA-Kern aus B4d gebaut). Es gelten `AGENTS.md` und `docs/konventionen.md`. Zusammenhang:
> `docs/roadmap.md` (B5), `docs/pilot-imc-2d.md`, `docs/testfaelle.md` (K6, K6 in 3D, K15).

## Ziel

Die Faserfacette mit Index-Matching-Kappe in 3D in der Auflösung der Quelle messen und für K6 in 3D
vier Größen berichten, die bisher fehlen: Reflexionsgrad gegen die Fresnel-Referenz, Fokus des
Reflexes, rücklaufender Fluss an der Facette und Konvergenz in der Abtastung und in den Iterationen.
Der Code steht auf dem Branch `b5a-pilot-3d-analysis` (`bwpm.pilot_imc3d`, erweitertes
`experiments/b3_scaling.py`, Tests `tests/test_pilot_imc3d.py`).

## Regeln

- Arbeite auf dem Branch `b5a-pilot-3d-analysis`. Nicht mergen, nicht auf `main` pushen.
- Verboten: Tests, Toleranzen, `reference/*.result.json`, `src/bwpm/*.py` ändern, damit ein Wert
  passt. Weicht ein Messwert von der Referenz ab, ist das ein Befund und kein Grund, die Referenz
  oder die Auswertung anzupassen. Fehler in der Auswertung selbst (nachweislich, mit Begründung) als
  eigener Commit, englisch, Imperativ.
- Keine Pfade, Hostnamen oder Tokens in Dateien oder Commits.
- Unklarheit oder nötige Entscheidung: anhalten und berichten.

## Schritt 1: Stand prüfen

```bash
cd ~/BWPM && git fetch && git checkout b5a-pilot-3d-analysis && git pull
python -c "from bwpm import fast; print(fast.CORE_AVAILABLE, fast.cuda_available())"   # True True
pytest -q tests/test_pilot_imc3d.py tests/test_core_cuda.py
```

Erwartet: beide Flags `True`, Tests grün. Ist der Kern ohne CUDA gebaut, wie in
`docs/agent-server-b4d.md` neu bauen.

## Schritt 2: Messläufe

Alle Läufe mit Kappe, Fenster 96 µm, CUDA-Kern; je Lauf eine JSON-Zeile. Die Datei
`experiments/out/b3_scaling_imc_cuda.json` wird je Aufruf überschrieben: nach jedem Lauf umbenennen
(z. B. `..._f2.5_it2.json`) und nicht einchecken.

```bash
cd experiments
python b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2 --core --backend cuda
python b3_scaling.py --fractions 2.5 --window-um 96 --iterations 1 --core --backend cuda
python b3_scaling.py --fractions 2.5 --window-um 96 --iterations 4 --core --backend cuda
```

Der erste Aufruf liefert die Konvergenz in der Abtastung (K15 meldet für f = 5 etwa 4 min und
11,3 GB GPU-Speicher), die beiden anderen die Konvergenz in der Zahl der Iterationen.

## Schritt 3: Berichten

Eine Tabelle je Lauf mit den Feldern der JSON-Zeile:

| Feld | Bedeutung | Erwartung aus dem groben Lauf (siehe unten) |
| --- | --- | --- |
| `reflectance` | reflektierter Fluss / einfallender Fluss am Eingang | 4,15 % |
| `reflectance_ref` | flussgewichtete Fresnel-Referenz Luft → Kappe | 4,16 %; beide stimmen auf unter 1 % überein |
| `reflex_focus_um` | Fokus des Reflexes, Facette bei 0 | −82,4 µm |
| `reflex_focus_ref_um` | Spiegelbild des gebauten Luftfokus | −82,42 µm |
| `backward_over_forward_at_facet` | rücklaufender / vorlaufender Fluss 1 µm vor der Facette | einige 10⁻⁵ |
| `backward_over_forward_ref` | grobe Referenz (Fresnel senkrecht, Kern und Mantel nach Fluss) | wie oben, Näherung |
| `incident_focus_check_um` | Fokus des einfallenden Feldes minus gebauter Fokus | < 10⁻⁶ µm (Test der Auswertung) |
| `strehl`, `best_focus_um` | wie K15 | 0,9969; Fokus nahe der Facette |
| `seconds`, `peak_rss_mb` | Laufzeit, Hauptspeicher | – |

Die Erwartungen stammen aus einem Lauf im Container mit λ/(1,2 n) (274² × 243). Das Gitter ist dafür
zu grob, sie sind **keine Abnahmekriterien**. Abnahme legt der Maintainer nach den Serverläufen fest.

Zusätzlich berichten:

1. Ändern sich Reflexionsgrad, Reflexfokus und rücklaufender Fluss zwischen f = 2,5, 3,5 und 5
   (Stelle, an der sich die Werte noch ändern)?
2. Ändern sie sich zwischen 1, 2 und 4 Iterationen?
3. Weicht `backward_over_forward_at_facet` merklich von der 2D-Zahl 5,6·10⁻⁵ (K6) ab? Bitte nur die
   Zahlen nennen. Eine Erklärung (Kern-Anteil des Flusses, Winkelabhängigkeit der Fresnel-Koeffizienten)
   ist erst nach einem Test belastbar, den dieser Auftrag nicht enthält.
4. Umgebung wie in `agent-server-b4d.md`: GPU, nvcc, GCC, Python.

Die Tabellen als Kommentar am zugehörigen Issue und als Abschnitt „Messlauf B5a“ in
`docs/testfaelle.md` (unter „K6 in 3D“) auf dem Branch einchecken. Sonst nichts ändern.

## Nicht Teil dieses Auftrags

- Kopplung mit Raytatouille (B5b), Leistungsnormierung gegenüber Raytatouille.
- Klärung der Strehl-Definition gegenüber wende2023 (0,97 gegen 0,997).
- Änderungen am Kern, an Sollwertdateien oder an der Python-Referenz.
