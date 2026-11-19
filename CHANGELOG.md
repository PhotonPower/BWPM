# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt

- Projektgerüst: Regeln für Agenten, Plan (`docs/roadmap.md`), Konventionen wie Raytatouille,
  Quellen, Testfälle, CI (ruff, mypy, pytest unter Linux und Windows).
- Spezifikation B0 (`docs/spezifikation.md`): Formeln von FPWPM und BWPM aus wende2023 in den
  Konventionen von Raytatouille, offene Punkte.
- `bwpm.interface`: Fresnel-Koeffizienten (Convention A), Transmissions- und Reflexionsmatrix nach
  wende2023 Gl. (2), (3), (8), (9); Nachweis, dass sie mit Convention A exakt der vektorielle
  Fresnel-Operator sind (K1).
- Entscheidungen E2–E4 (Randmaske, Prüfweg für den Rückwärtslauf, Brechzahlstufen) in Plan
  und Spezifikation; Testfälle K7, K8.
- B1: Vorwärtspropagation der FPWPM in 2D (`bwpm.fpwpm`: `Grid2D`, `Propagator2D`,
  `absorbing_mask`, `power_flux`, `longitudinal_component`) und Brechzahlstufen aus einer Geometrie
  (`bwpm.geometry`: Mittelpunkt, Flächenanteil). Testfälle K2, K3 (exakt), K7 (Regel für die
  Randzone), K8 (Gültigkeit an gekippten Flächen: Richtung bis 45°, Leistung nur auf einige Prozent).
- B2: BWPM in 2D (`bwpm.bidirectional.Bidirectional2D`): Reflexe im Vorwärtslauf, Rückwärtslauf,
  Iterationen, Energieresiduum. K4: exakt gegen die Transfer-Matrix des diskretisierten Stapels
  (Tabelle 1 der Quelle bestätigt); K5 Parabolspiegel (Fokus 49,87 µm bei f′ = 50 µm); K9
  gekippter Spiegel (Richtung konvergiert, Leistung nicht – Grenze des Verfahrens).
- Entscheidung E5: Gültigkeitsgrenzen aus B1/B2 (gekippte Flächen, gefangene Wellen).
- Fusionspilot in 2D (`bwpm.pilot_imc`, K6): Faserfacette mit Index-Matching-Kappe, Kette Field
  Tracing (PFT, Winkelspektrum) → BWPM → Field Tracing; Extra `fusion` mit `fieldtracing`.
  Rückreflex 82,4 µm aus der Zwischenbildebene verschoben, an der Facette 4,2 % → 5,6·10⁻⁵.
- B3: BWPM in 3D (`bwpm.volume3d`: `Grid3D`, `Bidirectional3D`, `absorbing_mask_2d`,
  `memory_estimate`) mit Reflexspeicher nur an Grenzflächenpixeln; K10–K12; Messskript
  `experiments/b3_scaling.py` für den Server. Streifende Wellen (k_z = 0) regularisiert.
- Messlauf B3 eingetragen (3D-Strehl 0,997, Laufzeit und Speicher); `Bidirectional3D(workers=…)`
  mit `scipy.fft` für mehrere Kerne; Speichermodell um die Materialkennungen ergänzt.
- B4a: Entwurf `docs/b4-entwurf.md` (Aufbau, Algorithmus im Kern, Ports, Adapter, Teilstufen) und
  Sollwertdateien für den Kern (`tools/export_reference.py`, `reference/*.result.json` im Format
  `raytatouille-result`, `bwpm.results` aus Fieldtracing übernommen).
- B4b: Build mit scikit-build-core, CMake-Presets und pybind11 (E7); C++-Kern `bwpm._core` mit k_z,
  Fresnel-Koeffizienten und Grenzflächenmatrizen; Fassade `bwpm.fast`; CI baut den Kern unter
  Linux und Windows und prüft zusätzlich die Installation ohne Compiler.
- B4c: CPU-Kern `bwpm._core.Solver` (C++20, pocketfft, OpenMP) für FPWPM und BWPM in 2D und 3D
  mit dünnem Reflexspeicher; `bwpm.fast.run_2d`/`run_3d`; Messskript `--core`. Trifft alle
  Sollwertdateien; einkernig 3,7-mal schneller als NumPy. pocketfft (BSD-3) unter
  `cpp/third_party/pocketfft`.
- Messlauf ohne Kappe eingetragen (3D-Rückreflex der blanken Facette 3,93 %; Threads in NumPy
  bringen kaum Gewinn).
- B4d (erste Fassung, Test auf dem Server offen): CUDA-Kern `bwpm._core.CudaSolver` (cuFFT,
  T und R mit Ausbreitungsfaktor im Kernel, dünner Reflexspeicher auf dem Gerät), CMake-Option
  `BWPM_ENABLE_CUDA`, Fassade `backend="auto"|"cpu"|"cuda"`, Tests K15, Messskript `--backend`.
- Regel für streifenden Einfall (#25): |ñ²k₀² − k⊥²| ≤ 10⁻¹²·|ñ²k₀²| gilt in Python, C++ und CUDA
  als exakt streifend; vorher hing das Ergebnis an diesen Punkten vom Compiler ab.
- B4d gemessen (K15): Faserfacette bei λ/(5 n) in 240 s auf der RTX 6000 Ada (130–150-mal
  schneller als NumPy); `BWPM_STATIC_LIBSTDCXX` für conda-Pythons; Messskript zählt die
  Brechzahlpaare ebenenweise (vorher > 40 min nach dem Lauf).
- B4d auf dem Server bestätigt (alle Tests grün mit CUDA und gcc-toolset-15).
- B4e: Fusion-Adapter `bwpm.adapter` (`describe`, `prepare` mit Inhalts-Hash, `apply`, `estimate`,
  Ports mit `port_power` in 3D, `plane_wave`) und Konformitätstests K16 für alle Backends:
  Transfer-Matrix, Airy, Energieerhaltung, Reziprozität, Determinismus, zentrierter Port.
- B5a: Auswertung des 3D-Piloten (`bwpm.pilot_imc3d`: Fluss, flussgewichtete Fresnel-Referenz, Fokus auf der Achse) mit Tests; `experiments/b3_scaling.py` meldet Reflexionsgrad gegen Referenz, Fokus des Reflexes und rücklaufenden Fluss an der Facette (K6 in 3D, Messlauf auf dem Server offen, `docs/agent-server-b5a.md`).
- B5a: Bericht `docs/pilot-imc-3d.md` zum Fusionspilot in 3D (Messlauf und Diagnose auf dem Server); K6 in 3D vom Maintainer abgenommen.
- B6 (#29): erste Born-Näherung als Sollwert für Streuung an schwachen Schichten
  (`bwpm.born.born_forward_2d`), geprüft gegen Airy und eine RCWA des Lamellengitters
  (`tests/rcwa_reference.py`); K17 vergleicht die FPWPM je Winkel und Polarisation, Messskript
  `experiments/born_angles.py`. Befund: paraxial < 10⁻³; bei großen Winkeln TE-Amplitude zu klein
  (0,80 bei 46°), für TM fehlt der Faktor cos²θ; konvergiert in dy = dz (Grenze des Verfahrens).
- B6: Ursache der Abweichung analytisch in erster Ordnung (Schrägfaktor der einfallenden statt
  der gestreuten Welle, Fresnel an Ein- und Austritt): TE cos θ (3 − cos θ)/2, TM
  (3 − cos θ)/(2 cos θ); Test gegen die FPWPM auf 10⁻³.
- B6: vektorielle Born-Referenz in 3D (`bwpm.born.born_forward_3d`, Querblock der dyadischen
  Green-Funktion), geprüft gegen 2D; K17 in 3D: FPWPM komponentenweise skalar mit dem Faktor r(θ),
  an einer dünnen Zufallsschicht bestätigt. Hinweise für Scattuccino S9b-2 (skalare Referenz,
  fehlendes E_z) in `docs/hinweise-scattuccino-s9b2.md`, Skript `experiments/vector_born_wm.py`.
- Entscheidung E9: Gültigkeitsgrenze bei Streuung in große Winkel (quantitativ bis ≈ 11° für s,
  6,6° für p); Vorschlag für schwache Fluktuationen in #37. B6 abgeschlossen.
- K17, Punkt 8 (Hinweis aus Scattuccino S9b-4): Die Extinktion des kohärenten Feldes in der
  FPWPM folgt der gestreuten Quernorm (∫ σ r²/cos θ dΩ), nicht dem vektoriellen Born; Test je
  Realisierung, μs-Spalte in `experiments/vector_born_wm.py`, Nachtrag in den Hinweisen.
- Nachtrag zu Scattuccino S9b-7 (PR #33) in `docs/hinweise-scattuccino-s9b2.md`: K17-Vorhersage
  für 1 − g des vollen Feldes (neue Spalte in `experiments/vector_born_wm.py`); E9 erklärt nur
  etwa 0,01 der gemessenen g-Abweichung von 0,042.
- Entscheidung E10: Variante B3 (#37) angenommen, neue Stufe B7 „Schwache Fluktuationen“
  (B7a Python-Referenz, B7b Kern) mit Abnahme K18.
- B7a, 2D: Born-Inkrement für schwache Fluktuationen δn (`Propagator2D(fluctuation=...)`,
  Spezifikation §10, Quelle chen2020). K18: K17-Sollwerte bei 29–63° auf 0,05–0,23 % getroffen
  (reine FPWPM 7–180 %), Mehrfachstreuung gegen RCWA auf < δn, Fluss erster Ordnung in dz,
  δn = 0 bitgleich.
- B7a, 3D: Born-Inkrement in `Bidirectional3D(fluctuation=...)` (gemeinsamer Kern
  `bwpm.fluctuation`); K18 in 3D gegen vektorielles Born an einer Zufallsschicht auf 0,05–0,9 %
  bis 75° (reine FPWPM 4–57 %), zweiter Ordnung im Schritt.
- Entscheidung E11 (#46): linearisiertes Potenzial V = 2ñ₀k₀²δn im Born-Inkrement statt
  chen2020 Gl. (1); Fehler gegen RCWA 0,1–0,8 % statt 0,3–2,1 %, neuer Test in K18. B7a
  abgenommen.
- B7b, CPU-Kern: Born-Inkrement in `bwpm._core.Solver` (Argument `fluctuation`) und in
  `bwpm.fast.run_2d/run_3d`; Sollwertdateien `fpwpm2d_fluctuation`, `bwpm3d_fluctuation`. CUDA noch
  ohne δn (`backend="auto"` nutzt dann die CPU).
- Korrektur am Born-Inkrement (gefunden in B7b): Die Quelle nimmt nur propagierende, nicht
  streifende Wellen (Regel aus §3, #25); an exakt streifenden Gitterpunkten schaukelte sich ein
  zufälliges Feld sonst auf (K18 Punkt 6).
- B7b, CUDA-Kern: Born-Inkrement in `bwpm._core.CudaSolver` (cuFFT für drei Komponenten, δn je
  Scheibe hochgeladen); `backend="auto"` nutzt CUDA jetzt auch mit δn. Geprüft auf einer RTX 2080
  (K18 Punkt 8, mit Laufzeiten). B7b mit CUDA 13.4 und 13.5 abgenommen; B7 abgeschlossen.
- Build: Unter Windows mit CUDA 13 übergibt CMake `/Zc:preprocessor` an MSVC (CCCL verlangt den
  standardkonformen Präprozessor).
- Windows mit CUDA: Der Build hinterlegt den CUDA-Ordner in `bwpm/_cuda_dll.py`; `bwpm.fast`
  fügt ihn vor dem Import von `bwpm._core` zum DLL-Suchpfad hinzu (ersatzweise `CUDA_PATH`). Bisher
  fand Python `cufft64_12.dll` nicht, weil es unter Windows kein RPATH gibt.
- Entscheidung E12: Rückstreuung an δn als Stufe B8; Entwurf `docs/b8-entwurf.md` (Quelle:
  Zusatzmaterial zu chen2020, im Volltext geprüft).
- B8a/B8b: Rückstreuung an δn in der BWPM (R1–R3, Spezifikation §11) in Python
  (`Bidirectional2D`, `Bidirectional3D`), CPU- und CUDA-Kern; Scheibenintegral mit zwei
  Gauß-Punkten statt Mittelpunkt (F1, #53), auch vorwärts. K19: Schicht und Bragg-Gitter gegen die
  exakte Lösung auf 5·10⁻⁴ bzw. 1,3·10⁻³. Behoben: `Bidirectional2D` ignorierte δn (Lücke aus B7).
  Sollwertdateien für δn neu erzeugt.
