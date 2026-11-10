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
- B5a: Bericht `docs/pilot-imc-3d.md` zum Fusionspilot in 3D (Messlauf und Diagnose auf dem Server, Abnahme offen).
