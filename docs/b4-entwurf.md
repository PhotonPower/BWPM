# Entwurf B4: C++/CUDA-Kern, Ports und Adapter

> Entwurf zur Freigabe (2026-11-05). Grundlage: Python-Referenz aus B0–B3, Fusionsplan
> (Fieldtracing `docs/rtt-fusion/plan.md`, Abschnitte 6.4–6.6), Pilot `docs/pilot-imc-2d.md`,
> Messlauf B3 (`docs/testfaelle.md`). Aufbau angelehnt an FourierGrating (gleicher Maintainer):
> C++20, CMake-Presets, scikit-build-core, pybind11, CUDA für sm_89.

## 1. Ziel

Die BWPM soll Volumina in der Auflösung der Quelle (λ/(5 n), etwa 1100² × 1000 in 3D) in Minuten
statt Stunden rechnen und als Löser-Modul der Fusion über Ports ansprechbar sein. Die
Python-Referenz bleibt unverändert die Prüfgrundlage; der Kern muss die Sollwertdateien unter
`reference/` (Format `raytatouille-result`) in der dort angegebenen Toleranz treffen.

Messlauf B3 als Ausgangspunkt: 800² × 708, zwei Iterationen, NumPy auf einem Kern: 2 h 26 min,
etwa 3 s je Ebene und Durchlauf. Ziel des Kerns auf der RTX 6000 Ada: < 0,1 s je Ebene und
Durchlauf bei 1152² (Faktor > 30 gegenüber NumPy).

## 2. Aufbau (Vorschläge zur Entscheidung)

| Nr. | Frage | Vorschlag | Begründung |
| --- | --- | --- | --- |
| D1 | Paketaufbau | ein Paket `bwpm`; kompiliertes Modul `bwpm._core` (optional), die Python-Referenz bleibt `bwpm.*`; neue Fassade `bwpm.fast` wählt `_core`, sonst die Referenz | Referenz und Kern im selben Repository, Sollwerte nebeneinander prüfbar; Installation ohne Compiler bleibt möglich |
| D2 | Build | CMake ≥ 3.24 mit Presets `cpu-release`, `cpu-debug` (Sanitizer), `cuda-release` (sm_89); Python über scikit-build-core | wie FourierGrating |
| D3 | Bindings | pybind11 | wie FourierGrating, hp-FEM, CliffordBem |
| D4 | FFT | CPU: pocketfft (C++-Header, BSD-3; dieselbe Bibliothek wie NumPy/SciPy, daher Abweichung zur Referenz im Bereich der Rundung); GPU: cuFFT | MIT-verträglich (kein FFTW/GPL), keine Laufzeitabhängigkeit auf der CPU |
| D5 | Genauigkeit | complex128 überall; complex64 auf der GPU erst als Option mit eigener Toleranz nach Messung | Sollwerte sind double; K4 zeigt Rundungsgenauigkeit |
| D6 | Operatoren | auf der GPU T und R **im Kernel berechnen** statt zwischenspeichern (geschlossene Formeln, wenige Dutzend Rechenoperationen je Punkt); auf der CPU Zwischenspeicher wie in Python | Speicher: der Cache kostet Arten × Paare × N × 64 B (bei 1152² und 6 Paaren 2 GB), die Rechnung im Kernel ist billiger als das Lesen |
| D7 | Parallelität CPU | OpenMP über die Brechzahlpaare einer Ebene und in den FFTs (pocketfft mit Threads) | Messlauf: die FFTs je Paar dominieren |
| D8 | Zielplattform GPU | CUDA ≥ 12.6, sm_89 (RTX 6000 Ada); weitere Architekturen über `CMAKE_CUDA_ARCHITECTURES` | Rechner des Maintainers; wie FourierGrating |

## 3. Algorithmus im Kern

Je Ebene und Durchlauf (vorwärts bzw. rückwärts):

1. eine 2D-FFT des Querfelds (zwei Komponenten, als Batch);
2. für jedes vorkommende Brechzahlpaar: Multiplikation mit der 2 × 2-Matrix (T oder R mit
   Ausbreitungsfaktor nach Tabelle 1), auf der GPU mit T, R im Kernel; inverse FFT; Übernahme der
   Pixel des Paars (Maske I_r,s) – alle Paare einer Ebene als ein cuFFT-Batch;
3. Randmaske in denselben Kernel;
4. Reflexquellen nur an Grenzflächenpixeln (Index und Wert, wie `bwpm.volume3d`).

Brechzahlen als Tabelle, Kennungen als uint8 (≤ 255 Materialien) oder uint16. Zwischenebenen nur
auf Anforderung. Gefangene Wellen: Iterationszahl begrenzen und am Energieresiduum abbrechen
(Grenze aus K4/K11); kein Filter im ersten Schritt.

## 4. Ports und Adapter (Fusion)

Port nach dem Piloten (`docs/pilot-imc-2d.md`, Abschnitt „Was der Pilot festlegt“):

- Ebene z = const, Richtung (+z oder −z), Medium am Port, Vakuumwellenlänge;
- Querfeld (E_x, E_y) auf einem zentrierten Gitter x_m = (m − N_x//2)·dx, y_n = (n − N_y//2)·dy;
  Spektren auf (0, 0) bezogen;
- Leistung als Poynting-Fluss (`power_flux`); Umrechnung auf leistungsnormierte Strahlen von
  Raytatouille an der Port-Grenze (Fieldtracing ADR 0004).

Adapter `bwpm.adapter` nach Fusionsplan 6.4:

| Aufruf | Inhalt |
| --- | --- |
| `describe()` | Verfahren, Version, Konventionen (e^(−iωt), Fresnel Convention A, s = TE), Ports (Gitterfeld ±z), linear ja, reziprok ja, deterministisch (CPU bitgleich, GPU in Toleranz), Ressourcen, Gültigkeitsgrenzen aus E5 |
| `prepare(structure, wavelengths)` | Brechzahlvolumen (Kennungen, Tabelle, Gitter), Inhalts-Hash |
| `apply(handle, port_in, field)` | Felder an den Ports z₀ (rücklaufend) und z₁ (vorlaufend), optional Zwischenebenen |
| `estimate(handle, request)` | Laufzeit und Speicher aus dem Modell des Messlaufs (`memory_estimate`, Zeit je Gitterpunkt) |
| `operator(...)` | nicht angeboten: für große Volumina zu teuer; lineare Wiederverwendung über den Cache der Fusion |

## 5. Konformitätsprüfung (Fusionsplan 6.5)

Für den Kern und den Adapter, gegen unabhängige Werte:

- ebene Welle durch eine Platte: Airy und Transfer-Matrix (K4, K11), TE und TM, schräg, Silber;
- Energieerhaltung bei verlustfreien ebenen Strukturen; Reziprozität (Vertauschen von Ein- und
  Ausgang an einem symmetrischen Stapel);
- Konventionen: Zeitfaktor, s/p-Basis und Phasenbezug an den Ports (Vergleich mit FourierGrating
  für einen Schichtstapel als Querprobe, sobald dessen Adapter existiert);
- Determinismus: CPU bitgleich zwischen Läufen, GPU in der erklärten Toleranz.

## 6. Teilstufen

| Teilstufe | Inhalt | Abnahme |
| --- | --- | --- |
| B4a | Sollwertdateien aus der Python-Referenz (`tools/export_reference.py`, `reference/`) | Dateien aktuell, Rundreise, unabhängige Werte (✅) |
| B4b | Build-Gerüst (CMake, Presets, scikit-build-core, pybind11, CI mit Compiler), C++ `interface` (T, R) | `interface_matrices` in Toleranz (✅, K13) |
| B4c | CPU-Kern 2D/3D: Schritt, Reflexquellen, Iterationen, Maske (pocketfft, OpenMP) | alle Sollwertdateien; Laufzeit gegen NumPy (✅, K14) |
| B4d | CUDA-Kern (cuFFT, T/R im Kernel, Batches) | Sollwertdateien in GPU-Toleranz; Messlauf IMC bei λ/(5 n) auf dem Server (gemessen: 240 s, K15) |
| B4e | Adapter `bwpm.adapter` und Konformitätstests | Abschnitt 5 |

## 7. Risiken

- **Rundung CPU/GPU:** cuFFT und pocketfft unterscheiden sich in der Rundung; Toleranzen je
  Sollwertdatei relativ zum Feldmaximum.
- **Batch-Speicher auf der GPU:** viele Paare in Faserebenen (Kern, Mantel, Kappe) vergrößern den
  Batch; Obergrenze und Aufteilung nach Speicher.
- **Windows-Build:** MSVC und CUDA; FourierGrating zeigt, dass es geht (MinGW-Sonderweg dort).

## 8. CUDA bauen und prüfen (B4d)

Im Container gibt es keinen CUDA-Compiler; der Gerätecode wurde dort nur mit clang 18 gegen die
Header von CUDA 12.9 übersetzt (Host- und Gerätecode ohne Fehler und Warnungen), nicht gelinkt
und nicht ausgeführt. Auf dem Rechner des Maintainers:

```bash
pip install -e . -C cmake.define.BWPM_ENABLE_CUDA=ON      # nvcc ≥ 12.6 im PATH oder CUDACXX
python -c "from bwpm import _core; print(_core.HAS_CUDA, _core.cuda_device_count())"
BWPM_REQUIRE_CUDA=1 pytest tests/test_core_cuda.py -v
python experiments/b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2 --core --backend cuda
```

Auf dem Rechner des Maintainers funktionierte (Bericht in PR #23): nvcc 12.9, gcc-toolset-15 als
Host-Compiler, `-C cmake.define.BWPM_STATIC_LIBSTDCXX=ON` (conda-Python mit älterer libstdc++).
GCC 14 (AlmaLinux 10) wird von nvcc erst ab CUDA 12.8 unterstützt; mit älterem CUDA hilft
`-C cmake.define.CMAKE_CUDA_FLAGS=-allow-unsupported-compiler` oder ein älterer Host-Compiler
(`-C cmake.define.CMAKE_CUDA_HOST_COMPILER=...`).

Erste Fassung: Paare einer Ebene nacheinander (je ein Kernel für Matrix × Spektrum mit T/R im
Kernel, eine inverse cuFFT, ein Kernel für die Übernahme der Pixel); Kennungen je Ebene
hochgeladen. Der Batch aller Paare einer Ebene folgt nach der ersten Messung.
