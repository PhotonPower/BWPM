# Auftrag für einen Agenten: CUDA-Kern (B4d) auf dem Server bauen, testen, messen

> Für einen Programmieragenten auf dem Rechner des Maintainers (AlmaLinux 10, RTX 6000 Ada,
> 1 TB RAM). Es gelten `AGENTS.md` und `docs/konventionen.md`. Zusammenhang: `docs/b4-entwurf.md`
> (Abschnitt 8), PR #23, Issue #22.

## Ziel

1. `bwpm._core` mit CUDA bauen, sodass `from bwpm import _core` ohne Umgebungsvariablen lädt.
2. Die Tests K15 (`tests/test_core_cuda.py`) grün bekommen – **ohne Toleranzen, Tests oder
   Sollwertdateien zu ändern**.
3. Laufzeit und Speicher messen (CPU-Kern gegen CUDA-Kern) und die Ergebnisse berichten.

## Regeln

- Arbeite auf dem Branch `b4d-cuda`. Nicht mergen, nicht auf `main` pushen.
- Erlaubt: Änderungen an `CMakeLists.txt`, `CMakePresets.json`, `pyproject.toml` (Build) und
  Fehlerbehebungen in `cpp/src/solver_cuda.cu` / `cpp/include/bwpm/solver_cuda.hpp`, jeweils als
  eigener Commit mit Begründung in der Commit-Message (englisch, Imperativ).
- Verboten: Tests, Toleranzen, `reference/*.result.json`, Python-Referenz (`src/bwpm/*.py`
  außer `fast.py` bei echten Fehlern) ändern, um Tests grün zu machen. Weicht der CUDA-Kern ab,
  ist das ein Befund, kein Grund für eine größere Toleranz.
- Keine Pfade, Hostnamen oder Tokens in Dateien oder Commits.
- Wenn etwas unklar ist oder eine Entscheidung braucht: anhalten und berichten.

## Schritt 1: Umgebung erfassen

```bash
cd ~/BWPM && git fetch && git checkout b4d-cuda && git pull
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
which nvcc && nvcc --version | tail -2
gcc --version | head -1
python --version && which python
find / -name "libcufft.so*" 2>/dev/null | head
```

Festhalten: GPU, Treiber, CUDA-Version von nvcc, GCC-Version, Python (miniconda), Fundorte von
cuFFT. Kommt nvcc aus conda, liegen die Bibliotheken oft unter `$CONDA_PREFIX/lib` oder
`.../targets/x86_64-linux/lib`.

## Schritt 2: Bauen

```bash
rm -rf build
pip install -e . -C cmake.define.BWPM_ENABLE_CUDA=ON -v 2>&1 | tee /tmp/bwpm_build.log | grep -E "BWPM|error" 
python -c "from bwpm import _core, fast; print(_core.HAS_CUDA, _core.cuda_device_count(), fast.cuda_available())"
```

Erwartet: `True 1 True`. Die Zeile `BWPM CUDA: cuFFT …; RPATH …` im Log nennt die gelinkte
cuFFT-Bibliothek und den eingetragenen Suchpfad.

**Fehlerbilder und Abhilfe** (der Reihe nach probieren, jedes Ergebnis notieren):

| Fehler | Ursache | Abhilfe |
| --- | --- | --- |
| `nvcc` meldet nicht unterstützten Host-Compiler | GCC 14 mit CUDA < 12.8 | `-C cmake.define.CMAKE_CUDA_FLAGS=-allow-unsupported-compiler` oder `-C cmake.define.CMAKE_CUDA_HOST_COMPILER=<älterer g++>` |
| `find_package(CUDAToolkit 12.6)` scheitert | nvcc nicht gefunden oder zu alt | `CUDACXX=<pfad>/nvcc` setzen; CUDA ≥ 12.6 nötig |
| `ImportError: libcufft.so.11` | RPATH zeigt nicht auf cuFFT | 1. `readelf -d <…>/_core*.so \| grep -i path` mit dem Fundort aus Schritt 1 vergleichen; 2. CMake korrigieren (RPATH aus `IMPORTED_LOCATION` von `CUDA::cufft`, siehe `CMakeLists.txt`); 3. als Ausweg `-C cmake.define.BWPM_CUFFT_STATIC=ON` |
| `ImportError` mit anderer CUDA-Bibliothek | weitere dynamische Abhängigkeit | `ldd <…>/_core*.so \| grep "not found"`; Pfad wie oben ergänzen |
| Fehler in `solver_cuda.cu` beim Übersetzen | nvcc strenger als clang (Code wurde nur mit clang 18 geprüft) | Fehler beheben, eigener Commit |

`<…>/_core*.so` findet man mit
`python -c "import bwpm, glob, os; print(glob.glob(os.path.join(os.path.dirname(bwpm.__file__), '_core*.so')))"`.

Bei `BWPM_CUFFT_STATIC=ON`: Größe von `_core*.so` notieren.

## Schritt 3: Tests

```bash
pytest -q                                   # alles: CPU-Kern, Referenz, Python
BWPM_REQUIRE_CUDA=1 pytest tests/test_core_cuda.py -v
```

Erwartet: alle grün. K15 vergleicht den CUDA-Kern mit allen fünf Sollwertdateien (relativ
1e−10, absolut 1e−11 des Maximums) und mit dem CPU-Kern auf zufälligen 3D-Volumina.

**Schlägt K15 fehl:** nichts an den Tests ändern. Stattdessen die größten Abweichungen je Fall
messen und berichten:

```bash
python - <<'PY'
import numpy as np
from pathlib import Path
from bwpm import fast
from bwpm.results import load_json
from bwpm.fpwpm import Grid2D, Propagator2D
from bwpm.volume3d import Bidirectional3D, Grid3D
root = Path(".")
for name, it in [("fpwpm2d_interface", 0), ("bwpm2d_tilted_silver", 2), ("bwpm2d_slab_exact", 40)]:
    d = load_json((root / "reference" / f"{name}.result.json").read_text("utf-8"))[1]
    p = d["parameters"]
    prop = Propagator2D(Grid2D(p["ny"], p["dy_mm"], p["dz_mm"], p["wavelength_um"]),
                        list(d["inputs"]["indices"]), d["inputs"]["labels"],
                        d["inputs"]["index_before"], d["inputs"].get("mask"))
    for backend in ("cpu", "cuda"):
        r = fast.run_2d(prop, d["inputs"]["field"], it, backend=backend)
        ref = d["expected"].get("transmitted")
        if ref is None:
            ref = d["expected"]["planes"][-1] if "planes" in d["expected"] else d["expected"]["forward"][-1]
        print(name, backend, float(np.max(np.abs(r.transmitted - ref)) / np.max(np.abs(ref))))
PY
```

Hinweise zur Fehlersuche, falls nur CUDA abweicht: Abweichung ~1e−16 … 1e−13 ist Rundung
(cuFFT gegen pocketfft); Abweichungen ~1e−3 oder größer deuten auf einen Fehler in Layout,
Vorzeichen oder Normierung (Stride/Batch des cuFFT-Plans, Zuordnung r/s im Rückwärtslauf,
Normierung 1/N der inversen FFT, Frequenzen k_x/k_y). Den CPU-Kern `cpp/src/solver.cpp` als
Vorbild nehmen; er trifft alle Sollwerte.

## Schritt 4: Messung

Erst die CUDA-Läufe, dann ein CPU-Vergleich bei der kleinsten Auflösung:

```bash
nvidia-smi --query-gpu=memory.used --format=csv -l 5 > /tmp/gpu_mem.log &
python experiments/b3_scaling.py --fractions 2.5 3.5 5 --window-um 96 --iterations 2 --core --backend cuda
kill %1
python experiments/b3_scaling.py --fractions 2.5 --window-um 96 --iterations 2 --core --backend cpu
python experiments/b3_scaling.py --fractions 2.5 --window-um 96 --iterations 2 --bare --core --backend cuda
```

Ergebnisse: `experiments/out/*.json` (nicht einchecken). Maximum aus `/tmp/gpu_mem.log` notieren.

Zum Vergleich (NumPy, ein Kern, 2 Iterationen, mit Kappe): λ/(2,5 n) 572² × 506 in 51 min,
λ/(3,5 n) 800² × 708 in 2 h 26 min; Strehl 0,99694 / 0,99695; Rückreflex 4,185 % / 4,180 %.
Ohne Kappe (λ/(2,5 n)): 19,6 min, Rückreflex 3,93 %. Die CUDA-Läufe müssen dieselben Strehl-
und Rückreflexwerte liefern (Abweichung erst ab der 6. Stelle zulässig).

## Schritt 5: Bericht

Als Kommentar in PR #23 (oder als Datei für den Maintainer, falls kein Schreibzugriff):

1. Umgebung (Schritt 1) und welche Build-Variante funktioniert hat (RPATH oder statisch), mit
   den nötigen `-C`-Optionen.
2. Commits auf `b4d-cuda` mit Kurzbegründung.
3. Testergebnis K15 (bestanden/nicht bestanden, bei Fehlschlag die gemessenen Abweichungen).
4. Tabelle der Messung:

| Lauf | Gitter | Backend | Laufzeit | Host-Speicher | GPU-Speicher (max) | Strehl | Rückreflex |
| --- | --- | --- | --- | --- | --- | --- | --- |

5. Beobachtungen (z. B. Zeit je Ebene, Auffälligkeiten).
