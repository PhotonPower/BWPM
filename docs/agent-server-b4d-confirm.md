# Auftrag für einen Agenten: Bestätigungslauf B4d auf dem Server

> Für einen Programmieragenten auf dem Rechner des Maintainers (RTX 6000 Ada). Es gelten
> `AGENTS.md`, `docs/konventionen.md` und die Regeln aus `docs/agent-server-b4d.md`. Kontext:
> PR #23, Issue #25 (Regel für streifenden Einfall, Commit `805aaa1`).

## Ziel

Bestätigen, dass nach der Regel für streifenden Einfall (#25) **alle** Tests auf dem Server grün
sind – mit dem C++-Kern unter GCC und mit CUDA. Nichts messen, nichts optimieren.

## Regeln

- Branch `b4d-cuda`; nicht mergen, nicht auf `main` pushen.
- Keine Tests, Toleranzen, Sollwertdateien oder die Python-Referenz ändern.
- Erlaubt sind nur Build-Korrekturen (`CMakeLists.txt`, `pyproject.toml`), jeweils als eigener
  Commit mit Begründung. Alles andere: anhalten und berichten.
- Keine Pfade, Hostnamen oder Tokens in Dateien oder Commits.

## Schritte

1. Stand holen. Der Branch wurde neu geschrieben (Rebase auf `main`), daher hart zurücksetzen:

   ```bash
   cd ~/BWPM && git fetch && git checkout b4d-cuda && git reset --hard origin/b4d-cuda
   git log --oneline -1          # erwartet: 805aaa1 oder neuer
   rm -rf build
   ```

2. Bauen wie im letzten Lauf (gcc-toolset-15, statische libstdc++ wegen conda):

   ```bash
   CXX=<g++-15> CC=<gcc-15> pip install -e . -C cmake.define.BWPM_ENABLE_CUDA=ON \
     -C cmake.define.BWPM_STATIC_LIBSTDCXX=ON -C cmake.define.CMAKE_CUDA_HOST_COMPILER=/usr/bin/g++
   python -c "from bwpm import _core, fast; print(_core.HAS_CUDA, _core.cuda_device_count(), fast.cuda_available())"
   ```

   Erwartet: `True 1 True`.

3. Alle Tests, dann die CUDA-Tests ausdrücklich:

   ```bash
   pytest -q
   BWPM_REQUIRE_CUDA=1 pytest tests/test_core_cuda.py tests/test_core_interface.py tests/test_core_solver.py -v
   ```

   Erwartet: alles grün, insbesondere die vier Tests aus #25
   (`test_equals_python_reference_on_a_grid[1.0-1.5]`, `[1.5-1.0]`, `test_bwpm3d_tilted`,
   `test_random_volumes_equal_python[2]`), der neue `test_grazing_rule_is_independent_of_rounding`
   und `test_core_cuda.py::test_bwpm3d_tilted`.

4. Zur Kontrolle, dass die Regel unter GCC 14 greift, zusätzlich einmal **ohne** CUDA mit dem
   System-Compiler bauen und den CPU-Kern testen:

   ```bash
   rm -rf build
   pip install -e . -C cmake.define.BWPM_STATIC_LIBSTDCXX=OFF   # System-GCC 14
   pytest -q tests/test_core_interface.py tests/test_core_solver.py
   ```

   Scheitert der Import an `GLIBCXX`, diesen Schritt mit gcc-toolset-15 und
   `BWPM_STATIC_LIBSTDCXX=ON` wiederholen und das notieren.

5. Danach wieder den CUDA-Build aus Schritt 2 herstellen, damit der Rechner im erwarteten Zustand
   bleibt.

## Bei Fehlschlägen

- Test, Fehlermeldung und die gemessene größte Abweichung notieren (`pytest -v` zeigt sie im
  `assert`).
- Abweichungen ~10⁻¹⁶ … 10⁻¹³ sind Rundung; ~10⁻⁹ an wenigen Punkten deutet auf weitere
  Streifpunkte hin (dann k⊥ und |ñ²k₀² − k⊥²|/|ñ²k₀²| an den betroffenen Gitterpunkten
  ausgeben); größere Abweichungen sind ein Fehler.
- Nicht reparieren außer Build; berichten.

## Bericht

Als Kommentar in PR #23:

1. Commit-Stand und Build-Befehl, der funktioniert hat.
2. Ergebnis von Schritt 3 und 4: Anzahl bestanden / fehlgeschlagen / übersprungen, je Lauf.
3. Bei Fehlschlägen die Angaben aus dem vorigen Abschnitt.
4. Ein Satz Fazit: „B4d auf dem Server bestätigt“ oder was fehlt.
