# BWPM

Schnelle bidirektionale vektorielle Wellenpropagation (BWPM) nach M. Wende, J. Drozella und
A. M. Herkommer, Opt. Express 31(18), 28874 (2023): Das Querfeld wird schrittweise längs z mit
FFT propagiert; Grenzflächen wirken in der Frequenzdarstellung über Transmissions- und
Reflexionsmatrizen, vor- und rücklaufende Felder werden in Iterationen überlagert. Damit lassen
sich reflektierende und katadioptrische Mikrooptik, Faserfacetten und Rückreflexe in Volumina
rechnen, die für rigorose Löser zu groß sind.

**Ziel:** Referenzimplementierung in Python, dann C++/CUDA-Kern mit Python-Anbindung, und
Kopplung mit Raytatouille und Field Tracing über Gitterfeld-Ports („Raytatouille Fusion“). Deshalb
gelten hier von Anfang an die Konventionen von Raytatouille.

**Stand:** B0–B2 abgenommen (BWPM in 2D), Grenzen in E5, Fusionspilot in 2D
([Bericht](docs/pilot-imc-2d.md)), B3 (3D) abgenommen, B4a–B4d (Sollwerte, C++-Kern, CUDA-Kern) abgenommen, siehe
[docs/roadmap.md](docs/roadmap.md).

## Installieren und testen

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements-dev.txt
pip install -e .            # mit Fusionspilot: pip install -e ".[fusion]"
pytest
ruff check . && ruff format --check . && mypy
```

Python ≥ 3.10; Laufzeitabhängigkeiten nur NumPy und SciPy. `pip install -e .` baut zusätzlich den
C++-Kern `bwpm._core` (C++20-Compiler nötig); ohne Compiler:
`pip install -e . -C cmake.define.BWPM_BUILD_CORE=OFF`.

## Dokumentation

- [Plan und Stufen](docs/roadmap.md)
- [Spezifikation](docs/spezifikation.md)
- [Konventionen](docs/konventionen.md) (verbindlich, wie Raytatouille)
- [Testfälle](docs/testfaelle.md)
- [Fusionspilot in 2D](docs/pilot-imc-2d.md)
- [Entwurf B4: C++/CUDA-Kern, Ports, Adapter](docs/b4-entwurf.md)
- [Quellen](docs/quellen.md)
- [Regeln für Programmieragenten](AGENTS.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
