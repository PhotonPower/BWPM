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

**Stand:** B0 und B1 abgenommen, B2 (BWPM in 2D) umgesetzt, siehe
[docs/roadmap.md](docs/roadmap.md).

## Installieren und testen

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements-dev.txt
pip install -e .
pytest
ruff check . && ruff format --check . && mypy
```

Python ≥ 3.10; Laufzeitabhängigkeiten nur NumPy und SciPy.

## Dokumentation

- [Plan und Stufen](docs/roadmap.md)
- [Spezifikation](docs/spezifikation.md)
- [Konventionen](docs/konventionen.md) (verbindlich, wie Raytatouille)
- [Testfälle](docs/testfaelle.md)
- [Quellen](docs/quellen.md)
- [Regeln für Programmieragenten](AGENTS.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
