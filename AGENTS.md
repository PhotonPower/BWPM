# Regeln für Programmieragenten

Diese Datei liest jeder Agent zu Beginn jeder Sitzung. Sie gilt zusammen mit
`docs/roadmap.md`, `docs/konventionen.md`, `docs/spezifikation.md` und `docs/quellen.md`.

## Projekt in einem Satz

BWPM ist die Referenzimplementierung der schnellen bidirektionalen vektoriellen
Wellenpropagation nach Wende et al. (2023) – erst in Python, dann als C++/CUDA-Kern – und das
erste Modul, das von Anfang an für die Kopplung mit Raytatouille und Field Tracing gebaut wird
(„Raytatouille Fusion“, Plan in Fieldtracing `docs/rtt-fusion/plan.md`).

## Bauen und testen

```bash
pip install -r requirements-dev.txt && pip install -e .
pytest
ruff check . && ruff format --check . && mypy
```

Vor jedem PR müssen lokal grün sein: pytest, ruff (check und format), mypy (strict). Die CI
testet zusätzlich Python 3.10 und Windows. Teure Referenzrechnungen (> 5 s) bekommen den Marker
`@pytest.mark.slow`; `pytest -m slow` führt sie aus. Große Läufe gehören auf den Rechner des
Maintainers; im Repository stehen nur Skripte und Ergebnisse, keine Pfade.

## Arbeitsregeln

1. **Ein Auftrag = ein Ziel einer Stufe** aus `docs/roadmap.md`. Nichts darüber hinaus bauen.
2. **Erst Test, dann Code.** Sollwerte kommen aus einer analytischen Lösung oder einer
   unabhängigen exakten Rechnung (Airy-Formel, Fresnel, Winkelspektrum, Quadratur), nie aus dem zu
   prüfenden Code. Bestehende Tests und Sollwerte nur mit Freigabe des Maintainers ändern; einen
   Test anzupassen, damit er grün wird, ist verboten.
3. **Näherungen prüft man über ihre Ordnung:** Konvergenz in der Schrittweite dz und dx, in der
   Zahl der Iterationen, in der Feinheit der Brechzahlstufen.
4. **Konventionen sind Gesetz** (`docs/konventionen.md`), identisch mit Raytatouille und
   Fieldtracing: mm, µm, rad, k in 1/mm, Felder ∝ exp(i(k·r − ωt)), Fresnel Convention A,
   s = TE, p = TM. Jede öffentliche Funktion nennt Einheiten und Vorzeichen im Docstring.
5. **Formeln mit Quelle:** jede Formel mit Verweis auf die Spezifikation (Abschnitt) und auf die
   Gleichung der Quelle aus `docs/quellen.md`. Gleichungsnummern fremder Quellen nur aus dem
   Volltext.
6. **Determinismus:** Zufall nur mit explizitem Seed (`np.random.default_rng(seed)`).
7. **Fremde Repositories nur lesen** (Raytatouille, Fieldtracing, FourierGrating, hp-FEM,
   CliffordBem, SpeckleBem, BioOptoMD). Vorschläge dafür entstehen als Dokumente unter `docs/`.
8. **Keine Geheimnisse im Repository:** keine Tokens, Passwörter oder persönlichen Pfade in
   Dateien, Commits oder Remote-URLs.
9. **Unklarheit = Frage, nicht Annahme.** Offene Punkte als Issue mit Label `frage`; Abweichungen
   von einer Quelle als Issue mit Label `paper`.

## Sprache und Stil

- Code, Kommentare, Docstrings, Commit-Messages englisch; Projektdokumentation deutsch.
- ruff (Zeilenlänge 100), mypy strict. Funktionen und Variablen `snake_case`, Klassen
  `PascalCase`, Konstanten `UPPER_CASE`.
- Commits im Imperativ, ein Thema pro Commit, z. B. `interface: add reflection matrix`.

## Definition of Done (PR-Checkliste)

- [ ] CI grün (Lint, Typen, Tests unter Linux und Windows)
- [ ] Referenztests mit unabhängigem Sollwert, bei Näherungen mit Konvergenztest
- [ ] Testfall in `docs/testfaelle.md` eingetragen oder aktualisiert (Status, gemessene Werte)
- [ ] Neue Quellen in `docs/quellen.md`
- [ ] Docstrings mit Einheiten, Vorzeichen und Verweis auf Spezifikation bzw. Quelle
- [ ] Eintrag in `CHANGELOG.md` unter „Unreleased“
