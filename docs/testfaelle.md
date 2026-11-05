# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ offen.

| Nr. | Prüfung | Referenz | Abnahme | Stufe | Status | Test |
| --- | --- | --- | --- | --- | --- | --- |
| K1 | Transmissions- und Reflexionsmatrix der Quelle (Gl. (2), (3), (8), (9)) an einer ebenen Grenzfläche | vektorieller Fresnel-Operator in Convention A (unabhängig aufgebaut) | Abweichung < 1e−13 relativ, auch absorbierend, Totalreflexion, senkrechter Einfall | B0 | ✅ | `test_interface.py` |
| K2 | homogenes Medium | Winkelspektrum | exakt (Rundung) | B1 | ○ | |
| K3 | ebene Grenzfläche, schräger Einfall, FPWPM | Fresnel-Transmission | Konvergenz in dz | B1 | ○ | |
| K4 | Platte schräg beleuchtet, BWPM | Airy-Formel für R und T, TE und TM | R + T = 1, Konvergenz über Iterationen wie wende2023 Abb. 2 | B2 | ○ | |
| K5 | Parabolspiegel (Silber) | Brennweite f′ | Fokuslage, wende2023 Abschnitt 2.1.4 als Plausibilität | B2/B3 | ○ | |
| K6 | Faserfacette mit Index-Matching-Kappe | Fresnel an der Facette; wende2023 Abschnitt 2.2.2 als Plausibilität | Rückreflex, Fokusversatz | B5 | ○ | |

## K1: Grenzflächenmatrizen (B0)

Messwerte und Fälle: `docs/spezifikation.md`, Abschnitt 3.
