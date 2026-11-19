# Quellen

Status: **Volltext geprüft** heißt, dass Formeln und Gleichungsnummern am Volltext verglichen
wurden; **bibliografisch** heißt, nur Autoren, Zeitschrift, Band und Seiten sind belegt.
Gleichungsnummern im Code nur aus Quellen mit Volltextprüfung (AGENTS.md, Regel 5).

| Kürzel | Quelle | Status | Verwendet für |
| --- | --- | --- | --- |
| wende2023 | M. Wende, J. Drozella, A. M. Herkommer, „Fast bidirectional vector wave propagation method showcased on targeted noise reduction in imaging fiber bundles using 3D-printed micro optics“, Opt. Express 31(18), 28874–28890 (2023), doi:10.1364/OE.497244, CC BY 4.0 | Volltext vom Maintainer bereitgestellt (2026-11-05) und gelesen: Gl. (1)–(11), Tabelle 1, Abschnitte 2.1–2.2; Gl. (2), (3), (8), (9) nachgerechnet (`docs/spezifikation.md`, Abschnitt 3) | BWPM insgesamt |
| wende2022 | M. Wende, J. Drozella, A. Toulouse, A. M. Herkommer, „Fast algorithm for the simulation of 3D-printed microoptics based on the vector wave propagation method“, Opt. Express 30(22), 40161–40173 (2022) | bibliografisch (aus dem Literaturverzeichnis von wende2023); Volltext fehlt | FPWPM, Herleitung von T |
| fertig2010 | M. Fertig, K.-H. Brenner, „Vector wave propagation method“, J. Opt. Soc. Am. A 27(4), 709–717 (2010) | bibliografisch (wie oben) | VWPM |
| fertig2011 | M. W. Fertig, „Vector wave propagation method“, Dissertation, Universität Mannheim (2011) | bibliografisch (wie oben) | bidirektionale VWPM |
| brenner1993 | K.-H. Brenner, W. Singer, „Light propagation through microlenses: a new simulation method“, Appl. Opt. 32(26), 4984–4988 (1993) | bibliografisch (wie oben) | skalare WPM |
| byrnes | S. J. Byrnes, „Multilayer optical calculations“, arXiv:1603.02720v5 | von Raytatouille übernommen (dort Volltext geprüft) | Fresnel Convention A |
| bornwolf1999 | M. Born, E. Wolf, *Principles of Optics*, 7. Aufl., Cambridge University Press (1999), Kap. 13 (Streuung an inhomogenen Medien, erste Born-Näherung) | bibliografisch | Sollwert K17 (`bwpm.born`) |
| chew1995 | W. C. Chew, *Waves and Fields in Inhomogeneous Media*, IEEE Press (1995) | bibliografisch | Spektraldarstellung der skalaren und der dyadischen Green-Funktion (K17, Faktor cos²θ für TM) |
| moharam1995 | M. G. Moharam, E. B. Grann, D. A. Pommet, T. K. Gaylord, „Formulation for stable and efficient implementation of the rigorous coupled-wave analysis of binary gratings“, J. Opt. Soc. Am. A 12(5), 1068–1076 (1995) | bibliografisch | RCWA des Lamellengitters (`tests/rcwa_reference.py`, Prüfung der Born-Referenz) |
| li1996 | L. Li, „Use of Fourier series in the analysis of discontinuous periodic structures“, J. Opt. Soc. Am. A 13(9), 1870–1876 (1996) | bibliografisch | Faktorisierung für TM in der RCWA |
