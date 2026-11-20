# 20261120-002 glueball data

## Question

What do lattice simulations say, with errors and traceable sources, about the pure SU(N)
Yang-Mills glueball spectrum (3+1D SU(3), 3+1D SU(N) up to large N, 2+1D SU(N)) and the string
tension scale? These are the target numbers for any Result B model.

## Method

- Downloaded the arXiv LaTeX e-prints of hep-lat/9901004, hep-lat/0510074, 2007.06422,
  hep-lat/0404008, 2106.00364, hep-lat/9804008, 1609.03873.
- `extract_tables.py` parses the continuum (and large-N) tables directly from the LaTeX and writes
  `data/*.csv` (schema in `data/README.md`, provenance in `data/SOURCES.md`). No numbers typed by
  hand. Spot-checked rows against the printed tables and row counts against table sizes.
- `ratios.py` (stdlib) computes headline ratios from the CSVs. Errors: stat and sys in quadrature,
  inputs treated as uncorrelated. For MP1999/Chen2006, m/sqrt(sigma) uses r0*sqrt(sigma) = 1.160(6)
  from AT2020.

Rerun: `python work/20261120-002-glueball-data/ratios.py` (via `tools/pin.py`).

## Result

All entries: lattice, continuum-extrapolated by the cited authors, infinite volume (authors' finite
volume checks), pure gauge. Observations, not proofs.

| setting / source | m(2++)/m(0++) | m(0-+)/m(0++) | m(0++)/sqrt(sigma) |
|---|---|---|---|
| 3+1D SU(3), AT2020 Tab 17 | 1.437(11) | 1.549(16) | 3.405(21) |
| 3+1D SU(3), MP1999 Tab VII (sqrt(sigma) via AT2020 r0*sqrt(sigma)) | 1.390(41) | 1.504(47) | 3.63(10) |
| 3+1D SU(3), Chen2006 Tab XXI (same conversion) | 1.401(44) | 1.502(47) | 3.59(10) |
| 3+1D SU(3), LTW2004 Tab 14 | 1.346(37) | n/a | 3.550(70) |
| 3+1D SU(inf), AT2021 Tab 38 | 1.497(8) | 1.534(11) | 3.072(14) |
| 3+1D SU(inf), LTW2004 Tab 14 | 1.451(48) | n/a | 3.307(53) |
| 2+1D SU(3), AT2016 Tab B4 | 1.658(5) | 2.114(10) | 4.368(7) |
| 2+1D SU(3), Teper1998 Tab 28 | 1.647(32) | 2.148(61) | 4.329(41) |
| 2+1D SU(inf), AT2016 Tabs B10-B11 | 1.680(4) | 2.186(8) | 4.116(6) |
| 2+1D SU(inf), Teper1998 Tab 34 | 1.692(46) | 2.219(80) | 4.065(55) |

(AT2021 Tab 34 SU(3) gives the same three numbers as AT2020.) MP1999's own Table VIII ratios,
1.39(4) and 1.50(4), match the recomputed values. In 2+1D the sigma is the fundamental string
tension sigma_f, and J=2 is parity doubled, so m(2-+) = m(2++) within errors there.

Physical-unit anchors (3+1D SU(3)): r0^-1 = 410(20) MeV (MP1999, Chen2006 convention);
r0 = 0.472(5) fm and sqrt(sigma) = 1.160(6)/r0 = 485(6) MeV (AT2020); "sqrt(sigma) ~ 440 MeV" from
Regge slopes, stated as rough (AT2020 text). m(0++) ~ 1.65-1.73 GeV depending on convention.

## Verdict

Data set is complete for the requested papers. Robust targets for a model: in 3+1D,
m(2++)/m(0++) ~ 1.4-1.5, m(0-+)/m(0++) ~ 1.5, m(0++)/sqrt(sigma) ~ 3.1 (N=inf) to 3.4-3.6 (SU(3)),
all with weak N-dependence absorbed by 1/N^2. In 2+1D the ratios differ clearly (2++/0++ ~ 1.66,
0-+/0++ ~ 2.1-2.2), so a model has to get the dimension dependence right, not just one number.

## Open issues / discrepancies

1. 3+1D SU(3) m(0++)/sqrt(sigma): AT2020 3.405(21) vs MP1999 3.63(10) and Chen2006 3.59(10)
   (after r0 -> sigma conversion with AT2020's r0*sqrt(sigma)), about 2 sigma. AT2021 Table 39
   lists other older values (3.35-3.61). The conversion factor is from a different action, so part
   of this may be scale-setting rather than spectrum.
2. 3+1D SU(3) m(2++)/m(0++): AT2020 1.437(11) vs LTW2004 1.346(37) (2.4 sigma) and MP1999 1.39(4).
3. 3+1D SU(8) m(0++)/sqrt(sigma): LTW2004 3.55(12) vs AT2021 3.099(26); the newer large-N value
   3.072(14) is lower than LTW2004's 3.307(53). Prefer AT2021 (larger N range, more spacings), but
   flag the shift.
4. AT2020 vs AT2021 SU(3): only 2+- gs differs, 8.74(12)* vs 8.91(15).
5. Chen2006 Table VIII misquotes MP1999's 1-- as 9.50(4) (MP1999: 9.40(6)). AT2016 Table B4
   duplicates the SU(3) 2-- level-3 row. Both recorded, neither changes the headline ratios.
6. Ratio errors ignore correlations (common sqrt(sigma) or r0); true errors on same-paper ratios
   are somewhat smaller.
7. Not extracted (available in the sources if needed): cubic-irrep tables, fixed-spacing tables,
   k-string tension ratios, topological susceptibility, AT2021 N-independent lower bounds.


## Erratum (2026-11-20)
AT2016 table numbers were the LaTeX source order (18-29); the printed arXiv labels are B1-B12. Corrected in the CSV, SOURCES.md and ratios.py after independent verification (work/20261120-004-data-verify). Values unchanged.
