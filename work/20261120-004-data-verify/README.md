# 20261120-004 data verify

## Question

Are the lattice glueball numbers in `data/` (extracted by script in `work/20261120-002-glueball-data/`)
correctly transcribed from the seven cited papers, are the citations and table numbers in
`data/SOURCES.md` right, and do the headline ratios and reported discrepancies in
`work/20261120-002-glueball-data/README.md` hold up?

## Method

Independent of the extraction script: the source here is the arXiv **PDF** (not the LaTeX), fetched
2026-11-20 from `https://arxiv.org/pdf/<id>` (versions: hep-lat/9901004v2, hep-lat/0510074v1,
2007.06422v2, hep-lat/0404008v1, 2106.00364v3, hep-lat/9804008v2, 1609.03873v1).
`extract_tables.py` and `ratios.py` were not run, imported or read for logic.

1. Metadata: arXiv ID, title, authors and posting date read from page 1 of each PDF; table
   captions and printed table numbers located in the PDF text.
2. Values, two independent ways:
   - By eye: cited tables rendered as page images (or read from PDF text) and compared row by row
     with the CSV: MP1999 VII, VIII; Chen2006 VIII, XXI; AT2020 17, 18; LTW2004 14; AT2021 34, 35,
     38; Teper1998 28, 33, 34; AT2016 B1, B4, B10-B12 (CSV "Tables 18, 21, 27-29"). AT2016 B3,
     B5-B9 (CSV "Tables 20, 22-26") were viewed but read by eye only for headline and sampled rows.
     Star marks checked for LTW2004 and AT2021.
   - By machine (`token_check.py`): `pdftotext -layout` text of each paper; for every CSV row with an
     error, the printed token `value(err)` (last-digit notation, or literal `1.7(1.5)` form) must
     occur on the PDF page holding that table. This confirms the number string, not its row/column.
3. Random sample (`sample_rows.py`): all 808 rows of `data/glueball_*.csv`, minus the 0++/2++/0-+
   ground-state rows and scale rows (checked exhaustively), pooled; `random.seed(0)`,
   `random.sample(pool, 30)`. Each sampled row checked by eye on the page image.
4. Scales: AT2020 eqns (8), (9), ref. [21], and the last row of Table 17.
5. `verify_ratios.py` (stdlib, written from scratch) recomputes the headline ratios and the
   discrepancy tensions from the CSVs. Errors: stat and sys in quadrature, inputs uncorrelated.

Rerun (repo root, after putting the PDF texts somewhere outside the repo):
`python tools/pin.py python work/20261120-004-data-verify/verify_ratios.py`,
`python tools/pin.py python work/20261120-004-data-verify/sample_rows.py`,
`python tools/pin.py python work/20261120-004-data-verify/token_check.py <csv> <pdftotext.txt> "Table 28=69" ...`
(page maps used: MP1999 VII=12, VIII=13; Chen VIII=12, XXI=21; AT2020 17=46, 18=47; LTW 14=30;
AT2021 34=84, 35=85, 38=87; Teper 28=69, 33=73, 34=74; AT2016 18=28, 20..29=29..38; PDF page index).

## Result

### Metadata

| source | arXiv ID, title, authors | posted | table numbers in SOURCES.md | content / units |
|---|---|---|---|---|
| MP1999 | correct | Jan 1999 | VII, VIII correct | continuum r0*m and MeV with r0^-1 = 410(20) MeV; ratios |
| Chen2006 | correct (13 authors, all listed) | **Oct 2005** (journal PRD 73 (2006) 014516, as AT2021 ref. [10] cites it) | VIII, XXI correct | continuum r0*M, MeV with 410(20) MeV |
| AT2020 | correct | Jul 2020 | 17, 18 correct | 17: M/sqrt(sigma) continuum, last row r0 sqrt(sigma) = 1.160(6); 18: GeV with r0 = 0.472(5) fm, stat errors only |
| LTW2004 | correct | Apr 2004 | 14 correct | continuum m/sqrt(sigma), SU(inf) extrapolation; Table 15 is fixed a ~ 0.2/sqrt(sigma) (correctly excluded) |
| AT2021 | correct | Jun 2021 (v3 Dec 2021) | 34, 35, 38 correct | continuum M/sqrt(sigma); 38 is SU(inf) with >= lower bounds |
| Teper1998 | correct | Apr 1998 (v2 Nov 1998) | 28, 33, 34 correct | continuum m/sqrt(sigma); 33 in g^2 N units, 34 in sqrt(sigma) units |
| AT2016 | correct | Sep 2016 | **wrong numbering**, see errors | continuum M/sqrt(sigma_f), B1 in g^2 N units |

### Rows feeding the headline ratios and discrepancies (all checked by eye; all match)

| file | N | JPC | quantity | CSV value | paper value | match |
|---|---|---|---|---|---|---|
| AT2020 | 3 | 0++ / 2++ / 0-+ | m/sqrt_sigma | 3.405(21) / 4.894(22) / 5.276(45) | same (T17) | yes |
| AT2021 | 3 | 0++ / 2++ / 0-+ | m/sqrt_sigma | 3.405(21) / 4.894(22) / 5.276(45) | same (T34) | yes |
| AT2021 | 8 | 0++ / 2++ / 0-+ | m/sqrt_sigma | 3.099(26) / 4.660(20) / 4.755(58) | same (T35) | yes |
| AT2021 | inf | 0++ / 2++ / 0-+ | m/sqrt_sigma | 3.072(14) / 4.599(14) / 4.711(26) | same (T38) | yes |
| AT2021 | 2, 4, 5, 6, 10, 12 | 0++ / 2++ / 0-+ | m/sqrt_sigma | all 18 values | same (T34, T35) | yes |
| MP1999 | 3 | 0++ / 2++ / 0-+ | r0*m | 4.21(11)(4) / 5.85(2)(6) / 6.33(7)(6) | same (T VII) | yes |
| MP1999 | 3 | 0++ / 2++ / 0-+ | m [MeV] | 1730(50)(80) / 2400(25)(120) / 2590(40)(130) | same | yes |
| MP1999 | 3 | 2++, 0-+ | m/m(0++) | 1.39(4), 1.50(4) | same (T VIII) | yes |
| Chen2006 | 3 | 0++ / 2++ / 0-+ | r0*m | 4.16(11)(4) / 5.83(5)(6) / 6.25(6)(6) | same (T XXI) | yes |
| Chen2006 | 3 | 0++ / 2++ / 0-+ | m [MeV] | 1710(50)(80) / 2390(30)(120) / 2560(35)(120) | same | yes |
| LTW2004 | 3 | 0++ / 2++ | m/sqrt_sigma | 3.55(7) / 4.78(9) | same (T14) | yes |
| LTW2004 | 8 | 0++ / 2++ | m/sqrt_sigma | 3.55(12) / 4.73(22) | same | yes |
| LTW2004 | inf | 0++ / 2++ | m/sqrt_sigma | 3.307(53) / 4.80(14) | same | yes |
| LTW2004 | 2, 4, 4an, 6 | 0++ / 2++ | m/sqrt_sigma | all 8 values | same | yes |
| AT2016 | 3 | 0++ / 2++ / 0-+ | m/sqrt_sigma | 4.3683(73) / 7.241(17) / 9.234(43) | same (B4) | yes |
| AT2016 | inf | 0++ / 2++ / 0-+ | m/sqrt_sigma | 4.116(6) / 6.914(13) / 8.998(28) | same (B10, B11) | yes |
| AT2016 | 2, 4, 6, 8, 12, 16 | 0++ / 2++ / 0-+ | m/sqrt_sigma | all 18 values | same (B3, B5-B9) | yes |
| AT2016 | all, inf | 0++ ; scale | m/g2N ; sqrt_sigma/g2N | all 20 values | same (B1) | yes |
| Teper1998 | 3 | 0++ / 2++ / 0-+ | m/sqrt_sigma | 4.329(41) / 7.13(12) / 9.30(25) | same (T28) | yes |
| Teper1998 | inf | 0++ / 2++ / 0-+ | m/sqrt_sigma | 4.065(55) / 6.88(16) / 9.02(30) | same (T34) | yes |
| Teper1998 | inf | 0++ / 2++ / 0-+ | m/g2N | 0.808(11) / 1.365(33) / 1.787(60) | same (T33) | yes |

### Scale rows (`string_tension_scales_su3_4d.csv`)

| row | CSV | paper | match |
|---|---|---|---|
| MP1999 1/r0 | 410(20) MeV | T VII caption 410(20) MeV | yes |
| Chen2006 1/r0 | 410(20) MeV | T XXI caption 410(20) MeV | yes |
| AT2020 r0*sqrt_sigma | 1.160(6) | T17 last row and eqn (9) | yes |
| AT2020 r0 | 0.472(5) fm | eqn (8), = 1/418(5) MeV, from ref. [21] R. Sommer, PoS LATTICE2013 015, arXiv:1401.3270 | yes; recomputed 1/r0 = 418.1(4.4) MeV |
| AT2020 sqrt_sigma | 485(6) MeV | eqn (9) | yes; recomputed 1.160/0.472 fm = 485.0(5.7) MeV |
| AT2020 sqrt_sigma | 440 MeV, no error | text before eqn (8), "something like 440 MeV" from Regge slope | yes |

The 1.160(6) is AT2020's continuum fit of r0 sqrt(sigma) using r0/a interpolations from refs. [22, 23]
(Wilson action). Applying it to MP1999/Chen2006 (improved anisotropic action) is what AT2021 itself
does in its Table 39 (it lists Chen2006 as 3.59(15) after "our transformation"), so the conversion
is the authors' own convention, with a larger error than the 0.10 obtained here.

### Random sample (seed 0, 30 rows; all match)

| CSV line | N | JPC lvl | CSV value | paper | match |
|---|---|---|---|---|---|
| chen_2006:15 | 3 | 3++ 0 | 3670(50)(180) MeV | T XXI 3670(50)(180) | yes |
| morningstar_peardon_1999:18 | 3 | 3++ 0 | r0 m 8.99(4)(9) | T VII 8.99(4)(9) | yes |
| at2016:27 | 2 | 0-+ 2 | 11.50(5) | B3 11.50(5) | yes |
| at2016:31 | 2 | 0++ 5 | 9.910(32) | B3 9.910(32) | yes |
| at2016:34 | 2 | 2-+ 0 | 7.795(12) | B3 7.795(12) | yes |
| at2016:76 | 3 | 2-- 0 | 8.599(30) | B4 8.599(30) | yes |
| at2016:77 | 3 | 2+- 0 | 8.665(24) | B4 8.665(24) | yes |
| at2016:84 | 3 | 1++ 0 | 10.221(50) | B4 10.221(50) | yes |
| at2016:161 | 6 | 0-- 2 | 9.155(39) | B6 9.155(39) | yes |
| at2016:196 | 6 | 1-- 0 | 9.568(43) | B6 9.568(43) | yes |
| at2016:206 | 8 | 0++ 1 | 6.332(15) | B7 6.332(15) | yes |
| at2016:231 | 8 | 2+- 0 | 8.273(34) | B7 8.273(34) | yes |
| at2016:254 | 12 | 0-- 0 | 6.070(15) | B8 6.070(15) | yes |
| at2016:261 | 12 | 0++ 4 | 8.955(44) | B8 8.955(44) | yes |
| at2016:284 | 12 | 2+- 1 | 9.656(30) | B8 9.656(30) | yes |
| at2016:309 | 16 | 0++ 1 | 6.295(18) | B9 6.295(18) | yes |
| at2016:313 | 16 | 0++ 3 | 8.157(41) | B9 8.157(41) | yes |
| at2016:343 | 16 | 1++ 0 | 9.90(7) | B9 9.90(7) | yes |
| at2016:364 | inf | 0++ 5 | 9.087(20) | B10 9.087(20) | yes |
| at2016:382 | inf | 2-+ 4 | 10.914(39) | B11 10.914(39) | yes |
| teper_1998:43 | 3 | 2-- 0 | 8.75(17) | T28 8.75(17) | yes |
| teper_1998:48 | 5 | 2-- 1 | 9.79(45) | T28 9.79(45) | yes |
| teper_1998:57 | 4 | 1++ 0 | 9.91(36) | T28 9.91(36) | yes |
| teper_1998:79 | inf | 0-- 1 | 7.63(37) | T34 7.63(37) | yes |
| teper_1998:86 | inf | 1-- 0 | 9.36(60) | T34 9.36(60) | yes |
| at2021:14 | 2 | 2++ 1 | 7.22(6), mark + | T34 7.22(6)^+ | yes |
| at2021:46 | 5 | 2+- 0 | 8.49(15) | T34 8.49(15) | yes |
| at2021:78 | 6 | 3++ 0 | 7.44(8), mark * | T35 7.44(8)^* | yes |
| at2021:101 | 10 | 1-+ 0 | 8.49(28) | T35 8.49(28) | yes |
| at2021:118 | 12 | 3+- 0 | 6.89(10) | T35 6.89(10) | yes |

(File names shortened: `glueball_sun_3d_athenodorou_teper_2016.csv` = at2016, etc.)

### Machine token check (all rows)

803 of 808 rows have an error and were token-checked: MP1999 38/38, Chen2006 41/41, AT2020 41/41,
LTW2004 21/21, AT2021 153/153, Teper1998 105/105, AT2016 404/404 found on the right page. The other
5 (scale rows without a table cell) are in the scale table above.

### Headline ratios (recomputed, `verify_ratios.py`)

| dataset | m(2++)/m(0++) | m(0-+)/m(0++) | m(0++)/sqrt(sigma) | vs 002 README |
|---|---|---|---|---|
| 4D SU(3) AT2020 T17 | 1.437(11) | 1.549(16) | 3.405(21) | same |
| 4D SU(3) AT2021 T34 | 1.437(11) | 1.549(16) | 3.405(21) | same (see below) |
| 4D SU(3) MP1999 | 1.390(41) | 1.504(47) | 3.63(10) | same |
| 4D SU(3) Chen2006 | 1.401(44) | 1.502(47) | 3.59(10) | same |
| 4D SU(3) LTW2004 | 1.346(37) | n/a | 3.550(70) | same |
| 4D SU(inf) AT2021 T38 | 1.497(8) | 1.534(11) | 3.072(14) | same |
| 4D SU(inf) LTW2004 | 1.451(48) | n/a | 3.307(53) | same |
| 3D SU(3) AT2016 B4 | 1.658(5) | 2.114(10) | 4.368(7) | same |
| 3D SU(3) Teper1998 | 1.647(32) | 2.148(61) | 4.329(41) | same |
| 3D SU(inf) AT2016 B10/B11 | 1.680(4) | 2.186(8) | 4.116(6) | same |
| 3D SU(inf) Teper1998 | 1.692(46) | 2.219(80) | 4.065(55) | same |

AT2020 = AT2021 for SU(3) is correct, not a duplication: the 0++, 2++, 0-+ entries of AT2021 Table 34
equal AT2020 Table 17 digit for digit, and AT2021 says "most of our SU(3) calculations ... have
recently been published separately [11]" with [11] = arXiv:2007.06422. Same underlying data.

### Discrepancies (all real, not transcription)

| comparison | difference | tension |
|---|---|---|
| m(0++)/sqrt(sigma) SU(3): MP1999 3.63(10) vs AT2020 3.405(21) | +0.224 | 2.1 sigma |
| same, Chen2006 3.59(10) vs AT2020 | +0.181 | 1.7 sigma |
| same, LTW2004 3.55(7) vs AT2020 | +0.145 | 2.0 sigma |
| m(2++)/m(0++) SU(3): LTW2004 1.346(37) vs AT2020 1.437(11) | -0.091 | 2.4 sigma |
| m(0++)/sqrt(sigma) SU(8): LTW2004 3.55(12) vs AT2021 3.099(26) | +0.451 | 3.7 sigma |
| m(0++)/sqrt(sigma) SU(inf): LTW2004 3.307(53) vs AT2021 3.072(14) | +0.235 | 4.3 sigma |

All input numbers are confirmed against the papers. The authors print the same comparisons:
AT2021 Table 39 (SU(3): this work 3.405(21); Chen2006 converted 3.59(15); LTW2004 3.55(7), 4.78(9))
and Table 40 (SU(8): this work 3.099(26); LTW2004 3.55(12), 4.73(22)). Also confirmed: Chen2006
Table VIII prints MP1999's 1-- as 9.50(4) (MP1999 T VII: 9.40(6)); AT2016 B4 prints the 2-- row
12.09(12) 5.7/3 twice; AT2020 2+- gs 8.74(12)* vs AT2021 8.91(15).

All settings: lattice Monte Carlo, continuum-extrapolated (and large-N where N = inf) by the cited
authors, 3+1D or 2+1D as labelled, SU(N). Observations, not proofs.

## Verdict

**Values verified: 0 transcription errors** in 808 rows (every row token-matched on the right page;
every headline row, every scale row and 30 random rows matched by eye; full tables compared by eye
row by row except AT2016 B3, B5-B9, where only headline and sampled rows were read by eye). Headline ratios and discrepancy figures reproduce exactly.

Errors in `data/SOURCES.md` metadata (overseer to fix):

1. **AT2016 table numbers.** SOURCES.md says table numbers equal "the printed numbering of the arXiv
   version". For 1609.03873v1 that is false: the printed numbers are appendix numbers. Mapping:
   CSV Table 18 = printed **B1**; Table 19 = B2; Tables 20-26 = **B3-B9** (SU(2), 3, 4, 6, 8, 12, 16);
   Tables 27, 28, 29 = **B10, B11, B12**. The same wrong numbers are in the `table` column of
   `data/glueball_sun_3d_athenodorou_teper_2016.csv` and in the 002 README ("AT2016 Tab 21",
   "Tabs 27-28"). The LaTeX labels in the CSV are fine. The other six papers' table numbers are correct.
2. **Chen2006 year.** The arXiv version was posted Oct 2005; "2006" is the journal year
   (Phys. Rev. D 73 (2006) 014516, per AT2021's bibliography). SOURCES.md labels it "2006 (arXiv
   version)". Minor; the source_id can stay.

## Open issues

- Error on the MP1999/Chen2006 m(0++)/sqrt(sigma) conversion: 0.10 here (quadrature of r0 m errors and
  1.160(6)); AT2021 Table 39 quotes 0.15 for Chen2006. Origin of the larger error not checked.
- Dagger marks in MP1999 Table VII and the AT2016 parity-doubling statement were not re-checked.
- Journal references for MP1999, LTW2004, Teper1998, AT2016 were not verified (only AT2020 and
  Chen2006, via the AT2021 bibliography).
