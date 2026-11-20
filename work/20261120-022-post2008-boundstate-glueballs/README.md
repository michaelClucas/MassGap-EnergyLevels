# 022: Does a Result B already exist in the post-2008 literature?

## Question
Result B (OVERSEER_PROMPT.md): a model built on an electron-level mechanism that, with k fitted parameters and k < number of reproduced quantities, predicts glueball ratios within lattice uncertainties, matches or beats existing models, and states what plays electron, lattice and potential. Skeptic 019 checked 7 models from 1976-2008 (none passes; best: Kaidalov-Simonov 3/10 within 1 sigma, chi^2/n 98; MSB helicity B 0/9, chi^2/n 75; both ~3.3 with a 5% budget). Here: the post-2008 bound-state models, first of all the functional Bethe-Salpeter (BSE) glueballs, which are a positronium-type bound-state equation and so count as electron-level.

Setting: model numbers are phenomenology (continuum, infinite volume, 4D SU(3), 4D SU(inf), 2+1D SU(inf)). Lattice comparison data: continuum-extrapolated AT2020 SU(3), AT2021 SU(3) and SU(inf), AT2016 2+1D SU(inf), all read from `data/`. Nothing here is a proof.

## Method
- Search: arXiv API keyword queries on abstracts (glueball AND Bethe-Salpeter / Dyson-Schwinger / Coulomb gauge / holographic+lattice / constituent gluon, helicity, flux tube; author queries Huber, Buisseret, Chevalier, Szczepaniak, Reinhardt), sorted by date, through 2026-11. INSPIRE was not used.
- Every arXiv ID below was verified via the arXiv API (title and authors match). Numbers were read from the LaTeX source tables of each e-print (downloaded to a scratch folder, not committed), not from abstracts.
- `models.csv`: one row per model number, with table reference.
- `score.py` (log `score.log`, output `comparison.csv`): loads lattice data with `load()` imported from `work/20261120-019-skeptic-c12/model.py`; pull formulas copied from 019 `compare.py` (ratio m/m(0++), lattice ratio errors in uncorrelated quadrature, pull_5% = (r - lv)/hypot(le, 0.05 lv)); summary metrics (n, number with |pull| <= 1, chi^2/n) copied from 019 `rank.py`. Models in sqrt(sigma) units are also compared directly on m/sqrt(sigma). FBSE MeV values are converted with r0 and r0 sqrt(sigma) from `data/string_tension_scales_su3_4d.csv` (the FBSE papers quote masses with the same r0 = 0.472 fm as AT2020).

## Per-model results

Counts are model predictions compared with a lattice state present in `data/` (0++ ground state is the normalization in ratio tests). "lat" = lattice errors only; "5%" = 5% model budget added.

| model (arXiv) | mechanism; electron / lattice / potential | fitted k (for ratios) | lattice set | test | n | within 1 sigma | chi^2/n lat | chi^2/n 5% |
|---|---|---|---|---|---|---|---|---|
| FBSE Huber-Fischer-Sanchis-Alepuz 2020/2021/2025 (2004.00415, 2110.09180, 2503.03821; Table 1 each) | 2-gluon (+ghost) BSE, Landau gauge; electron = dressed gluon/ghost propagators from a self-contained 3PI DSE system; no lattice (continuum); potential = BSE kernel (1-particle exchanges, + two-loop diagrams in 2025) | 0. Scale: one lattice number (peak of the lattice gluon dressing function). Discrete choices: truncation, decoupling-family member (results invariant), Schlessinger extrapolation | AT2020 SU(3) | ratio | 9 | 1 (4++, a * state) | 104 | 5.44 |
| same | | scale from lattice propagator | AT2020 | m(0++)/sqrt(sigma) = 3.71 vs lattice | 1 | 0 (pull +14.6, +9%) | 213 | 3.2 |
| same, with own errors (extrapolation only) added | | | AT2020 | ratio | 9 | 5 (6 within 2 sigma; 2-+ -3.0, 3++ -3.5, 0++* -2.2) | 3.1 | - |
| FBSE Sanchis-Alepuz-Fischer-Kellermann-von Smekal 2015 (1503.06051, Table 1) | same class, earlier truncation | 0 | AT2020 | ratio | 1 | 0 (0-+/0++ = 2.76, pull +74) | 5529 | 235 |
| Helicity constituent, Chevalier-Mathieu 2025 (2503.15146, Table III, DGA) | spinless Salpeter, 2 massless helicity-1 gluons; electron = transverse gluon; potential = (9/4) sigma r - 3 alpha_s/r | 0 new (sigma, alpha_s carried over from MSB 0802.0088; no instanton term) | AT2020 | ratio | 8 | 0 | 67 | 3.35 |
| Constituent + Casimir, Buisseret-Chevalier-Mathieu-Semay 2025 (2509.09454, Table III) | 2- and 3-gluon spinless Salpeter, helicity gluons, adjoint string (Casimir scaling), Y-junction | 3 (alpha_bar, C_gg, C_ggg; least squares to SU(inf) lattice) | AT2021 SU(inf) | absolute m/sqrt(sigma) | 11 | 0 | 124 | 2.58 |
| same | | | AT2021 SU(inf) | ratio | 10 | 0 | 80 | 3.08 |
| same, extended to SU(3) by its Casimir law (0 extra parameters) | | | AT2020 SU(3) | absolute / ratio | 11 / 10 | 1 / 1 | 73 / 61 | 4.25 / 4.53 |
| Casimir law alone, M(3)/M(inf) = 3/sqrt(8) (2509.09454 eq. model1) | N-dependence only | 0 | AT2021 N=3 vs inf | per-state ratio | 12 | 5 | 7.6 | 0.46 |
| Holography, Klebanov-Strassler, Dymarsky-Melnikov 2022 (2206.14826, Tables 4-5) | supergravity modes of a SUSY cascading theory; Schrodinger-type radial problem; no electron map | 0 | AT2021 SU(inf) | ratio | 5 | 0 | 161 | 3.36 |
| 2+1D constituent, Buisseret-Mathieu-Semay 2013 (1301.3247, Tables I-II) | 2- and 3-gluon massive transverse constituents, adjoint string | 1 (m_g/sqrt(sigma) = 0.7, set by hand) | AT2016 2+1D SU(inf) | absolute / ratio | 16 / 15 | 2 / 1 | 1332 / 1291 | 3.68 / 4.72 |
| Holography, NN Einstein-dilaton, Chen-Chen-Zhou 2025 (2507.06729, Table 1) | bottom-up gravity dual; no electron map | 0++ and 0++* are inputs | AT2020 | ratio | 1 (0++* is a fit input: 1.54 vs 1.72) | 0 | - | - |

Reference (019, pre-2008): MSB helicity B 0/9, chi^2/n 75 (5%: 3.26); Kaidalov-Simonov 3/10, 98 (3.34).

Checked but not scored (no testable table, or not pure-gauge lattice):
- Meyers-Swanson 2012, 1211.4648: BSE with lattice gluon propagator as input; results only in a figure; coupling g(1 GeV) fitted to glueball masses; authors state the couplings needed for quarks, ghosts and glueballs differ by ~10x. Does not meet.
- Kaptari-Kampfer 2020, 2004.06523: rainbow-ladder BSE, pseudoscalar only, effective parameters adjusted to lattice. No 0++, so no ratio. Does not meet.
- Ballon-Bayona et al. 2017, 1708.08968: dilaton-gravity; 2 parameters fitted to 0++, 0++* of Meyer 2004 (not in `data/`), other masses only in figures. Not scored.
- Rinaldi-Vento graviton soft wall (2002.11720), Sonnenschein-Weissman rotating strings (1507.01604): fit to old lattice / experimental f0 states. Not scored.
- Coulomb-gauge variational (Szczepaniak, Reinhardt): only post-2008 paper found is 1208.0882 (finite-T thermodynamics, no new spectrum table).
- Chevalier-Mathieu review 2026, 2609.16790: compiles the FBSE and lattice numbers above; says the FBSE tensor is high and 2-+ low relative to lattice; no new model.
- Carmo Terin 2026, 2601.19727 (RGZ BSE): single 0++ mass; cannot satisfy k < n with more than one quantity. Not scored.

## Verdict per model
- **FBSE (Huber-Fischer-Sanchis-Alepuz): does not meet Result B, but it is the closest.** It satisfies the parameter count (0 fitted parameters for 9 ratios; absolute scale from one lattice propagator number) and has a clean electron-level map (positronium-type BSE: dressed gluon = electron, kernel = potential, continuum, no lattice). It fails "within lattice uncertainties": 1/9 ratios within 1 lattice sigma, chi^2/n 104; its own extrapolation errors (~10%, which the authors call a lower bound) are ~5-10x the lattice errors; even with them, 2-+ (-3.0 sigma), 3++ (-3.5, both lattice and model mark 3++ *) and 0++* (-2.2) miss. m(0++)/sqrt(sigma) is 9% high. No C = - states, no N-dependence. Against the old models it does not win: chi^2/n 104 vs 75 (MSB) and 5.4 vs 3.3 at 5% budget.
- **BCMS 2025 constituent + Casimir: does not meet.** k = 3, n = 11, 0/11 within 1 sigma at SU(inf). Its N-dependence law is parameter-free and roughly right (5/12 within 1 sigma, chi^2/n 0.46 at 5%), but it predicts N-independent ratios, while the lattice m(2++)/m(0++) moves from 1.437(11) at N=3 to 1.497(8) at N=inf (AT2021). Same constituent-gluon class as C12; 1+- and 1-- degenerate at 6.0 sqrt(sigma) vs lattice 5.76 and 7.26.
- **Chevalier-Mathieu 2025 helicity: does not meet** (0/8, lowest lattice chi^2/n of the post-2008 set at 67, essentially equal to MSB 2008 on which it is built).
- **Holographic models: do not meet**, and have no electron/lattice/potential map beyond "Schrodinger problem in a radial potential". KS: 0/5. Bottom-up fits use the 0++ levels as input.
- **2+1D constituent (BMS 2013): does not meet** (2/16; 0-+ three-gluon state off by 33%).

## Final verdict
No published post-2008 model checked here meets Result B, so the "known result" outcome does not apply. Best model: functional BSE (Huber-Fischer-Sanchis-Alepuz, 2503.03821 Table 1), which meets every criterion except accuracy. Distance from the bar: it would need its 0++*, 2-+ and 3++ ratios to move by +24%, +26% and +25% (to lattice), and its uncertainty to shrink from ~10% to the lattice ~1-2%; at present 1/9 ratios within 1 lattice sigma (lattice chi^2/n 104). The best chi^2 among all post-2008 models is 61-80 per prediction (BCMS, Chevalier-Mathieu), the same as the best pre-2008 constituent models (75, 98). Every bound-state model, old or new, sits at the 5-10% level; none reaches the ~1% lattice errors. The C12 statement "none of the published models checked fits within lattice errors" now extends to post-2008 functional, constituent and holographic models.

## Open issues
- FBSE N-dependence and C = - (three-gluon) states are not computed in the literature; whether the BSE ratios move with N as the lattice does (m(2++)/m(0++) 1.437 -> 1.497) is an open, testable sub-question.
- FBSE errors cover only the time-like extrapolation; truncation errors are not quantified beyond "marginal shifts" from two-loop terms.
- Ratio errors use uncorrelated quadrature (019 convention), which makes the within-1-sigma counts generous.
- The FBSE 2021 -> 2025 tables are not a pure rescale (4-+ changed from 3240(300) to 5300(600) MeV); only 2025 is scored.
- Literature search was arXiv-abstract keyword search only (no INSPIRE citation crawl); a model with only figures, or published outside arXiv, could be missed. Bottom-up holographic fits (1708.08968, Rinaldi-Vento) were not digitized.
- KS holography: multiplet 8 contains a 0-+ degenerate with the 0++; the authors' J^PC assignment was used.

## Reproduce
From the repo root (1 thread, < 100 MB, < 5 s):
```
python tools/pin.py .venv/Scripts/python.exe work/20261120-022-post2008-boundstate-glueballs/score.py > work/20261120-022-post2008-boundstate-glueballs/score.log
```
