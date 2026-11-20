# 019: Skeptic/phenomenologist audit of C12 (constituent-gluon Result-B attempt, work 013)

## Question
Does C12 survive? Specifically: (1) are the competitor-model numbers in `013/competitors.csv` correct against the primary sources, with honest parameter counts; (2) does the comparison and ranking (Kaidalov-Simonov best, 3/10 within 1 sigma) hold with verified numbers against Athenodorou-Teper 2020 (AT2020) ratios from `data/`; (3) is the model fit reproducible and could the negative result be an implementation bug; (4) is "Result B not met" fair, and is "no electron-level model in the literature fits within lattice errors" supported?

Setting throughout: continuum-extrapolated lattice data, infinite volume, 4D SU(3) (AT2020 Table 17, via `data/`), plus 013's 4D SU(inf) and 2+1D SU(3)/SU(inf) fits. Model numbers are phenomenology, not proofs.

## Method
- Downloaded arXiv PDFs (IDs verified from the title pages): 0802.0088 (Mathieu-Buisseret-Semay, "Gluons in glueballs: spin or helicity?"), hep-ph/9911291 (Kaidalov-Simonov), hep-th/0003115 (Brower-Mathur-Tan), 0810.4453 (Mathieu-Kochelev-Vento review "The physics of glueballs"), 0803.0815 (Mathieu-Semay three-gluon), hep-ph/0308268 (Szczepaniak-Swanson 2003), hep-ph/9511422 (Szczepaniak-Swanson-Ji-Cotanch 1996). Text extracted with pdftotext; flux-tube Table V alignment cross-checked on the ar5iv HTML of 0810.4453. PDFs are not committed.
- Re-extracted every row of `013/competitors.csv` into `competitors_verified.csv` with a `status_vs_013` column.
- `compare.py`: 013's script unchanged except it reads `competitors_verified.csv` and writes `comparison_verified.csv` (log `compare.log`). `rank.py`: ranks by several metrics including a common subset (`rank.log`).
- `fit.py`, `model.py`: byte copies of 013's; rerun here (`run.log`, `fits.csv`).
- `limits.py`: tests of `model.py` in limits with known answers (`limits.log`).

## 1. Competitor numbers (verified table)

| model | source checked | rows in 013 | verdict |
|---|---|---|---|
| MSB helicity Model B | 0802.0088 Table III | 9 | all **correct** (Table III is cleanly readable in the PDF; the 013 "garbled" note was an ar5iv artifact). Omitted: 4++ = 3.768 GeV (added). |
| Kaidalov-Simonov | hep-ph/9911291 Table 6 | 11 | all **correct**. Omitted 3-- = 4.03 (no AT2020 3-- state, so no effect). |
| "SSJC Coulomb gauge PRL 76 (1996)" | 0802.0088 Table III col. "CGQCD [8]"; ref [8] = Szczepaniak-Swanson PLB 577 (2003) 61, hep-ph/0308268, Table 1 | 9 | values **correct, attribution wrong**: they are Szczepaniak-Swanson 2003, not SSJC 1996 (which shows only a figure, no table). Relabelled `SS2003-CoulombGauge`. Omitted 4++ = 3.99 (added). |
| Brower-Mathur-Tan AdS7 | hep-th/0003115 Tables 1-2 | 6 | 5 correct, **one wrong**: 0++ level 1 must be the T4 0++ degenerate with 2++ at m = 4.7007 (Table 1; text: dilaton 0++ "remains degenerate with the tensor"), not S4 n=1 = sqrt(46.986) = 6.8547. Pull on 0++*/0++ goes from +50.9 to +1.2. Added 0-+* = sqrt(72.489) = 8.514. |
| Isgur-Paton flux tube | 0810.4453 Table V (primary PRD 31 (1985) 2910 not on arXiv, not checked) | 7 | all **correct vs the review**. |
| Jaffe-Johnson bag | 0810.4453 text | 5 | all **correct vs the review**; omitted (TE)^3 0++ and 3+- at 1.45 GeV (added; the 0++ one is the bag's 0++*). Primary PLB 60 (1976) not checked. |
| hard wall D/N | 0810.4453 Table VIII | 6 | all **correct**. |

Parameter counts (what each paper says):
- Kaidalov-Simonov: "no fitting parameters"; alpha_s = 0.3 chosen; sigma_f = 0.238 GeV^2 (Table 6 caption) / 0.228 (text) chosen "for comparison with lattice", vs 0.18 in Table 5. Their ratios are **not** scale-free: 2++/0++ = 1.71 in Table 5 and 1.64 in Table 6, so the sigma choice is an effective tunable input for ratios. "0 glueball-fitted" is generous; 1 hidden choice is fairer.
- MSB Model B: sigma = 0.185 GeV^2 (meson value), alpha_s = 0.45 chosen near refs [7,8], I = 0.45 GeV (instanton term -P I delta_J0) set to reproduce lattice 0++ and 0-+. For ratios: 1 fitted (I) + 1 chosen (alpha_s). 013's count is right in substance.
- SS2003: "no free parameters"; scale and RG parameters fixed by the Wilson-loop static potential. Correct.
- BMT: ratios parameter-free; Fig. 2 sets the scale on 2++. Correct.
- Isgur-Paton: parameters hand-chosen; the review says splittings are insensitive to them. Correct.
- Bag: B sets the scale only (alpha_s = 0 in Jaffe-Johnson). Correct.

## 2. Comparison recomputed (AT2020 ratios m/m(0++) from `data/`)

| model | ratio params | n | within 1 sigma | chi^2/n lat | chi^2/n 5% |
|---|---|---|---|---|---|
| MSB helicity B | 1 fitted + 1 chosen | 9 | 0 | 75 | 3.26 |
| Kaidalov-Simonov | 0 fitted, sigma and alpha_s chosen | 10 | 3 (0++*, 3++, 2--) | 98 | 3.34 |
| this work rel k=3 | 3 fitted (2 of them on ratios here) | 11 | 2 (fitted) / 0 predicted | 70 | 5.43 |
| SS2003 Coulomb gauge | 0 | 9 | 0 | 274 | 12.5 |
| flux tube (IP) | hand-chosen | 5 | 0 | 447 | 12.1 |
| bag (JJ) | 0 | 6 | 0 | 461 | 19.0 |
| BMT AdS7 | 0 | 6 | 0 | 1094 | 49.9 |
| hard wall Dirichlet | 0 (0++ input) | 2 | 0 | 20 | 0.59 |

Changes vs 013: MSB 0/8, chi^2/n 84 -> 0/9, 75 (4++ added, pull -1.2); BMT 1705 -> 1094 (0++* fixed; still 0 within 1 sigma, closest pull +1.2); bag 580 -> 461 (0++*, 3+- added); SS2003 relabelled, 304 -> 274. KS unchanged.

Ranking: "KS best" holds only for the count within 1 sigma (3/10). By lattice chi^2/n MSB is better (75 vs 98); with the 5% budget they tie (3.26 vs 3.34). On the common set {2++, 0-+, 0++*}: MSB 0/3 chi^2/n 90, KS 1/3 119. KS's 3 hits include 3++, which AT2020 marks * (J assignment uncertain). Fair statement: KS and MSB are the two best checked models and are statistically indistinguishable at the 5% level; neither is within lattice errors. The constituent-gluon conclusion of 013 is unaffected.

## 3. Our model: reproduction and bug search
- `fit.py` rerun: identical `fits.csv` (2 cells differ only as "-0.0" vs "0.0"). Best fit rel k=3: m = 3.1323 sqrt(sigma), kappa = 0.7187, C0 = -3.3044. Counts reproduced: 4D SU(3) 0/11 (chi^2/pred 94.0), SU(inf) 0/12 (196.4), 2+1D SU(3) 0/13 (764.5), SU(inf) 0/15 (1305.3). Minor: the unreported 2+1D nr variant has 1/15 at SU(inf); does not change anything.
- `limits.py`, all 34 checks pass:
  - NR linear, 2m + p^2/m + s r: E - 2m = (s^2/m)^(1/3) x {2.33811 (Airy zero), 3.36125 (l=1), 4.24818 (l=2)} to < 2e-4, at m = 1 and at the fit point m = 3.13, s = 9/4.
  - Relativistic 2 sqrt(p^2+m^2) + r, L = 0, m = 0, 1, 3.13: matches an independent sine-basis solver (exact |p| in the box) to < 2e-4 for n = 0, 1. Massless: 3.1574 vs sqrt(2) x 2.2322 = 3.1568 (literature value for |p| + r quoted from memory, not re-fetched).
  - Tensor operator: S = 0 gives 0; sum_J (2J+1)<T> = 0 for (L,S) = (1,1), (2,2), (3,1); 3P_J pattern -2 : 1 : -0.2 has the spin-1/2 ratios -4 : 2 : -2/5; sign within S = 1 same as for quarks.
  - Smeared OGE: analytic V', V'' match finite differences; grad^2 V_V = 4 pi kappa delta_s exactly; delta_s normalized.
  - Colour: C_A/C_F = 9/4 (N=3), 2 (N -> inf); adjoint Coulomb strength N alpha_s; sum S_i.S_j values; hyperspherical K = 0 averages <cos a> = 32/(15 pi), <1/cos a> = 16/(3 pi).
- Not testable as a bug: the Breit-Fermi coefficients for spin-1 constituents are borrowed from spin-1/2 (013 lists this as a discrete choice); first-order PT with spin terms of 0.5-1 sqrt(sigma). These are modelling choices, and 0802.0088 shows the spin formalism's spurious 1-+ persists in a different implementation (their Model A puts 1-+ at 2.50 GeV, degenerate with 0-+). So the negative result is not an implementation artifact.

## 4. Verdict
- "Result B not met": **correct.** 0 predicted states within 1 lattice sigma in every setting, reproduced; implementation passes every limit test; the failure modes (spurious 1-+, C = - ordering) are also in the published independent implementations (0802.0088 Model A; 0803.0815 Table IV: 1-- 3999 MeV below 1+- 4626 MeV, verified).
- "Kaidalov-Simonov does better (3/10 within 1 sigma)": numbers correct, but "best" is metric-dependent (MSB ties or beats it), and KS ratios depend on a chosen sigma.
- "No electron-level model in the literature fits within lattice errors": **not supported as worded.** It is supported only for the 7 published models checked, all from 1976-2008 and none tuned to AT2020. Not checked: functional bound-state (Bethe-Salpeter) calculations, e.g. Huber-Fischer-Sanchis-Alepuz arXiv:2004.00415 (abstract claims "very good quantitative agreement" with lattice for 0++, 0-+ and excitations; ID verified, numbers not checked), their J = 2 follow-up, improved holographic QCD, later Coulomb-gauge and flux-tube work. A Bethe-Salpeter bound state is arguably an electron-level (positronium-type) mechanism, so this gap matters for the Result-B search.
- Fixes needed in 013: SSJC attribution, BMT 0++* level, omitted rows; the "garbled MSB Table 3" caveat can be dropped.

## Replacement CLAIMS.md wording (C12)
| C12 | Result-B attempt: constituent-gluon bound states (electron = massive spin-1 gluon, potential = adjoint Cornell, hyperfine = Breit-Fermi) | 4D SU(3), SU(inf); 2+1D SU(3), SU(inf); continuum lattice data from data/ | 3 fitted parameters (4D) / 2 (2+1D) plus ~10 discrete choices. 0/11, 0/12, 0/13, 0/15 predictions within 1 lattice sigma; chi^2/prediction 94, 196, 765, 1305 (lattice errors); reproduced, and the solver passes limit tests (skeptic 019). Structural failures: spurious light 1-+, no 1+- state, too-weak N-dependence; the first two also appear in published implementations (0802.0088 Model A, 0803.0815). Not new (same class as Kaidalov-Simonov, Mathieu-Buisseret-Semay). Against AT2020 ratios with competitor numbers verified from primary sources (019), none of the 7 published models checked (bag, flux tube, Kaidalov-Simonov, MSB helicity B, Szczepaniak-Swanson 2003, Brower-Mathur-Tan, hard wall; 1976-2008) has all compared states within 1 lattice sigma; Kaidalov-Simonov has the most within 1 sigma (3/10, chi^2/n 98), MSB helicity B the lowest chi^2/n (75, 0/9); with a 5% model budget both are at chi^2/n ~3.3, this model 5.4. Post-2008 models (functional Bethe-Salpeter, holographic) not checked. Result B not met. | refuted (as a Result-B candidate); competitor numbers verified, skeptic passed with rewording (019) |

STATE.md line "Result B" should drop "No electron-level model in the literature fits within lattice errors" in favour of "none of the 7 published models checked fits within lattice errors; post-2008 bound-state models not checked".

## Open issues
- Huber-Fischer-Sanchis-Alepuz (2004.00415 and J = 2 follow-up) should be compared on AT2020 ratios with their own error bars; it is the strongest unchecked bound-state candidate.
- Isgur-Paton and Jaffe-Johnson numbers verified only against the 2008 review, not the 1976/1985 originals.
- Lattice ratio errors use uncorrelated quadrature (013 convention), which overstates them slightly; this makes every "within 1 sigma" count generous, not stricter.
- KS ratios depend on their sigma choice; a fair comparison would rerun KS at the sigma that minimizes chi^2 and count it as 1 parameter.

## Reproduce
From repo root (1 thread, < 1 GB, ~6 min total):
```
python tools/pin.py .venv/Scripts/python.exe work/20261120-019-skeptic-c12/fit.py > work/20261120-019-skeptic-c12/run.log
python tools/pin.py .venv/Scripts/python.exe work/20261120-019-skeptic-c12/compare.py > work/20261120-019-skeptic-c12/compare.log
python tools/pin.py .venv/Scripts/python.exe work/20261120-019-skeptic-c12/rank.py > work/20261120-019-skeptic-c12/rank.log
python tools/pin.py .venv/Scripts/python.exe work/20261120-019-skeptic-c12/limits.py > work/20261120-019-skeptic-c12/limits.log
```
