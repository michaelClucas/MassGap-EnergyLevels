# 013: Result-B attempt, constituent-gluon bound states ("glueball atom")

## Question
Can an atomic-bound-state model (constituents in a confining potential, fine/hyperfine structure from one-gluon exchange) reproduce lattice glueball spectra with fewer fitted parameters than reproduced quantities, within lattice errors, and beat or match existing models? (Result B in `OVERSEER_PROMPT.md`.)

## Novelty, stated first
This model is a reimplementation of the spin-formalism constituent-gluon potential model: Cornwall-Soni type massive gluons, adjoint Cornell potential with Casimir scaling (as in Kaidalov-Simonov hep-ph/9911291 and Brau-Semay / Mathieu-Buisseret-Semay arXiv:0802.0088 "Model A"), Breit-Fermi spin terms with Gaussian (Godfrey-Isgur-style) smearing. Nothing here is new. The only question was whether it fits the 2020-2021 lattice data within errors. It does not.

## Mechanism and electron-analogy map
| atomic physics | here |
|---|---|
| electron | constituent gluon: massive spin-1 quasi-particle, mass m (4D), spinless single-polarization quasi-particle (2+1D) |
| nucleus + Coulomb potential | the other gluon(s); potential V = sigma_A r - kappa erf(r/r_s)/r (4D), sigma_A r + (g^2N/2pi) ln r (2+1D) |
| lattice | none (continuum two- and three-body problem) |
| fine / hyperfine structure | Breit-Fermi spin-spin, spin-orbit, tensor from V_V = OGE, V_S = sigma_A r, first-order perturbation theory |
| level spacing / gap | glueball mass = lowest bound-state energy of two gluons in the 0++ channel |

Setting: continuum, infinite volume, 3+1D SU(3) and SU(inf), 2+1D SU(3) and SU(inf). Units sqrt(sigma) = 1 (sigma = fundamental string tension), so sigma is not a parameter.

Hamiltonians (`model.py`):
- Two gluons: H = 2 sqrt(p^2 + m^2) (or 2m + p^2/m, "nr") + V + C0. Colour singlet is symmetric, so Bose symmetry forces L + S even; P = (-1)^L, C = +.
- Three gluons (d-coupled, C = -): hyperspherical K = 0 approximation, pairwise Casimir potentials (sigma_A/2, kappa/2 per pair), kinetic 3 sqrt(P^2/3 + m^2), spin-symmetric S = 1, 3, giving 1--, 3-- (4D) and 0-- (2+1D).
- Radial equations on a finite-difference grid; the relativistic operator is built from the eigendecomposition of p^2. Grid check (`conv.log`): 4D levels move < 4e-4; 2+1D 0++ moves 0.007 (about 1 lattice sigma, from the L = 0 2D centrifugal -1/(4r^2)).

## Parameter count (all choices counted)
Continuous, fitted: m/sqrt(sigma), kappa = N alpha_s, C0/sqrt(sigma) (4D k = 3); m, C0 (2+1D; kappa fixed by the lattice value of sqrt(sigma)/(g^2 N) from AT2016 Table B1, an input not a fit).
Discrete choices (not fitted but each one shapes the answer): (1) Casimir scaling sigma_A = (C_A/C_F) sigma; (2) Breit-Fermi coefficients borrowed from spin-1/2 quarks (not derived for spin-1 constituents); (3) smearing width r_s = 1/m; (4) first-order PT, no L-mixing by the tensor force; (5) spin formalism, not helicity formalism; (6) d-coupled S-wave three-gluon states only; (7) constant per gluon (3-gluon states get 1.5 C0); (8) level assignment = energy order within each J^PC; (9) large N at fixed lambda = g^2N with m, C0 fixed in sqrt(sigma) units; (10) 2+1D gluons spinless with intrinsic parity -. Effective count is 3 fitted + 10 discrete choices; I quote k = 3 (4D) and k = 2 (2+1D) for chi^2 purposes, which flatters the model.

## Fits (`fit.py`, output `run.log`, `fits.csv`)
Lattice: AT2020 Table 17 (4D SU(3)), AT2021 Table 38 (4D SU(inf)), AT2016 Tables B4/B10-B12 (2+1D). All read from `data/`.
Model-error budget: 5% of each mass added in quadrature (a convention; nothing derives it).

| variant | fitted on | parameters | result |
|---|---|---|---|
| rel, k = 2 (m, kappa; C0 = 0) | 0++, 2++ | m runs to the 0.3 bound, kappa = 1.26 | no physical solution: `scan_k2.log` shows that with C0 = 0, fitting 0++ forces m(2++) - m(0++) >= 1.57 for every m >= 0.4 (lattice 1.49), and at m = 0.3 the 1/m^2 spin terms give negative masses (2-+ = -4.7). Discarded. |
| nr, k = 2 | 0++, 2++ | m = 0.95, kappa = 1.42 | cannot even fit the 2 inputs (chi^2 = 1.3e3) |
| **rel, k = 3** (best) | 0++, 2++, 0-+ | m = 3.13, kappa = 0.72 (alpha_s = 0.24), C0 = -3.30 | fits the three inputs; predictions below |
| nr, k = 3 | 0++, 2++, 0-+ | m = 3.44, kappa = 0.89, C0 = -3.68 | similar, slightly worse |
| 2+1D rel, k = 2 | 0++, 2++ | m = 2.46, C0 = -2.35 | predictions below |

The working 4D fit needs m ~ 3.1 sqrt(sigma) (~1.5 GeV) and a constant -3.3 sqrt(sigma) cancelling most of 2m. That is a heavy-constituent model with a large ad hoc self-energy, not a light-gluon atom.

### 4D SU(3), rel k = 3, predictions (m/sqrt(sigma))
| state | lattice | model | pull (lat) | pull (5%) |
|---|---|---|---|---|
| 0++* | 5.855(41) | 5.984 | +3.1 | +0.44 |
| 2-+ | 6.32(9) | 5.924 | -4.4 | -1.21 |
| 2++* | 6.788(40) | 6.947 | +4.0 | +0.46 |
| 0-+* | 7.29(13) | 7.092 | -1.5 | -0.51 |
| 4++ | 7.60(12) | 7.084 | -4.3 | -1.30 |
| 3++ | 7.71(9) | 7.168 | -6.0 | -1.37 |
| 2-+* | 8.18(8) | 7.602 | -7.2 | -1.39 |
| 1-- | 8.31(10) | 6.513 | -18 | -4.2 |
| 1-+ | 8.48(12) | 5.835 | -22 | -6.0 |
| 1+-, 3+-, 2+-, 2--, 4+-, 1+-* | 6.07 ... 9.02 | no state | - | - |

0 of 11 predictions within 1 lattice sigma; chi^2/prediction = 94 (lattice errors), 6.1 (5% budget). Six lattice states, including the third-lightest glueball 1+-, have no counterpart. Structural failures: (a) spurious low 1-+ (L = 1, S = 1) at 5.8 vs lattice 8.48, the known defect of the spin formalism that motivated the helicity formalism of 0802.0088; (b) lightest C = - state is 1-- (S-wave three gluons), while lattice has 1+- at 6.07 and 1-- at 8.31.

### N-dependence (4D SU(inf), same parameters, lambda fixed)
0 of 12 within 1 sigma, chi^2/pred = 196 (lattice), 5.2 (5%). The model's only N-dependence is sigma_A/sigma: 9/4 -> 2, which lowers m(0++) by 3% (3.40 -> 3.30); the lattice drops 10% (3.405 -> 3.072, 16 sigma miss). Letting kappa or m run with N would add parameters.

### 2+1D SU(3), rel k = 2
0 of 13 within 1 sigma; chi^2/pred = 765 (lattice), 4.4 (5%). 0++* = 7.00 vs 6.486(13) (40 sigma). Closest: 3-gluon 0-- = 6.485 vs 6.391(14) (6.7 sigma, 1.5%), with no extra parameter. 2-+ = 2++ is parity doubling, built into any model (not counted as a success). The model has no 0-+, 0+-, 1^PC states. SU(inf) with same m, C0: 0/15 within 1 sigma, chi^2/pred 1305 (lattice), 3.0 (5%).

## Comparison with existing models (`compare.py`, `comparison.csv`, `competitors.csv`)
Ratios m/m(0++) against AT2020 (removes each model's scale; ratio errors uncorrelated-quadrature). Published numbers from: MSB 0802.0088 Table 3 (Model B); Kaidalov-Simonov hep-ph/9911291 Table 6; SSJC PRL 76 (1996) as quoted in 0802.0088 Table 3; Brower-Mathur-Tan hep-th/0003115 Tables 1-2; Isgur-Paton flux tube via review 0810.4453 Table 5; Jaffe-Johnson bag via 0810.4453 text; hard wall via 0810.4453 Table 8. Numbers were transcribed by a lookup subagent from ar5iv HTML; the MSB Table 3 extraction was garbled and cross-checked against its Table 4 ratios and 0810.4453 Table 4. Not independently re-verified here.

| model | ratio params fitted to glueballs | states compared | within 1 sigma | chi^2/n lat | chi^2/n 5% |
|---|---|---|---|---|---|
| Kaidalov-Simonov | 0 (sigma from Regge, alpha_s = 0.3 chosen) | 10 | 3 (0++*, 3++, 2--) | 98 | 3.3 |
| MSB helicity B | 1 constant on 0++, 0-+ (+ alpha_s from earlier fit) | 8 | 0 | 84 | 3.6 |
| this work rel k = 3 | 3 on 0++, 2++, 0-+ | 11 (2 fitted) | 2 (both fitted), 0 predicted | 70 | 5.4 |
| SSJC Coulomb gauge | 0 (fitted to static potential) | 8 | 0 | 304 | 13.7 |
| flux tube (IP) | hand-chosen | 5 | 0 | 447 | 12.1 |
| bag (JJ) | 0 for ratios | 4 | 0 | 580 | 19.5 |
| hard wall (Dirichlet) | 0 for ratios | 2 | 0 | 20 | 0.59 |
| BMT AdS7 | 0 | 5 | 0 | 1705 | 62 |

Our lower lattice-only chi^2/n than KS and MSB comes from including fitted states; on predictions alone it is 94. With 3 fitted parameters against 0 or 1 for KS and MSB, it is worse than both. No model in the table reproduces the lattice spectrum within lattice errors; the 2020 errors (~0.5-1.5%) are far below any of these models' accuracy.

## Verdict
Result B: **not met.** With k = 3 fitted parameters (plus ~10 discrete choices), the model reproduces 0 predicted quantities within lattice errors in 4D SU(3) (0/11), 4D SU(inf) (0/12), 2+1D SU(3) (0/13) and 2+1D SU(inf) (0/15). It needs number reproduced > k; it gets 0. It misses ordering-level features (spurious 1-+, 1-- below 1+-, no 1+- at all), and its N-dependence at fixed lambda is a third of the lattice's. It does not beat Kaidalov-Simonov (zero glueball-fitted parameters, 3/10 within 1 sigma) or MSB helicity. It is a reimplementation of a known model class, so even a good fit would not have been new. With a 5% model-error budget all of these models, ours included, are "roughly right" (chi^2/n ~ 3-6); that level of agreement is what the constituent picture gives, and it was known by 2000.

What the attempt does show: the atomic analogy needs (i) a heavy constituent (m ~ 3 sqrt(sigma)) plus a large negative constant to get 2++/0++ = 1.44, which undercuts the "light gluon atom" picture, and (ii) the gluon spin cannot be treated as an electron-like spin (spurious states), which is the helicity-formalism point of 0802.0088.

## Open issues
- Helicity formalism (MSB) instead of spin formalism would remove 1-+; it is already published and also misses lattice errors (0/8).
- Three-gluon states are crude (K = 0 hyperspherical, symmetric relativistic kinetic, no P-wave), so the 1+- channel is untested rather than refuted; MSS arXiv:0803.0815 Table 4 (full three-body) also puts 1-- (3999 MeV) below 1+- (4626 MeV), the same ordering failure.
- Tensor-induced L mixing and second-order spin effects ignored; spin terms are large (0.5-1 sqrt(sigma)), so PT is not controlled.
- Competitor numbers come from a subagent lookup (ar5iv HTML); a skeptic should verify MSB Table 3 and KS Table 6 before any reuse. Kaidalov-Simonov print sigma_f = 0.238 (caption) vs 0.228 (text); ratios unaffected.
- 2+1D L = 0 grid convergence is ~1 lattice sigma for 0++ (not decisive given 40 sigma misses).

## Reproduce
From repo root (1 thread each, < 1 GB RAM, ~10 min total):
```
python tools/pin.py .venv/Scripts/python.exe work/20261120-013-resultB-constituent-gluon/fit.py > work/20261120-013-resultB-constituent-gluon/run.log
python tools/pin.py .venv/Scripts/python.exe work/20261120-013-resultB-constituent-gluon/scan_k2.py
python tools/pin.py .venv/Scripts/python.exe work/20261120-013-resultB-constituent-gluon/compare.py
python tools/pin.py .venv/Scripts/python.exe work/20261120-013-resultB-constituent-gluon/conv.py 3.1323 0.7187 -3.3044 2.4635 -2.3482
```
