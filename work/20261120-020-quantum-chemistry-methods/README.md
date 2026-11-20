# 20261120-020 Quantum-chemistry methods (HF, CI, CC, DMRG/MPS/PEPS) as gap techniques

## Question

The last electron-level technique family in `OVERSEER_PROMPT.md` not yet assessed for Result C: Hartree-Fock, configuration interaction, coupled cluster, DMRG / matrix product states / PEPS. Can any of them give a lower bound on the Yang-Mills gap (proof technique), or what quantitative evidence do they give?

## Method

1. For each method, state what it rigorously provides (cited theorems, or proved here).
2. Literature: prior uses on lattice gauge theory, IDs checked against the arXiv API ([arXiv-API]) or Crossref ([Crossref]) on 2026-11-20; [memory] = not checked.
3. Explicit calculation (`qc_methods.py`, imports `src/ks_su2.py` unchanged) on the SU(2) Kogut-Susskind single plaquette and open 2-plaquette ladder at g^2 = 0.5, 1, 2, 3: each method's gap estimate vs the converged ED gap of the same lattice (plaquette j_max = 20, ladder j_max = 6; these equal the C5 certified values to the printed digits). Plus 1D spin-chain checks for the area-law direction.

Setting of every KS number: 2+1D Hamiltonian lattice, SU(2), a = 1, finite lattice, fixed g^2, H = (g^2/2) sum_l j_l(j_l+1) - (2/g^2) sum_p Tr U_p. Nothing here is a continuum or infinite-volume statement.

Rerun: `python tools/pin.py .venv/Scripts/python.exe work/20261120-020-quantum-chemistry-methods/qc_methods.py` (1 thread, a few minutes, well under 2 GB). Output: `ci_truncation.csv`, `hf_ladder.csv`, `cc_plaquette.csv`, `mps_ladder.csv`, `temple.csv`, `mps_chains.csv`, `summary.txt`, `run.log`.

## (a) What each method rigorously provides

Let E_0 <= E_1 <= ... be the eigenvalues of H on the sector, Delta = E_1 - E_0.

- **Rayleigh-Ritz / CI (cited: min-max, Reed-Simon IV Thm XIII.1; Hylleraas-Undheim, Z. Phys. 65, 759 (1930), doi:10.1007/BF01397263 [Crossref]; MacDonald, Phys. Rev. 43, 830 (1933), doi:10.1103/PhysRev.43.830 [Crossref]).** For any trial subspace V, the Ritz values satisfy E_k(V) >= E_k for every k. So Delta(V) = E_1(V) - E_0(V) = Delta + (E_1(V) - E_1) - (E_0(V) - E_0), a difference of two non-negative errors: **no bound in either direction** unless E_0 is known exactly (then Delta(V) + (E_0(V) - E_0) >= Delta, an upper bound). Proved (one line from HUM).
- **Hartree-Fock (proved here for gauge theory).** HF is Ritz over a nonlinear set (product states), so E_HF >= E_0 (upper bound). Koopmans-type or excited-determinant gaps have no bound in either direction. Gauge-specific fact: a gauge-invariant product over **links** of SU(2) wavefunctions is constant (proved: invariance of phi_1(g U_1) phi_4(U_4 g^-1) for all g forces a one-dimensional character of SU(2), which is trivial). So link-HF = the strong-coupling vacuum. The usable HF ansatz is a product over **plaquette holonomies** (class functions chi(U_p)); on the open ladder the cross term from the shared link has zero expectation in such states (the expectation of a generator vanishes in a conjugation-invariant state), so E_HF = 2 E_0(plaquette) exactly (proved).
- **Coupled cluster (cited: standard; e.g. Bartlett-Musial, Rev. Mod. Phys. 79, 291 (2007) [memory]).** Projected equations, not a Rayleigh quotient: **no variational bound**, E_CC can lie below E_0. EOM-CC excitation energies are eigenvalues of a truncated non-Hermitian similarity transform: no bound.
- **DMRG / MPS (cited).** A normalized MPS is a trial state: E_MPS >= E_0. Excited-state DMRG (state-averaged or orthogonalized) gives Ritz values on a subspace: each >= E_k, gap unbounded both ways (same HUM argument). Hastings, J. Stat. Mech. P08024 (2007), arXiv:0705.2024 [arXiv-API]; Arad-Kitaev-Landau-Vazirani arXiv:1301.1162 [arXiv-API]: a 1D local Hamiltonian with gap Delta has area-law ground state and an MPS approximation with bond dimension poly(L) (exp(O(1/Delta)) prefactors). Direction: gap => efficient MPS. Verstraete-Cirac, PRB 73, 094423 (2006), arXiv:cond-mat/0505140 [arXiv-API]: MPS approximate ground states efficiently also for some critical chains. So MPS efficiency does not imply a gap (counterexamples below).
- **PEPS in 2D (cited).** No 2D analogue of Hastings's theorem is proved. PEPS with bond dimension 2 can be critical: Verstraete-Wolf-Perez-Garcia-Cirac, PRL 96, 220601 (2006), arXiv:quant-ph/0601075 [arXiv-API]. Spectral gap of 2D translation-invariant Hamiltonians is undecidable in general: Cubitt-Perez-Garcia-Wolf, Forum Math. Pi 10, e14 (2022), arXiv:1502.04573 [arXiv-API] (general class; says nothing specific about YM).
- **Parent Hamiltonians (cited, [memory]: Fannes-Nachtergaele-Werner, CMP 144, 443 (1992)).** Injective MPS have a gapped frustration-free parent Hamiltonian. That is a gap of a different operator; for YM it is the RK/H_FF route already assessed in C6 and C8 (the coherent state exp(-beta S_W/2) is the gauge-theory version).
- **Temple / Kato / Lehmann (cited: Temple, Proc. R. Soc. A 119, 276 (1928), doi:10.1098/rspa.1928.0098 [Crossref]; Kato, J. Phys. Soc. Japan 4, 334 (1949) [memory]; Lehmann-Maehly [memory]).** Two-sided eigenvalue bounds from a trial vector's mean eta and variance sigma^2, but the lower bound E_k >= eta - sigma^2/(beta - eta) requires an a priori beta <= E_{k+1}. See (d).

## (b) Prior uses on lattice gauge theory

| paper | method | setting | gap content |
|---|---|---|---|
| Schuette, Wichmann, Weichmann, arXiv:hep-lat/9805011 [arXiv-API] | coupled cluster (exp of loop functions) | KS Hamiltonian, lattice | ground state / method |
| Wichmann, Schuette, Metsch, Wethkamp, PRD 65, 094511 (2002), arXiv:hep-lat/0112015 [arXiv-API] | CC + excitations | SU(2) KS, 3+1D, infinite lattice | SU(2) glueball masses; no error bound |
| Carlsson, McKellar, PRD 68, 074502 (2003), arXiv:hep-lat/0303016 [arXiv-API] | analytic variational (Ritz) | SU(N <= 5) KS, 2+1D | glueball masses; upper bounds on energies only |
| Byrnes, Sriganesh, Bursill, Hamer, PRD 66, 013002 (2002), arXiv:hep-lat/0202014 [arXiv-API] | DMRG | massive Schwinger, 1+1D | extrapolated gap (numerical) |
| Banuls, Cichy, Jansen, Cirac, JHEP 11 (2013) 158, arXiv:1305.3765 [arXiv-API] | MPS | Schwinger, 1+1D, continuum extrapolation | mass spectrum (numerical) |
| Milsted, PRD 93, 085012 (2016), arXiv:1507.06624 [arXiv-API] | uniform MPS, tangent space | O(2), O(4) rotor chains (abstract: equivalent to KS matter-free non-Abelian LGT) | mass gap and beta function (numerical) |
| Kuhn, Zohar, Cirac, Banuls, JHEP 07 (2015) 130, arXiv:1505.04441 [arXiv-API] | MPS | SU(2), 1+1D | string breaking, no gap |
| Banuls, Cichy, Rep. Prog. Phys. 83, 024401 (2020), arXiv:1910.00257 [arXiv-API] | review | TN for LGT | review |
| Magnifico, Felser, Silvi, Montangero, Nat. Commun. 12, 3600 (2021), arXiv:2011.10658 [arXiv-API] | TTN | 3+1D QED, finite density | phases, no glueball gap |
| Cataldi, Magnifico, Silvi, Montangero, PRR 6, 033057 (2024), arXiv:2307.09396 [arXiv-API] | TTN | 2+1D SU(2) with matter, hardcore-gluon truncation | finite density, not the glueball gap |
| Ebner, Schafer, Seidl, Muller, Yao, PRD 110, 014505 (2024), arXiv:2401.15184 [arXiv-API] | ED, entanglement | 2+1D SU(2) KS plaquette chains | area law ground state, volume law excited states |
| Mathew et al., arXiv:2501.18301 [arXiv-API] | MPS (loop-string-hadron) | 1+1D SU(2) | dynamics toolbox |
| Cataldi, arXiv:2501.11115 [arXiv-API] (thesis) | TN review | Hamiltonian LGT | phases, scarring |
| Ludwig, Jakobs, Urbach, arXiv:2602.23213 [arXiv-API] | TN | 2+1D SU(2) quantum link model, hexagonal | static potential, Luscher term |

Literature finding (two web searches plus the arXiv API, not exhaustive): no 2020-2026 tensor-network paper found that computes the pure-gauge non-Abelian glueball gap in 2+1D or 3+1D with a continuum extrapolation. Every listed use gives upper bounds or extrapolated numbers, never a lower bound. Milsted-Osborne arXiv:1604.01979 (already in `literature/01`) is the TN programme for continuum YM.

## (c) Explicit calculation: the gap estimate lands on either side

Gap error = method gap - exact gap. Full rows in the CSVs; selected values:

**CI hierarchies (`ci_truncation.csv`).** "depth k" = all states within k plaquette moves of the strong-coupling vacuum (CIS, CISD, ... on the HF/strong-coupling reference); on the plaquette this equals the j_max truncation. "sel_E1 n" = selected CI keeping the n configurations of largest weight in the E1 eigenvector.

| lattice, g^2 | depth-CI gap error, k = 1..6 | sel_E1 gap error, n = 2, 4, 8, 12, 24, 32 |
|---|---|---|
| plaquette, 0.5 | +4.37 +2.2 +0.82 +0.21 +0.031 +0.0025 | (n = 2, 3, 4) +0.086, -2.44, -1.55, then + |
| ladder, 0.5 | +2.75 +2.14 +1.42 +0.90 +0.58 +0.41 | +4.74 +2.75 +2.41 +1.14 -2.96 -2.25 |
| ladder, 1.0 | +0.61 +0.37 +0.25 +0.097 +0.012 +8e-4 | -0.81 -2.38 -2.71 -3.02 +0.029 +0.0093 |
| ladder, 2.0 | +0.21 +0.078 +3.4e-3 +5.5e-5 +4.5e-7 +2.1e-9 | -3.35 (n=2), then + |

Numerical observation (these lattices, these couplings): the vacuum-built CI hierarchies overestimate the gap at every k sampled and converge monotonically. This is not a theorem: selected CI, an equally standard choice, underestimates by up to 3.0 (ladder g^2 = 1, n = 12: Delta(V) = 0.045 vs 3.062) and changes sign with n. Each individual Ritz value is still above its eigenvalue (HUM), which is all that is guaranteed.

**HF on the ladder (`hf_ladder.csv`).** E_HF - E_0 = 0.0295, 0.0145, 0.0012, 1.0e-4 at g^2 = 0.5, 1, 2, 3 (>= 0 as proved). Koopmans gap (= plaquette gap) minus exact ladder gap: +0.370, +0.207, +0.072, +0.026. Overestimates here; no bound.

**Coupled cluster, plaquette (`cc_plaquette.csv`).** exp(S)|0>, S = sum_{2j <= n} s_j chi_j (the Schuette-Wichmann ansatz on one plaquette), projected equations solved to residual < 3e-15. E_CC - E_0 alternates in sign with n: g^2 = 0.5: -16.0, +0.98, -1.11, +0.14, -0.093, ...; g^2 = 1: -1.00, +0.080, -0.015, ... So **E_CC < E_0** at n = 1, 3, 5, ... (non-variational, as expected). EOM-CC gap error alternates too: g^2 = 0.5: +4.20, -1.97, -0.57, -0.41, -0.24, +0.29, -0.099, +0.034, -0.012, +0.0016. Convergence is fast at g^2 >= 1 and slow at g^2 = 0.5.

**MPS / DMRG-centre on the ladder (`mps_ladder.csv`).** Cut through the shared link; keep D Schmidt multiplets per side (all electric sectors c), chosen from the E0-targeted, state-averaged or E1-targeted density matrix of the exact states (the best environments a two-site DMRG could reach), then Ritz on the product span. Gap error, averaged target, D = 3, 4, 6, 8, 12, 16, 24: g^2 = 0.5: -1.35, +2.22, +1.60, +0.67, +0.0028, +0.32, +0.020 (non-monotone, both signs); g^2 = 2: +0.27, +0.29, +0.074, -0.0050, +0.0013, ... E1-targeted runs underestimate at most D (g^2 = 1: -0.81 ... -0.0020). E0-targeted overestimate. The sign of the gap error follows the targeting choice, not the physics.

Status: numerical, finite lattice, fixed g^2, 2+1D SU(2) KS.

## (d) What can be made rigorous, and the circularity

**Temple/Kato/Lehmann need the gap they are meant to prove (proved).** A lower bound on Delta needs a lower bound on E_1. Temple for E_1 from trial phi_1 gives E_1 >= eta_1 - sigma_1^2/(beta - eta_1) only if beta <= E_2. Any certified beta is a lower bound on E_2, i.e. a statement of the same kind as the one sought, one level up; using it requires a lower bound on E_3 through Temple again, and so on. The regress closes only if some step gets a lower bound from an independent source. In C5 that source is the electric floor on the truncation complement, d = (g^2/2) c_Q - 4 N_p/g^2 (Feshbach/Schur lemma, `work/20261120-006-certified-gap-bounds/`), and that floor is exactly what fails at weak coupling: c_Q must exceed ~8 N_p/g^4, which needs j_max ~ g^-2 on every link and a dimension exponential in the volume (C5, C9, C11).

Numerical check (`temple.csv`): with beta = the CI value of E_2 (the only number available without an a priori bound; it is >= E_2 by HUM, so on the wrong side), 3 of 88 Temple "bounds" are false: plaquette g^2 = 0.5, k = 2: "E_1 >= -0.364" but E_1 = -1.667; k = 3: "E_1 >= -1.611"; ladder g^2 = 0.5, k = 2: "E_1 >= -6.778" but E_1 = -7.397. With the exact E_2 as beta, 0 of 88 fail, as the theorem says. The false cases are at the weakest coupling and smallest truncation, where the CI E_2 is furthest above the true E_2.

**Area law: gap => MPS, not MPS => gap (proved / numerical).**
- Ferromagnetic Heisenberg ring H = -sum S_i.S_{i+1}: the product state |up...up> (bond dimension 1) is an exact ground state, and the gap is 1 - cos(2 pi/L) -> 0 (proved, standard; checked by ED at L = 8, 12: 0.29289, 0.13397).
- Critical transverse-field Ising chain (h = 1, open): TT-SVD truncation to D = 8 has relative energy error 1.8e-12, 5.5e-10, 4.3e-9 at L = 8, 12, 16, while the finite-L gap 0.369, 0.251, 0.190 falls ~1/L (gapless in infinite volume, exact solution, cited). Gapped h = 2 for comparison: error 1e-14 at D = 8, gap 2.06. Small bond dimension at finite size occurs in both.
- Consequence: a converged DMRG/PEPS calculation at bond dimension D is evidence about the state, not a certificate of a gap. In 2D (the YM case) not even the forward direction (gap => efficient PEPS) is proved.

## Per-method result table

| method | rigorous output | gap lower bound? | KS calculation here | status |
|---|---|---|---|---|
| HF (plaquette-product) | E_HF >= E_0 | no | E_HF = 2 E_0(plaq) proved; Koopmans gap +0.03 to +0.37 off | proved / numerical |
| CI (any truncation) | E_k(V) >= E_k each k (HUM) | no; gap error has either sign | depth-CI over, selected-CI under by up to 3.0 | cited / numerical |
| Coupled cluster, EOM-CC | none (non-variational) | no | E_CC below E_0 at odd n; EOM gap error alternates | numerical |
| DMRG / MPS | E_MPS >= E_0; Ritz bounds per level | no; gap error has either sign | sign follows targeting, non-monotone in D | cited / numerical |
| Area law (Hastings, AKLV) | 1D gap => efficient MPS | converse false | ferro Heisenberg D = 1 gapless; critical TFIM D = 8 | cited / proved / numerical |
| PEPS (2D) | E_PEPS >= E_0 | no; critical PEPS exist | none run | cited |
| Temple / Kato / Lehmann | two-sided bounds given beta <= E_{k+1} | only with an a priori bound on E_{k+1}: circular | 3/88 false bounds with uncertified beta | proved (logic) / numerical |
| Parent Hamiltonian of MPS/PEPS | gap of a different FF operator | reduces to C6, C8 | not repeated | cited |

## Verdict

Quantum-chemistry and tensor-network methods give **upper bounds on individual energies and numerical gap estimates**, never a gap lower bound. Their gap estimate is a difference of two upper bounds and lands on either side of the true gap; on the SU(2) KS plaquette and ladder (2+1D, finite lattice, fixed g^2) this was shown explicitly for selected CI, CC/EOM-CC and MPS truncations. The only rigorous upgrade (Temple/Kato/Lehmann) needs an a priori lower bound on the next eigenvalue, which is circular unless closed by an independent floor; the floor used in C5 is the one that fails at weak coupling. The area law runs from gap to MPS only, and has no proved 2D form. Prior art (CCM glueballs 1998-2002, variational SU(N) 2003, DMRG/MPS Schwinger and rotor chains, TN 2+1D SU(2) at finite density) is consistent with this: all numerical, none bounded below. For Result C: this technique family is covered, and adds no lower-bound mechanism beyond C5; it inherits C5's weak-coupling floor obstruction. Novelty: none (the HUM, non-variational CC and area-law-direction facts are textbook); the contribution is the explicit KS demonstration.

## Open issues

- No MPS/DMRG run on a longer plaquette chain; the ladder "DMRG centre" uses the exact states to pick environments, so it is the best case of the method, not a DMRG simulation.
- Kato, Lehmann-Maehly, Fannes-Nachtergaele-Werner and Bartlett-Musial references cited from memory, not checked.
- The 2020-2026 literature sweep for TN glueball spectra was two searches plus API checks; a dedicated search could still find a 2+1D SU(2) TN glueball paper (would be prior art for numerics only, not for a bound).
- Lehmann-Maehly with a Goerisch-type certified beta (from a solvable comparison operator H >= H_0) is the standard way to break the regress in quantum chemistry. For KS the natural comparison operator is the electric term plus -||W||, which is the C5 floor; no better comparison operator is known at weak coupling.
