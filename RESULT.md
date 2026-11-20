# Result

Nothing here claims the Yang-Mills Millennium problem is solved, refuted, or that any gap is proved in the continuum.

## Question

Can the Yang-Mills mass gap (pure SU(N), 4D, continuum, infinite volume; Clay formulation: an Osterwalder-Schrader/Wightman theory with spectral gap Delta > 0) be explained with the mathematics behind electron energy levels and gaps: bound-state quantization and its lower-bound methods, band theory, interaction-driven gaps (BCS, Mott, Peierls, excitonic), spin-chain spectral-gap methods (AKLT, Knabe, martingale, Haldane), and quantum-chemistry methods?

## Answer (one line)

**Result C, partial obstruction:** every electron-level technique the program examined proves a gap only where Yang-Mills already has one (strong coupling, fixed lattice spacing in finite volume, or a frustration-free model that is not Yang-Mills); each failure is shown by explicit calculation or counterexample for the route as instantiated here, one route (Gribov region / Payne-Weinberger) stays open on an unproved log-concavity hypothesis, and SDP/bootstrap lower-bound hierarchies were not examined.

## What was shown

Full table, with block type (proved / shown for our version only / argued / open), evidence, claim ids and work folders: `work/20261120-025-resultC-draft/claim.md`. Claims C1-C13 and their status: `CLAIMS.md`. Summary:

| # | technique family | where it gives a gap | what blocks it for continuum YM | block type |
|---|---|---|---|---|
| 1 | Bound-state ED + certified bounds (min-max, Feshbach/Schur, strong-coupling PT); KS as tight binding on the group manifold | Certified two-sided gap intervals for the untruncated SU(2) Kogut-Susskind Hamiltonian on small 2+1D lattices at fixed coupling | The lower bound needs an electric floor on the truncation complement that clears E1; at weak coupling this needs a truncation that grows with volume and 1/g (heuristic scaling), so the basis is exponential while asymptotic freedom sends g -> 0. Small open lattices tend to a cutoff value, not a mass. | shown for our version only; scaling argued |
| 2 | Spin-chain gap methods (AKLT, Knabe, Gosset-Mozgunov, martingale), Rokhsar-Kivelson parent Hamiltonian, Bakry-Emery | Frustration-free RK parent of lattice SU(2): gap >= 1/2 - 2(d-1) beta (strong coupling, = Shen-Zhu-Zhu condition); 2D all beta via Piazza-Shen (cited) | KS is not frustration-free with finite-range terms (natural decomposition proved, general conjectured); the RK correction to KS is extensive, so no gap transfers; Bakry-Emery curvature wins only at strong coupling; RK dynamics is z = 2, not Lorentz-invariant YM. | proved (no transfer, strong-coupling only); conjectured / argued (rest) |
| 3 | Configuration-space QM: Gribov region, Payne-Weinberger | Kinetic-only bound lambda_1 >= pi^2/diam(Lambda)^2 on a convex region | Applying it to the real Hamiltonian needs log-concavity of J psi_0^2 (open Moncrief/Mondal hypothesis); magnetic potential is non-convex in the Gribov region at arbitrarily small fields; the kinetic bound vanishes as volume grows (small T^3). | **open** (conditional on the hypothesis) |
| 4 | Band theory: theta as Bloch momentum, band gap | Correct band data (vacuum energy vs theta, susceptibility, flux energies) | Each Bloch fiber is the full interacting theory, so the in-sector gap is untouched (proved, finite-dimensional version in Lean); bandwidth grows with volume while the fiber gap does not. | proved |
| 5 | Interaction-driven gaps: BCS, Mott, Peierls, excitonic, dual superconductor, gap equations, large N | Controlled only for vector large N, 3D compact U(1), strong coupling | No Fermi surface (bosons), the coupling at the gap scale is fixed, Savvidy background unstable; 1/N exists but the planar limit is not a gap-equation saddle; solved-for "masses" are gauge-dependent or non-spectral (explicit counterexamples). | argued, with toy counterexamples |
| 6 | Quantum chemistry: HF, CI, CC/EOM-CC, DMRG/MPS, Temple/Kato | Upper bounds on energies; numerical gap estimates | Wavefunction methods add no lower-bound machinery beyond certified ED + Schur (row 1); Temple needs an a priori bound equal in difficulty to the gap. RDM/SDP and bootstrap hierarchies do give lower bounds and were not examined. | proved (upper-bound structure); **open** (SDP/bootstrap) |
| 7 | Stepping stones: O(3) rotor chain, spin-1 Heisenberg, AKLT | Certified finite-chain gaps; AKLT thermodynamic-limit gap | Same blocks as rows 1-2 reappear; no technique reaches the O(3) continuum gap, whose known value rests on integrability (TBA), not a proof. | shown for our version only; argued |
| 8 | Bound-state glueball models (Result-B search), incl. functional Bethe-Salpeter | None within lattice errors | Our constituent-gluon "glueball atom" (3 fitted parameters) reproduces 0 predicted states within 1 lattice sigma; none of 12 published models (1976-2026) does either. | numerical (Result B not met) |

## Where it holds and where it doesn't

- Holds: finite lattices, fixed lattice spacing, 2+1D SU(2) KS (rows 1, 6), d-dimensional lattice SU(2) for the RK/Bakry-Emery statements (row 2), lattice Coulomb gauge SU(2) and small-volume T^3 (row 3), abstract Hilbert-space statements and toy models (rows 4, 5), 1+1D chains (row 7), 4D SU(3)/SU(inf) and 2+1D continuum-extrapolated lattice data (row 8).
- Does not hold or not shown: no statement here is a theorem about continuum 4D YM. Rows marked "shown for our version only" or "argued" leave other instantiations of the same technique possible. No 3+1D lattice computation, no SU(3) Hamiltonian, no lattice larger than 2x2 plaquettes.

## Numbers

All values regenerate with `python -m src.reproduce` and are listed with script and output in `work/20261120-025-resultC-draft/number_map.csv` (69 numbers; full run 69/69 passed, `work/20261120-027-full-reproduce/`). Lattice inputs from `data/` (verified against the papers, `work/20261120-004-data-verify/`).

| quantity | setting | value | source |
|---|---|---|---|
| m(2++)/m(0++) | 4D SU(3), lattice continuum | 1.437(11) | Athenodorou-Teper 2020, `data/` |
| m(0-+)/m(0++) | 4D SU(3) | 1.549(16) | same |
| m(0++)/sqrt(sigma) | 4D SU(3) / SU(inf) / 2+1D SU(3) | 3.405(21) / 3.072(14) / 4.368(7) | Athenodorou-Teper 2020, 2021, 2016 |
| KS plaquette gap, g^2 = 1 | 2+1D SU(2), single plaquette | 3.2690992 (exact rational interval) | C5, work 006/010 |
| KS ladder gap, g^2 = 1 | open 2-plaquette ladder | 3.0620133 (certified, width < 1e-8) | C5, work 006/024 |
| KS 2x2 torus gap, g^2 = 3 | all flux sectors / zero flux | [2.02927, 2.03424] / [4.12793, 4.12802] | C9, work 012/016/024 |
| KS 2x2 torus, g^2 = 1 | floor vs ground-state bound | d = -7.861 < E0 upper -6.99998: no lower bound possible | C9, work 012 |
| Strong-coupling plaquette gap | single plaquette | 3g^2/2 + 56/(15 g^6) - 20144/(3375 g^14) + ... | C2, work 005; arithmetic in Lean |
| RK parent gap | lattice SU(2), d dims | >= 1/2 - 2(d-1) beta; Hessian -8 sharp in d = 2 | C6, work 007/010 |
| U(1) RK gap rate | 2D free-boundary boxes | inf_L gap <= 1.2014e-3 at beta = 6 | C8, work 011/016 |
| Pendulum-well bands | N = 1e6 wells, V0 = 4 | bandwidth 2.7102 > fiber gap 1.86542 | C4, work 008/009 |
| Large-N gap equation vs exact | 2D XY at beta = 2 | predicts m = 1.0563e-2; Froehlich-Spencer: m = 0 at large beta | C10, work 014/018 |
| Temple with uncertified input | KS plaquette/ladder | 3 of 88 "bounds" false | C13, work 020/023 |
| AKLT gap (Gosset-Mozgunov) | infinite spin-1 chain | >= 0.3011 (dense n = 10); published certified >= 0.34976 (Rai et al. 2026) | C11, work 015/018/024 |
| Glueball-atom fit | 4D SU(3) / SU(inf) / 2+1D SU(3) / SU(inf) | chi^2/prediction 94 / 196 / 765 / 1305; 0 within 1 sigma | C12, work 013/019 |
| Functional Bethe-Salpeter | 4D SU(3) | 0 fitted, 1/9 ratios within 1 sigma, chi^2/n 104 | work 022 |

## Comparison with existing approaches

- Rigorous lattice gaps exist only at strong coupling (Osterwalder-Seiler; Yarotsky for KS with Gauss law; Shen-Zhu-Zhu via Bakry-Emery; Cao-Nissim-Sheffield) and in 2D (Piazza-Shen log-Sobolev at all beta). Every positive statement in rows 1-2 lands inside those regimes.
- The main known routes toward the continuum (Balaban's constructive RG, continuum stochastic quantization) are not electron-level methods and were not examined.
- Glueball models (bag, flux tube, constituent gluon, holographic, functional Bethe-Salpeter) all sit at the 5-10% level against ~1% lattice errors; the program's own model does not improve on them (row 8).
- The one published family that gives rigorous lower bounds without frustration-freeness or a truncation floor is SDP/bootstrap (variational 2-RDM, lattice-YM positivity bootstrap, gap bootstraps for the TFIM and qubit Hamiltonians). None bounds a YM gap yet; it is the leading open direction.

## Literature check

The obstruction is mostly a synthesis of known facts plus explicit calculations; no novelty is claimed for it as a whole. Searches 2020-2026 (work 001, 007, 011, 022, 029) found no published result contradicting any row. Prior art for each row is listed in `claim.md` ("Literature check") and `literature/01`-`09`. What is ours: the certified two-sided KS gap intervals on small lattices (unrefereed preprints arXiv:2611.11534 and arXiv:2503.15539 overlap in method; numerical overlap unconfirmed), the explicit counterexamples and toy calculations in rows 2-6, the C8 U(1) and SU(2) bounds, the O(3)-chain calibration, and the consolidated table. Claimed full proofs (arXiv:2506.00284, arXiv:2606.19362) were neither relied on nor refuted.

## Known weaknesses

- Certification is floating point (rechecked by LDL^T inertia counts with measured residuals), not interval arithmetic; only the single plaquette is in exact arithmetic. Some rotor-chain rows (L >= 4) and the AKLT n = 12 input still assume Lanczos completeness.
- The truncation-scaling statements (j_max ~ sqrt(N_p)/g, l_max ~ sqrt(N/g)) are heuristic.
- All Hamiltonian computations are 2+1D SU(2) on at most 4 plaquettes; the 2x2 torus is dominated by torelons.
- Several citations were checked against abstracts or bibliographic records only; Maehly 1952 does not resolve; the Gosset-Mozgunov prefactor and the Hamer-Kogut-Susskind coefficients were not checked against the papers.
- Lean covers finite-dimensional lemmas only: C4 abstract steps, the C5 Schur lemma (E0 form and k-eigenvalue count), the easy half of Courant-Fischer, the C2 arithmetic. Physics inputs and analytic statements are not formalized.
- Seven numbers (2x2 torus j_max = 2 builds, the long rotor certification) are checked against committed outputs, not regenerated by the reproduce script.
- Skeptics: every claim C1-C13 passed a skeptic with rewording; the consolidated claim passed a technical skeptic (work 028) and a literature skeptic (work 029), both "survives with required changes", all applied (work 030).

## Reproduce

From the repo root (Python env in `.venv`, pinned in `requirements.txt`; Lean via `elan`):

    python tools/pin.py .venv/Scripts/python.exe -m src.reproduce --quick   # ~1 min, 29 numbers
    python tools/pin.py .venv/Scripts/python.exe -m src.reproduce           # ~35-75 min at 1 thread, all 69 numbers
    cd lean && lake build                                                     # Lean lemmas, no sorry

The script runs on copies under `artifacts/`, checks every row of `work/20261120-025-resultC-draft/number_map.csv`, exits nonzero on any failure, and takes a lock so two runs cannot overlap.
