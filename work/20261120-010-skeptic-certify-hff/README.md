# 20261120-010: Skeptic audit of C5 (certified KS gap bounds) and C6 (frustration-free H_FF)

## Question

Can claims C5 (`work/20261120-006-certified-gap-bounds/`, `src/certify.py`) and C6 (`work/20261120-007-frustration-free-gauge/`) be broken? The auditor did not do either piece of work.

## Method

- Line-by-line reading of both READMEs and of `src/certify.py`, `certify_gaps.py`, `torus_j32_lower.py`.
- `plaq_exact.py` (new, does not import `src/`): single plaquette from scratch (Peter-Weyl basis chi_j, H = 2g^2 j(j+1) - (2/g^2)[adjacent j]), the same Feshbach lemma redone in **exact rational arithmetic** (Sylvester inertia of tridiagonal LDL^T with `fractions.Fraction`, bisection to 1e-13), 2 j_max = 40, g^2 in {0.2, ..., 10}. Then every one of the 120 plaquette intervals in 006 `results_plaq1.csv` is checked against the exact intervals. Log `plaq_exact.log`, output `plaq_exact.csv`, runtime ~2 s.
- `rerun_certify.py`: reruns 25 points through the 006 code path (plaquette j_max 2-8; ladder j_max 1, 2, 4; 2x2 torus j_max 1/2 and 1 only, from a basis built to j = 1, so no j_max >= 3/2 torus job), checks that all 006 rows for each (lattice, g^2) pairwise intersect, and measures the slack between the 006 floating-point plaquette bounds and the exact ones. Log `rerun_certify.log`, ~3 s.
- `hff_checks.py`: (1) single-plaquette H_FF in the gauge-invariant sector must have E0 = 0 exactly for every beta (this tests the coefficients of the explicit form together); (2) finite-difference Hessian of the Wilson action in normal coordinates on periodic lattices d = 2 (L = 2, 3), d = 3 (L = 2), random and pi-flux (all U_p = -1) configurations, compared with -8(d-1). Log `hff_checks.log`, ~1 min.

All runs: `python tools/pin.py .venv/Scripts/python.exe <script>`, 1 BLAS thread, < 200 MB.

Settings: C5 checks are 2+1D Hamiltonian lattice, SU(2), a = 1, finite lattice, fixed g^2. C6 checks are SU(2), d-dimensional spatial lattice, finite volume.

## C5: verdict **survives**, with wording changes

Lemma (Feshbach/Schur) checked step by step:
- Hypotheses: D = QHQ >= d on Q, X >= BB*, E < d. The step 2Re<u,Bv> >= -2||B*u|| ||v|| followed by minimizing over ||v|| gives -||B*u||^2/(d-E). Correct. Codimension: the nonnegative subspace V (+) Q has codimension <= k, so min-max gives lambda_k >= E. Multiplicity is handled by the codimension count. **proved (rechecked)**.
- Choice E = U (upper bound), L = f(U): f is nonincreasing and f(U) <= lambda_k(A) <= U, so f(L) >= f(U) >= L. Correct. The code adds the rounding margin to U before use, which keeps this valid when LAPACK slightly underestimates lambda_k(A). In the nested floor, `e1` has no margin, but a sign error there would need a rounding error larger than the margin subtracted from L, so it does not matter.
- Gap: `Delta_lo = lo[1] - hi[0]` (lower bound on E1 minus **upper** bound on E0), `Delta_hi = hi[1] - lo[0]`. Correct.
- Floor d: the complement needs only ||QWQ||, and ||QWQ|| <= ||W|| <= sum_p ||Tr U_p|| = 2 N_p (multiplication by a function with |Tr U| <= 2). So d = (g^2/2) c_Q - 4N_p/g^2 is valid. It is wasteful, since ||W|| = 2N_p is attained only at U = 1, but it is not wrong. c_Q is minimized over all flux sectors, which is conservative. The enumeration stops correctly (the analytic bound J(J+1) + J^2/k + J increases with J).
- Shell coupling: Tr U_p moves each loop spin by +-1/2, and on all three lattices the four links of a plaquette are distinct (checked on 2x2), so QWP lands in the j_max + 1/2 shell. The exact BB* holds. Crude shell: Q Tr U_p psi = Q Tr U_p Pi_p psi, then Cauchy-Schwarz over plaquettes. Correct.
- Sector: plaquette and ladder intervals are for all of G. Torus intervals are for the zero-flux sector only, as already stated in C1 and C5. Torelons in nonzero-flux sectors lie lower, so the torus interval is not the lowest excitation of the full gauge-invariant space.

Numerical evidence:
- Exact-rational plaquette bounds, e.g. g^2 = 1: Delta in [3.269099185521, 3.269099185521] (width < 1e-12); g^2 = 0.2: 3.871581222302. **All 120 plaquette intervals in 006 contain the exact intervals; 0 violations.**
- The float bounds beat the exact bounds by at least **2.5e-13** (smallest slack, g^2 = 10, j_max = 1/2, E0 lower). So at the quoted widths (1e-11 to 1e-12) the `10 n eps ||A||_inf` margin carries real weight. It held on the plaquette, but only by about 2e-13.
- Rerun: plaquette points reproduce 006 to <= 1e-12. Ladder and torus reruns differ from 006 only because the basis top (j_top) differs, which changes the nested floor and the exact/crude shell choice. Every rerun interval contains the tightest 006 interval (e.g. torus g^2 = 10, crude shell j_max = 1: [14.86050, 14.86055] contains [14.860549, 14.860553]). For each lattice and g^2, all 006 rows pairwise intersect (0 inconsistencies).

Weak points:
1. **The torus j_max = 3/2 interval [4.1257, 4.1281] at g^2 = 3 has no record in the repo.** It is not in `results.csv` or `run.log`, and `results_periodic2x2_j32.csv` does not exist. The 006 README says it came from a development run whose rerun died. It is quoted in CLAIMS. I did not rerun it (excluded by budget).
2. Rounding: the margin `10 n eps ||A||_inf` is a heuristic constant, not a proved LAPACK backward-error bound. It does not cover errors in the floating-point matrix elements (CG/intertwiner products on ladder and torus). Those errors are ~1e-16 relative, so they are almost certainly harmless, but "proved modulo rounding" is the honest status for ladder and torus. Plaquette is now proved outright (exact arithmetic).
3. The large-n inertia route (eigsh plus LDL with Frobenius residual) is used only in the unrecorded dev run.

Required CLAIMS.md wording for C5:
- Status: "proved for the single plaquette (exact rational arithmetic, skeptic 010); ladder and torus: computer-assisted, floating point with a heuristic rounding margin (not interval arithmetic); skeptic passed with rewording".
- Remove the torus g^2 = 3, j_max = 3/2 interval [4.1257, 4.1281], or mark it "unrecorded dev run, not reproduced". The recorded torus g^2 = 3 interval is [4.08345, 4.13804] (j_max = 1).

## C6: verdict **survives (math correct), weakened in scope and novelty wording**

Checked:
- Q_l = (1/2) psi grad_l psi^{-1}, so Q_l psi = 0. Q_l^*Q_l = E_l^2 + (beta^2/16)|grad_l S|^2 - (beta/8) Delta_l S, with E^2 = -Delta_{S^3}/4 on unit S^3 (spin j <-> degree 2j harmonics, n(n+2)/4 = j(j+1)). Rederived: **correct**. Frustration-free, gauge invariant (S and the bi-invariant operators commute with G^V), and unique ground state on the full space (connected compact manifold, psi > 0 smooth): **correct**.
- Explicit form: Delta(linear function) = -3 x linear on S^3, and |grad(u.a)|^2 = |a|^2 - (u.a)^2. Numerically, single-plaquette H_FF = 4E^2 - (3beta/2)cos + (beta^2/4) sin^2 has E0 = 0 to 1e-14 for beta in [0, 10]. **Correct**.
- Identity: with the repo's H_KS = (g^2/2) sum j(j+1) - (2/g^2) sum Tr U_p (no constant) and beta = 16/(3g^4), the constants cancel exactly, giving H_FF = H_KS/(g^2/2) + (beta^2/16) sum_l |grad_l S|^2. The CLAIMS line is right as written. The 007 README's "+ const" is only needed for its (1 - 1/2 tr) convention. The correction is >= 0, so H_FF >= H_KS/(g^2/2). It is bounded per link, ||(beta^2/16)|grad_l S|^2|| <= (beta^2/4)(d-1)^2 since |grad_l (1/2 tr U_p)| <= 1, but it is **extensive**. The identity therefore transfers nothing about gaps by itself. Its only content is that the ground states agree to first order in 1/g^4, which was checked. Fine as stated, provided nobody reads it as a gap comparison.
- Langevin correspondence: <f psi, H_FF f psi> = (1/4) Integral |grad f|^2 dpi_beta, so gap = (1/4) x the Poincare constant. **Correct**.
- Bakry-Emery: unit S^3 has Ric = 2. Check: Lichnerowicz gives lambda_1(-Delta) >= 3, which matches the spin-1/2 Casimir 3/4 x 4. The curves U e^{sX} with ||X||_op = 1 are unit-speed geodesics, so the block bound |1/2 tr(A X B Y)| <= 1 gives |Hess(1/2 tr U_p)(v,v)| <= (sum |v_l|)^2 <= 4 sum |v_l|^2. Each link lies in 2(d-1) plaquettes (verified on the generated lattices), so Hess S >= -8(d-1) and gap >= (1/4)(2 - 8(d-1)beta) = 1/2 - 2(d-1)beta. **Correct.** Numerics: lambda_min(Hess S) = -8.000 at pi-flux in d = 2, so the constant is **sharp in d = 2**. In d = 3, pi-flux gives -12, not -16, so the constant may be loose there (not proved).
- In SZZ normalization (action N beta_SZZ Re tr U_p, metric tr(XY^*)), beta_ours = 4 beta_SZZ for N = 2, and beta_ours < 1/(4(d-1)) becomes beta_SZZ < 1/(16(d-1)), which is the SZZ condition recorded in `literature/02-rigorous-lattice.md`. So the bound is **exactly** SZZ's Bakry-Emery condition. The convention mapping is derived here from the recorded condition, not read from the SZZ paper.
- "Uniform in volume": justified for every finite torus or box, because the constant uses only local data. The **infinite-volume** statement (GNS Hamiltonian, limit of ground states) is not written out (007 open issue). The CLAIMS setting "infinite volume limit uniform" overstates this.

Novelty:
- RK construction known, and the bound identical to SZZ: fair. Yarotsky for strong-coupling KS: fair (not re-audited here).
- **New problem for the program:** the repo's own `literature/02-rigorous-lattice.md` lists Piazza-Shen, arXiv:2610.35286. Its arXiv abstract (checked) states a log-Sobolev inequality for the 2D SU(N) Wilson lattice measure at every fixed beta, with a constant uniform in volume. If the Dirichlet form is the standard gradient form on SU(N)^edges (not confirmed from the abstract; boundary conditions also unstated), then LSI implies Poincare. Via 007's fact 6 that gives gap(H_FF) > 0 uniformly in volume for **all** beta in d = 2. That is the qualitative content of 007's "Target T". 007 `literature-notes.md` says no Langevin gap result beyond the SZZ regime was found, which contradicts the repo's own file 02. At most the explicit rate c min(1, 1/beta) and free-boundary specifics of Target T could be new. This must be checked before anyone works on Target T.

Required CLAIMS.md wording for C6:
- Setting: "d-dim lattice, finite volume (torus or box), bound uniform in volume; infinite-volume statement not written out".
- Claim: add "(the bound is exactly the Shen-Zhu-Zhu Bakry-Emery condition beta_SZZ < 1/(16(d-1)) in Hamiltonian form; the Hessian constant is sharp in d = 2); the identity is exact with the repo H_KS, the correction is >= 0, bounded per link and extensive, so it transfers no gap".
- Status: "proved (skeptic 010 passed with rewording)".
- STATE line 1 (Target T): add "check Piazza-Shen 2610.35286 first: LSI for 2D SU(N) at all beta, uniform in volume, likely already implies gap(H_FF) > 0 for all beta in d = 2".

## Open issues

- Rerun or drop the torus j_max = 3/2, g^2 = 3 interval (needs about 7 min and 4 GB, so it is a job for the overseer).
- Interval-arithmetic version for the ladder (small dims; mpmath/flint feasible) would remove the heuristic margin there.
- Read Piazza-Shen 2610.35286 itself: which Dirichlet form, which boundary condition, and how the constant depends on beta.
- Whether lambda_min(Hess S) >= -4d (pi-flux value) holds in d = 3. It would improve the d = 3 constant, but SZZ-equivalent bounds are not new anyway.
