# 20261120-006: Certified two-sided bounds on the SU(2) KS gap (untruncated, finite lattice)

## Question

The iteration-1 ED (`work/20261120-003-su2-ks-ed/`) gives E0, E1 of the Kogut-Susskind Hamiltonian truncated to link spins j <= j_max. Can these be turned into proved intervals [lo, hi] for E0, E1 and Delta = E1 - E0 of the **untruncated** Hamiltonian on the same finite lattice at fixed coupling? Such intervals are the input that Knabe / Lemm-Mozgunov finite-size gap criteria need.

Setting of every number below: 2+1D Hamiltonian lattice, gauge group SU(2), lattice spacing a = 1, finite lattice (single plaquette; open 2-plaquette ladder; periodic 2x2), fixed g^2, H as in `src/ks_su2.py`:

    H = (g^2/2) sum_l j_l(j_l+1) - (2/g^2) W,   W = sum_p Tr U_p   (fundamental trace)

A certified interval here is a statement about this finite lattice at this g^2. It says nothing about the continuum or infinite volume.

## Sectors (what "the gap" means)

H acts on the gauge-invariant space G (Gauss's law at every vertex). G splits into Z2 electric-flux sectors: the parity of sum 2j_l over the links crossing a noncontractible dual cycle is conserved by every Tr U_p, and Gauss's law makes it independent of which homologous cut is used.

- plaq1, ladder2 (open, planar): every Z2 cycle is a sum of plaquette boundaries, so there is one sector, S = G. The ladder gap covers both reflection parities (the ladder E1 is P-even at all couplings sampled in iteration 1).
- periodic2x2: S = zero-flux sector (both cut parities even), all momenta, no point-group resolution. Cuts: x-links at x = 0 and y-links at y = 0.

Delta = lambda_1(H|S) - lambda_0(H|S), eigenvalues counted with multiplicity. lambda_0 is simple (positivity improvement, iteration-1 README), so Delta > 0 iff lambda_1 > lambda_0. On the torus lambda_1 may belong to a degenerate multiplet at nonzero momentum.

**Completeness of the truncated basis (checked).** The truncated space P must be all of S with spins <= j_max, otherwise the electric floor on the complement fails. `src/certify.py` enumerates every link-spin configuration admitting an invariant at each vertex (SU(2) invariant exists iff sum 2j is even and 2 max j <= sum j, proved by the iterated triangle rule), filters by flux, and takes all intertwiner labels. `build_sector` asserts Tr U_p never leaves the enumerated set. In every case run, the enumerated sector has the same dimension as the iteration-1 BFS basis from the strong-coupling vacuum (`checks.txt`, "vacuum-connected"). So the iteration-1 basis was already the full zero-flux sector, on the torus too (23 / 622 / 9418 at j_max = 1/2, 1, 3/2).

## Method

Write S = P (+) Q, with P = span of spin networks in S with all spins <= j_max (finite-dimensional) and Q = its orthogonal complement. In block form H = [[A, B], [B*, D]] with A = PHP, B = PHQ, D = QHQ.

### Upper bounds (proved, cited)

lambda_k(H|S) <= lambda_k(A). Min-max principle with the trial subspace P (Reed-Simon IV, Thm XIII.1; Rayleigh-Ritz). For P too large for dense LAPACK, Lanczos vectors are orthonormalized and the Ritz values on their span are used; those are again min-max upper bounds for whatever subspace was found.

### Lower bounds: Schur-complement (Feshbach) bound (proved here)

**Lemma.** Suppose D >= d (as forms on Q) and X >= B B* (as forms on P), X finite. Let E < d and suppose lambda_k(A - X/(d-E)) >= E. Then lambda_k(H|S) >= E.

*Proof (proved).* For phi = u + v, u in P, v in Q in the form domain:
<phi,(H-E)phi> = <u,(A-E)u> + 2 Re<u,Bv> + <v,(D-E)v> >= <u,(A-E)u> - 2||B*u|| ||v|| + (d-E)||v||^2 >= <u,(A-E)u> - ||B*u||^2/(d-E) >= <u,(A - E - X/(d-E))u>, minimizing over ||v|| in the second step. By hypothesis A - X/(d-E) - E >= 0 on a subspace V of P of codimension <= k (span of its eigenvectors k, k+1, ...). Then H - E >= 0 on V (+) Q, of codimension <= k in S, and min-max gives lambda_k(H|S) >= E. QED.

**Choice of E (proved).** Let U >= lambda_k(A) be the upper bound above and f(E) = lambda_k(A - X/(d-E)). f is nonincreasing on E < d (X >= 0) and f(E) <= lambda_k(A). If U < d, set L = f(U) (or any computed lower bound of it). Then L <= f(U) <= lambda_k(A) <= U, so f(L) >= f(U) >= L, and the Lemma gives **lambda_k(H|S) >= L**. Output: E_k in [L_k, U_k], Delta in [L_1 - U_0, U_1 - L_0]. If U >= d, no lower bound (reported as -inf).

This is the standard Feshbach-Schur / Schur-complement estimate (cited for the mechanism: Bach, Frohlich, Sigal, Adv. Math. 137 (1998) 299 (isospectral Feshbach map); Griesemer-Hasler, J. Funct. Anal. 254 (2008)); the version above is self-contained and proved directly. It was chosen over Temple/Lehmann because those need an a priori lower bound on lambda_{k+1}, which is what we are trying to prove; here the only a priori inputs are d and X, both explicit.

### The inputs d and X

**Coupling X.** B = -(2/g^2) PWQ. Tr U_p multiplies four loop links by D^{1/2}, changing each spin by +-1/2 (Clebsch-Gordan series), so QWP maps P only into the shell P_+ = {states with max spin = j_max + 1/2}. Hence B B* = (4/g^4) W_{P,shell} W_{shell,P} **exactly**, computable once the sector basis is built up to j_max + 1/2. This is used for every j_max below the largest one built ("schur-exact-shell").
At the largest j_max built, a crude bound is used ("schur-crude-shell", proved): Q Tr U_p psi = Q Tr U_p Pi_p psi, with Pi_p projecting on states having a spin-j_max link in p, and ||Tr U_p|| <= 2 (multiplication by |Tr U| <= 2 on L^2). Cauchy-Schwarz over plaquettes: ||QWPpsi||^2 <= 4 N_p <psi, n psi>, n = diag(# plaquettes with a spin-j_max link). So X = (4/g^4) 4 N_p n.

**Floor d.** Two proved bounds, take the larger:
1. Crude: D >= (g^2/2) c_Q - (2/g^2) ||W|| >= (g^2/2) c_Q - 4 N_p/g^2. c_Q is a lower bound on sum_l j_l(j_l+1) over gauge-invariant configurations (any flux sector) with some spin > j_max. `casimir_floor` computes the exact minimum by enumeration (budget 2e6 configurations) and stops when the analytic bound below for larger max spin exceeds it; otherwise it uses the analytic bound. Analytic bound (proved): a link of spin J forces sum of the other spins at each endpoint >= J (Gauss: 2 max <= sum), and with k other legs Cauchy-Schwarz gives sum j(j+1) >= J^2/k + J. Both endpoints count if their other-leg sets are disjoint (not true on the 2x2 torus, which has parallel links), else the better one.
2. Nested ("nested"): Q = Q1 (+) Q2, Q1 = basis states with j_max < max spin <= j_top (the sector basis is complete up to j_top), Q2 the rest. Apply the Lemma with k = 0 to QHQ on Q1 (+) Q2, with the crude X and crude floor of Q2. This gives inf spec(QHQ) >= d_nested. Used only when dim Q1 <= 2500 (dense).

### Numerical honesty

- Dense eigenvalues (LAPACK evr, n <= 2500) are backward stable. 10 n eps ||A||_inf is subtracted from every lower and added to every upper bound. This is a floating-point margin, **not interval arithmetic**: the bounds are proved for the computed matrices up to rounding. The matrix elements themselves were validated pointwise to 1e-14 in iteration 1.
- For n > 2500 `src/certify.py` also has an inertia route (Sylvester's law with an LDL^T factorization and an explicitly computed residual, Weyl correction), sparse (SuperLU, diagonal pivoting) and dense (Bunch-Kaufman). It was tested on the torus at n = 9418 (sparse: inertia correct, residual 1.6e-12, but 273 s per factorization; dense: > 3 min) and is not used in the final run.

## Validation

Single plaquette: the iteration-1 j_max = 40 ED values (agreeing with the independent group-manifold solver to 1.2e-8) lie inside every certified interval at every g^2 and j_max (`checks.txt`: "inside every certified interval: yes"). Intervals shrink with j_max and collapse to ~1e-9 width once j_max is converged.

## Result

Tightest certified interval for Delta = lambda_1 - lambda_0 per lattice and coupling (j_max at which it was obtained in brackets). Each is a proved two-sided bound for the untruncated H on the stated sector, up to floating-point rounding (see above). "none" = no lower bound (floor d not above lambda_1(A)).

| g^2 | plaq1 (S = G) | ladder2 (S = G) | periodic2x2 (zero flux, all momenta) |
|---|---|---|---|
| 0.2 | [3.8715812, 3.8715812] (16) | none at j_max <= 6; upper 43.0 is meaningless | none |
| 0.5 | [3.6635501, 3.6635501] (8) | [3.2936633, 3.2936633] (6), width 1e-8 | none |
| 1.0 | [3.2690992, 3.2690992] (4) | [3.0620125, 3.0620125] (6), width 8e-10 | none; upper bounds only (Delta not bounded) |
| 1.5 | [3.1172414, 3.1172414] (4) | [2.9980864, 2.9980864] (4) | none |
| 2.0 | [3.4262886, 3.4262886] (4) | [3.3545853, 3.3545853] (4) | none (j_max = 3/2 run not completed) |
| 3.0 | [4.6356258, 4.6356258] (4) | [4.6092285, 4.6092285] (4) | [4.08345, 4.13804] (1); [4.1257, 4.1281] (3/2, dev run, see below) |
| 5.0 | [7.5297906, 7.5297906] (2) | [7.5237504, 7.5237504] (2) | [7.234803, 7.235425] (1) |
| 10.0 | [15.003733, 15.003733] (2) | [15.002971, 15.002971] (2) | [14.860549, 14.860553] (1) |

Plaquette and ladder interval widths are 1e-12 to 1e-10 unless stated; digits beyond those shown are in `results.csv`. All plaquette and ladder values agree with the iteration-1 converged ED to the printed digits, and the torus g^2 >= 5 intervals contain the iteration-1 j_max = 3/2 values (7.2351, 14.860551).

Why the torus fails at g^2 <= 2 at j_max <= 1: the crude floor d = (g^2/2) c_Q - 16/g^2 is below lambda_1(A) (c_Q = 7.5 at j_max = 1; d = -5 at g^2 = 2). At j_max = 3/2 (c_Q = 11) the floor clears lambda_1 for g^2 >= 2, but the dense route is too slow at dim 9418 and the sparse LDL route takes ~2-4 min per factorization; it was run once, at g^2 = 3, in a development run of the same code path (crude shell, crude floor d = 11.17, c_Q = 11, LDL residual 1.6e-12, 413 s): E0 in [-0.3931522, -0.3931387], E1 in [3.7325523, 3.7349347], so Delta in [4.1257, 4.1281]. That row is **not** in `results.csv`; `torus_j32_lower.py` reproduces it (`python ... torus_j32_lower.py 3.0 2.0`, ~7 min per coupling) and its CSV is merged by `certify_gaps.py`. The background rerun for g^2 = 2, 3 died with its parent shell and was not repeated within the time budget.


Full table (all j_max, E0/E1 intervals, d, c_Q): `results.csv`. Log: `run.log`, `checks.txt`.

## Verdict

- Proved (modulo floating-point rounding, not interval arithmetic), finite lattice, fixed g^2, 2+1D SU(2) KS, a = 1: the untruncated Hamiltonian has a gap in the stated sector with two-sided bounds as tabulated. Single plaquette g^2 in [0.2, 10] and open ladder g^2 in [0.5, 10]: Delta known to ~1e-9 or better, so the iteration-1 truncated numbers are now certified values of the full-H gap. Periodic 2x2 zero-flux sector: certified only for g^2 >= 3 (width 5e-2 at g^2 = 3 from j_max = 1, 2.4e-3 from the j_max = 3/2 dev run; 6e-4 at g^2 = 5; 4e-6 at g^2 = 10).
- The method is the Schur-complement (Feshbach) Lemma above with an exactly computed shell coupling and an explicit electric floor. Its only a priori input is the floor d on the complement, which is elementary (Gauss law + |Tr U| <= 2). Temple/Lehmann were not needed.
- Validation: on the plaquette, the independently known spectrum lies inside every interval at every j_max.
- What it does not do: nothing here bounds the gap uniformly in volume or as g -> 0. These certified finite-lattice gaps are the local-gap input a Knabe or Lemm-Mozgunov criterion would consume, but those criteria need a frustration-free (or controlled) Hamiltonian, which KS is not. The weak-coupling torus, where torons should drive Delta toward ~g^2, is still out of reach: the crude floor needs j_max >= 2 there, which needs momentum blocking.

## Open issues

1. Torus j_max = 3/2 lower bounds and j_max = 2: need momentum / point-group blocking (intertwiner-basis rotation at 4-valent vertices) to bring blocks to <= ~2500 so dense inertia counts are cheap. Not done here (task 3 not reached). The inertia code path exists.
2. Floating-point: not interval arithmetic. For a computer-assisted proof, redo the final eigenvalue/inertia steps with an interval library (e.g. mpmath intervals or `python-flint` arb matrices) on the small matrices involved; all dims used are <= ~800 for certified rows.
3. The crude floor ||W|| <= 2 N_p and the crude shell X lose most at weak coupling. A sharper floor would come from a multi-shell (iterated) Schur bound; the nested floor is the first step of that.
4. Gap per symmetry sector (momentum, ladder reflection parity, torus J^PC) is not certified separately; the intervals are for the whole sector S as defined above.
5. The Lemma needs P to be the full truncated sector. That was verified by enumeration here, but any new lattice must repeat the completeness check (`build_sector` asserts closure; compare with `vacuum_component_size`).
