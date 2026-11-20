# 20261120-024: Sylvester inertia certification of the eigenvalue-counting step in C5, C9, C11

## Question

The Schur-complement (Feshbach) lower bounds of work 006, 012 and 015 need lambda_k of a computed matrix M = A - X/(d - E) to be at least E. When that eigenvalue came from an iterative solver (eigsh), a skipped eigenvalue would make the bound too high (skeptic 018). Can the step be replaced by a count that does not depend on the solver finding every low eigenvalue?

Setting of every number: as in the source work. SU(2) Kogut-Susskind, 2+1D Hamiltonian lattice, a = 1, finite lattice (open 2-plaquette ladder; periodic 2x2 torus), fixed g^2 (C5, C9). 1+1D O(3) rotor chain, finite chain, M = 0 sector, fixed g; spin-1 AKLT open chain, S_z = 0 sector (C11). Nothing here is a continuum or infinite-volume statement except the cited Gosset-Mozgunov step, which is unchanged.

## Method

**What is checked.** For each recorded lower bound E of level k (0-based) the Lemma hypothesis of work 006 is tested directly: E < d and

    N(E) := #{eigenvalues of A - X/(d - E) below E} <= k.

This is exactly what the Lemma needs, independent of how E was found (standard choice, bisection, Lanczos). For C5 the matrices are rebuilt with `src/certify.py` and the floor d is itself re-checked (the nested floor is a k = 0 instance of the same Lemma; if only a weakened value verifies, the weakened d is used for E0, E1). For C9 A, X come from the work-012 block pickles (extracted by `extract_torus.py`) and d is recomputed (split floor; matches the recorded `d_floor` to 1e-8). For C11 rotor rows the depth-1 sub-certificates that enter d (eigsh in 015) are checked first, then E0, E1. For AKLT the statement is #eig(H|S_z=0 < eps_claimed) <= 2 (the 2 kernel states in S_z = 0).

**Why the factorization gives the count** (`src/inertia.py`). Sylvester's law of inertia: a congruence K' = P L D L^T P^T with L unit lower triangular (nonsingular) and P a permutation has the same numbers of negative, zero and positive eigenvalues as D. D has 1x1 and 2x2 blocks; a 2x2 block has one negative eigenvalue if det < 0 and two of the diagonal's sign if det > 0, so #neg(D) is read off exactly.
- Dense (n <= 5000 for KS, n <= 1500 for rotor): LAPACK sytrf (Bunch-Kaufman symmetric pivoting) via `scipy.linalg.ldl`. Symmetric pivoting is a congruence, so it gives inertia.
- Sparse: SuperLU with `diag_pivot_thresh = 0`, `SymmetricMode`, MMD ordering on A^T + A. The code checks `perm_r == perm_c`, so the factorization is P^T K P = L U with the same permutation on both sides, and L diag(U) L^T is a congruence of the permuted matrix. Plain LU with partial (row-only) pivoting is not a congruence and the signs of diag(U) do not give the inertia; that path is never used. Without pivoting this is not backward stable for indefinite K, so stability is not assumed: the residual is computed.
- Error control (floating point, not interval arithmetic). R = P L D L^T P^T - (A - sI) is formed explicitly. ||R||_2 <= ||R_computed||_F + gamma_m ||B||_1 ||B||_inf, B = |L| |D|^(1/2) (diagonal majorant for 2x2 blocks), m the longest inner product; the second term is the standard bound (Higham, Ch. 3) on rounding in forming the product. Weyl: if r <= delta at shift s + delta, then N(s) <= #neg(D); at s - delta, N(s) >= #neg(D). `inertia_below(A, s)` returns this bracket [n_lo, n_hi]; `at_most(A, s, k)` returns "verified" if n_hi <= k, "failed" if n_lo > k, else "inconclusive". If inconclusive (an eigenvalue within delta of the recorded bound) it retries at s - 3 delta and reports "verified at claimed - w". The C5 and C9 runs used an earlier, looser rounding term (2 n eps || |L||D||L^T| ||_F); that only enlarges delta, so their verdicts stand, but several of their "weakened" verdicts would likely be plain "verified" with the final bound.
- Self-test (`selftest.py`, `selftest.log`): random symmetric matrices n = 50-1500 (dense and sparse paths) and a 2D Laplacian with known spectrum, all counts exact; a shift placed on an eigenvalue returns the bracket [10, 11], never a wrong count.

## Where Lanczos was actually used

| claim / work | lower-bound eigenvalue step | Lanczos? |
|---|---|---|
| C5 / 006 ladder (dim <= 616 for P, 752 for Q1), torus j_max <= 1 (622) | dense LAPACK `eigvalsh` | no |
| C5 / 006 torus j_max = 3/2 dev row | sparse LDL^T inertia already (`certify.ldl_count`) | no; row superseded by C9 |
| C9 / 012 all flux sectors, refine.csv (P blocks <= 2472) | dense LAPACK | no. Upper bounds from j_max = 2 blocks use eigsh + Rayleigh-Ritz, which are min-max upper bounds for any subspace (no completeness needed) |
| C11 / 015 rotor, dim > 1500 (4332 ... 184756) | eigsh, incl. the depth-1 floor sub-certificates | yes |
| C11 / 015, 018 AKLT eps_12 (dim 73789) | eigsh | yes; eps_10 (8953) dense, no |

So the skeptic's concern applies to C11 only. C5 and C9 were nevertheless re-checked by inertia: it replaces "LAPACK eigenvalue minus a heuristic 10 n eps ||A|| margin" with "exact pivot count of a factorization whose residual is measured".

## Results (`inertia.csv`; per-part files `inertia_c5.csv`, `inertia_c9.csv`, `inertia_c11.csv`, `inertia_c11_aklt.csv`)

Columns: claim, work, system, g (g^2 for KS, g for the rotor), j_or_l_max, level_k, quantity, claimed_lower, dim, lanczos_used, n_below_lo/n_below_hi (bracket of N(E); empty = side not needed), allowed (= k), delta, residual, method, verdict, note.

| claim | system (rows) | Lanczos in source | verified | verified with weakening w | failed | inconclusive | not checked / too large |
|---|---|---|---|---|---|---|---|
| C5 (006) | open ladder, j_max 1-6, g^2 0.2-10: floor d, E0, E1 (77 rows, dims 11-752) | no | 24 | 53 (w <= 1.4e-7, at g^2 = 10 floor; <= 1e-8 on the quoted g^2 = 1, j_max = 6 E0, E1) | 0 | 0 | 0 |
| C5 (006) | 2x2 torus zero flux, j_max 1/2, 1, g^2 2-10 (11 rows, dims 23-622) | no | 8 | 3 | 0 | 0 | 0 |
| C9 (012) | 2x2 torus, all flux sectors x 4 momenta, g^2 1-3, + refine.csv g^2 = 1.5 (43 rows, dims 98-2472) | no (upper bounds: Ritz, not affected) | 39 | 4 (g^2 = 3 zero-flux E0 lower bounds of the 4 momentum blocks, w <= 3.2e-8) | 0 | 0 | 0 |
| C11 (015) | rotor L = 3 open/ring, l_max 5, g = 2 and 0.35 (dim 4332) | yes | 7 | 1 (L = 3 ring g = 2 E0, w = 3.5e-9) | 0 | 0 | - |
| C11 (015) | rotor L = 4 (8092), L = 5 (8953), L = 8 ring (12870): key rows | yes | - | - | - | - | 18 not checked (time) |
| C11 (015) | rotor L = 6 (73789), L = 10 (184756): key rows | yes | - | - | - | - | 10 too large to factor (5 GB / time) |
| C11 (015/018) | AKLT open chain eps_12 (dim 73789, S_z = 0) | yes | - | - | - | - | too large (the n = 10, dim 8953 sparse factorization did not finish in ~9 min) |

Rotor "key rows" = per (L, bc, l_max) the smallest-g row with a positive gap lower bound (the 015 summary table) and the g = 2 row. A further 69 finite rotor lower bounds at dim > 1500 (other g) were not attempted. Rotor rows with dim <= 1500 used dense LAPACK in 015 (no Lanczos) and were not re-checked.

Each "weakened" row: the inertia count puts an eigenvalue of A - X/(d - E) within the residual bound (~1e-9) of the recorded E, so the recorded value cannot be confirmed at full precision; E - w is verified. For C5 these weakenings come from the earlier, looser rounding term (see Method). Example rows: ladder g^2 = 1, j_max = 6: E0 >= -3.348020165 - 1.0e-8, E1 >= -0.2860068833 - 9.6e-9 (the certified Delta = 3.0620125 interval widens from 8e-10 to ~2e-8 at worst); torus g^2 = 3 zero flux k = 0: E0 >= -0.393139756 - 3.2e-8, E1 >= 3.734794588 verified exactly.

## Verdict

- **No certified lower bound failed.** Every checked row has at most k eigenvalues of A - X/(d - E) below E, in some cases only after lowering E by w <= 1.4e-7.
- **Proved modulo floating-point rounding of an explicitly residual-checked LDL^T factorization (no solver-completeness assumption):**
  - C5 ladder, all g^2 and j_max in `work/20261120-006/results.csv` with a finite lower bound; quoted intervals hold at their printed precision (Delta = 3.0620125 at g^2 = 1 holds with width <= ~2e-8 instead of 8e-10 at the weakest reading).
  - C5 torus j_max <= 1 rows (g^2 = 3, 5, 10).
  - C9: all certified intervals (all-sector g^2 = 3 [2.02927, 2.03424]; zero flux g^2 = 3 [4.12793, 4.12801], g^2 = 2 [2.6055, 2.6258], g^2 = 1.5 [0.7694, 1.7016]) are reconfirmed at printed precision. These never depended on Lanczos for the lower side.
  - C11 rotor L = 3 open and ring, l_max = 5, at g = 2 and at the smallest certified g = 0.35 (open [0.13109, 0.13204], ring [0.12446, 0.12824]). Their dependence on Lanczos completeness is removed.
- **Still conditional on Lanczos completeness:**
  - C11 rotor rows with dim 8092-184756: L = 4 (g down to 0.5 open, 0.7 ring), L = 5 (g = 1), L = 6 (g = 1 ring, 3 open), L = 8 and L = 10 (g = 2). This includes the claim's "only g = 2 at L = 8-10".
  - The AKLT eps_12 = 0.364356 input, hence the GM >= 0.3026 and Knabe >= 0.3008 numbers. The Lanczos-free AKLT statement remains GM >= 0.30107 from dense n = 10 (work 018).
- Suggested CLAIMS wording change (overseer): C5 ladder/torus and C9 "numerical-certified (floating point, heuristic rounding margin)" can become "certified up to floating-point rounding, eigenvalue counts by LDL^T inertia with measured residual (work 024)". C11: "Lanczos completeness assumed" can be dropped for L <= 3 only.

## Open issues

1. Not interval arithmetic. The residual R and the rounding term are evaluated in floating point with the standard rounding model; an enclosure (e.g. python-flint arb or mpmath intervals on the LDL^T product) would make the C5/C9 rows fully computer-assisted. Dims <= 2472 make that feasible.
2. Rotor rows not checked: see table (time; the sparse path took 15-400 s per count at dim 4332 and the L = 4 depth-1 floor needs 3 extra counts per row). Large rows (dim 73789, 184756) need momentum and reflection blocking before any factorization fits in 5 GB and the time budget; the sparse no-pivot SuperLU fill was already 2.9e6 for the dim-4332 ring.
3. The sparse path is unpivoted. It worked here (residuals 1e-11 to 1e-9), but a symmetric indefinite sparse solver with Bunch-Kaufman pivoting (e.g. MUMPS/PARDISO LDL^T) would be the robust choice; none is in `requirements.txt`.
4. The weakened verdicts (eigenvalue within ~1e-9 of the recorded bound) change certified intervals only in the 8th digit or beyond; quoted intervals in CLAIMS.md are unaffected at their printed precision, but the full-precision `results.csv` values should be read as valid to the stated weakening.
5. Courant-Fischer/Sylvester in Lean are still missing, so none of this is formalized.

## Reproduce (repo root; 1 thread each)

    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/selftest.py
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c5_ladder_torus.py      # ~4 min
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/extract_torus.py        # needs 012 artifacts (012/build.py), peak < 3.5 GB (monitored)
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c9_torus.py             # ~15 min
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c11_rotor.py 1501 20000 --key
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c11_aklt.py 10 12
    python work/20261120-024-inertia-certification/merge.py

Logs: `run_c5.log`, `run_c9.log`, `run_c11_rotor.log`, `run_c11_aklt.log`, `selftest.log`, `extract.log`.
