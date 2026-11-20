# 20261120-026: easy half of min-max and the k-eigenvalue Schur count (C5) in Lean

## Question
Can the C5 eigenvalue-count statement ("Schur condition on a codimension-k subspace of the P-space implies at most k eigenvalues of H below E") be proved in Lean 4 + mathlib v4.36.0, which has no Courant-Fischer theorem?

## Method
Lean v4.36.0, mathlib v4.36.0 (pins unchanged, no `lake update`). New file `lean/MassGap/MinMax.lean` (namespace `MassGap.MinMax`), imported from `lean/MassGap.lean`. It uses mathlib's `LinearMap.IsSymmetric.eigenvectorBasis` / `eigenvalues` (eigenvalues listed with multiplicity, sorted decreasing) and the existing `MassGap.Schur.feshbach_form_bound` (work 021).

Proof of the min-max half: the span W of eigenvectors with eigenvalue < E has dimension #{i | λ_i < E} (orthonormal, so linearly independent). On W, `re⟪x,Tx⟫ = Σ λ_i |c_i|² < E‖x‖²` unless x = 0, so W ∩ V = 0 and dim W + dim V ≤ n.

Proof of the count: `feshbach_form_bound` is applied with P-side space V and Q-side space K (as subtypes), blocks of H built with mathlib's orthogonal projections. This gives the form of H ≥ E on V ⊔ K. Since V ⊥ K, dim(V ⊔ K) = dim V + dim K ≥ n - k, and the min-max half gives the count.

Build (repo root):

    python tools/pin.py bash -c 'export PATH="$HOME/.elan/bin:$PATH"; export LEAN_NUM_THREADS=1; cd lean && lake build'
    python tools/pin.py bash -c 'export PATH="$HOME/.elan/bin:$PATH"; cd lean && lake env lean ../work/20261120-026-lean-minmax/axioms_check.lean'

## Results
Setting for all statements: finite-dimensional inner product space M over 𝕜 = ℝ or ℂ (any `RCLike`), n = dim M, T or H symmetric (self-adjoint). This is the abstract finite-dimensional statement used for the lattice Hamiltonian at fixed lattice and truncation; no lattice, volume or gauge group enters.

- `card_lt_add_finrank_le` [proved] (Courant-Fischer, easy half): if `E‖x‖² ≤ re⟪x, T x⟫` for all x in a subspace V, then #{i | λ_i < E} + dim V ≤ n (eigenvalues counted with multiplicity).
- `card_eigenvalues_lt_le` [proved]: if moreover n ≤ dim V + k, then T has at most k eigenvalues < E.
- `le_eigenvalue_of_card_lt_le` [proved]: if at most k < n eigenvalues are < E, the (k+1)-th lowest eigenvalue (`eigenvalues ⟨n-1-k⟩` in mathlib's decreasing order) is ≥ E.
- `form_ge_on_sup` [proved]: Q-space K, P-space Kᗮ, V ≤ Kᗮ. Hypotheses (the 006 ones, stated with the form of H): floor `d‖v‖² ≤ re⟪v,Hv⟫` on K; coupling bound `‖⟪u,Hv⟫‖² ≤ re⟪u,Xu⟫ ‖v‖²` and `re⟪u,Xu⟫ ≥ 0` for u ∈ V, v ∈ K (X ≥ B Bᴴ on V); E < d; Schur condition `E‖u‖² ≤ re⟪u,Hu⟫ - re⟪u,Xu⟫/(d-E)` on V. Conclusion: `E‖w‖² ≤ re⟪w,Hw⟫` on V ⊔ K. Uses `feshbach_form_bound`.
- `schur_count` [proved] (C5, k-eigenvalue version): same hypotheses plus dim Kᗮ ≤ dim V + k (V has codimension ≤ k in the P-space). Conclusion: H has at most k eigenvalues < E, with multiplicity.
- `schur_count_sorted` [proved]: for k < n, E ≤ (k+1)-th lowest eigenvalue of H. With k = 1 this is the E1 lower bound used in C5.

`X` is any linear map on M; only its form on V is used. `H` need not preserve K or Kᗮ.

## Axioms check
`axioms_check.lean` prints the axioms of every theorem above and of `feshbach_form_bound`. All depend only on `propext`, `Classical.choice`, `Quot.sound` (no `sorryAx`). See `run.log`.

## Verdict
The C5 eigenvalue-count lemma (finite dimension, any k) is formally proved with no `sorry`, closing the gap noted in work 021. The min-max half needed is proved from `eigenvectorBasis`.

## Open issues
- Only the "at most k below E" half of min-max is proved; the full Courant-Fischer equality (λ_k = max over codim subspaces of min of Rayleigh quotient) is not.
- Floating-point aspects of C5 (rounding margin, Lanczos completeness) are outside this formalization: the lemma takes the form inequalities as hypotheses.
- `MinMax.lean` alone typechecks in about 30 s on one thread.

## Build log excerpt (`run.log`)
    error: unknown short option '-j'          (lake has no -j; rerun without it, LEAN_NUM_THREADS=1)
    ⚠ [5496/5500] Built MassGap.MinMax (27s)
    ✔ [5498/5500] Built MassGap (23s)
    Build completed successfully (5500 jobs).
    'MassGap.MinMax.card_lt_add_finrank_le' depends on axioms: [propext, Classical.choice, Quot.sound]
    'MassGap.MinMax.schur_count' depends on axioms: [propext, Classical.choice, Quot.sound]
    'MassGap.MinMax.schur_count_sorted' depends on axioms: [propext, Classical.choice, Quot.sound]

MinMax.lean emits only the mathlib "consider adding `module`" import warnings (harmless).
