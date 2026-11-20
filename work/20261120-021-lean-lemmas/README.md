# 20261120-021: first Lean lemmas (C4 abstract step, C5 Schur bound, C2 coefficient)

## Question
Can the abstract finite-dimensional steps behind C4 (Bloch decomposition gives no in-sector gap bound), C5 (Schur-complement / Feshbach lower bound) and the C2 second-order coefficient be checked in Lean 4 + mathlib, with no `sorry`?

## Method
Lean v4.36.0, mathlib v4.36.0 (pins unchanged, no `lake update`). New files in `lean/MassGap/`, imported from `lean/MassGap.lean`:
- `BlochObstruction.lean` (namespace `MassGap.Bloch`)
- `SchurBound.lean` (namespace `MassGap.Schur`)
- `StrongCoupling.lean` (namespace `MassGap.StrongCoupling`)

`axioms_check.lean` (this folder) prints the axioms of every main theorem. All depend only on `propext`, `Classical.choice`, `Quot.sound` (no `sorryAx`).

Build (repo root):

    python tools/pin.py bash -c 'export PATH="$HOME/.elan/bin:$PATH"; cd lean && lake build'
    python tools/pin.py bash -c 'export PATH="$HOME/.elan/bin:$PATH"; cd lean && lake env lean ../work/20261120-021-lean-lemmas/axioms_check.lean'

Cold build of the three new files: about 5 min on 3 cores (BlochObstruction 158 s, SchurBound 78 s, StrongCoupling 47 s).

## Results

### 1. C4 abstract step (`BlochObstruction.lean`): proved
Setting: finite dimension; the spectrum lemma over any field, the gap lemmas for real symmetric matrices. A two-block matrix `fromBlocks A 0 0 D` stands in for `H = ⊕_θ H_θ` with two fibers. "Gap" means distance from the lowest eigenvalue to the next distinct eigenvalue (`IsGroundGap S e Δ`: `Δ > 0`, `e ∈ S`, `e + Δ ∈ S`, every `x ∈ S` is `e` or `≥ e + Δ`; multiplicities ignored).
- `spectrum_fromBlocks_diag` [proved]: `spectrum K (fromBlocks A 0 0 D) = spectrum K A ∪ spectrum K D` (A1, finite-dimensional).
- `band_data_do_not_bound_in_sector_gap` [proved]: for every `f : ℝ → ℝ → ℝ` with `f > 0` and every pair of ground energies `(e₁, e₂)`, there are real symmetric 2x2 blocks with exactly these ground energies and block-1 gap `< f e₁ e₂`. So band data (the blocks' ground energies) do not bound the in-sector gap from below.
- `twoBlock_counterexample` [proved]: `H(δ) = diag(e, e+δ) ⊕ diag(e, e+1)` is symmetric, both blocks have ground energy `e`, block gaps `δ` and 1, any `δ > 0`; `spectrum_twoBlock` gives spec H = {e, e+δ} ∪ {e, e+1}.
- `in_block_gap_of_band_gap` [proved] (A2 positive direction): if spec(A ⊕ D) ⊆ (-∞, b₀] ∪ [b₁, ∞) and A has exactly one level `e ≤ b₀` in the lowest band, every other level `x` of A has `x - e ≥ b₁ - b₀`.
- `overlapping_bands_example` [proved]: the converse fails; blocks diag(0,1), diag(2,3) both have gap 1, bands overlap.

### 2. C5 Schur-complement bound (`SchurBound.lean`): E0 version proved, eigenvalue-count version partial
Setting: finite dimension for the matrix form, `𝕜 = ℝ` or `ℂ` (any `RCLike`).
- `schur_lower_bound_iff` [proved]: for `H = [[P, B], [Bᴴ, Q]]` and real `E` with `Q - E·1` positive definite, `H - E·1 ≥ 0 ⟺ P - E·1 - B (Q - E·1)⁻¹ Bᴴ ≥ 0` (PosSemidef). Uses mathlib `Matrix.PosDef.fromBlocks₂₂`.
- `schur_lower_bound` [proved]: the `⇐` direction, the E0 statement (no eigenvalue of H below E).
- `feshbach_form_bound` [proved]: the Lemma of work 006 in form language, for arbitrary inner-product spaces (no finiteness used). Hypotheses: `re⟪v,Dv⟫ ≥ d‖v‖²`; `‖⟪u,Bv⟫‖² ≤ re⟪u,Xu⟫ ‖v‖²` and `re⟪u,Xu⟫ ≥ 0` (i.e. `X ≥ B Bᴴ`); `E < d`; `re⟪u,Au⟫ - re⟪u,Xu⟫/(d-E) ≥ E‖u‖²` for all `u` in a submodule `V`. Conclusion: `re⟪u,Au⟫ + 2 re⟪u,Bv⟫ + re⟪v,Dv⟫ ≥ E(‖u‖² + ‖v‖²)` for all `u ∈ V` and all `v`, i.e. the form of H is `≥ E` on `V × Q`.
- Not formalized: "Schur complement has at most k eigenvalues below E ⇒ H has at most k". The remaining step is: form `≥ E` on a codimension-k subspace ⇒ at most k eigenvalues below E (the easy half of min-max). Obstacle: mathlib v4.36.0 has no Courant-Fischer / min-max characterization of `Matrix.IsHermitian.eigenvalues` or `LinearMap.IsSymmetric.eigenvalues` (searched: no "Courant", "min-max"); only the extreme eigenvalues via Rayleigh quotients (`LinearMap.IsSymmetric.hasEigenvalue_iSup_of_finiteDimensional`, `hasEigenvalue_iInf_of_finiteDimensional`), plus `eigenvectorBasis` and `eigenvalues_antitone`. The step can be proved by hand from `eigenvectorBasis` (span of eigenvectors with eigenvalue < E meets the subspace only in 0, then a dimension count) but was not done within this budget. The product space `V × Q` would also need `WithLp 2 (EP × EQ)` to be an inner product space; the form lemma avoids this by stating the form explicitly.

### 3. C2 coefficient (`StrongCoupling.lean`): arithmetic proved, inputs cited
Setting: SU(2) KS single plaquette, lattice, finite volume, `x = g²`. Inputs taken from 003/005, not formalized: unperturbed energies `2x j(j+1)` and matrix elements `⟨j±1/2|V|j⟩ = -2/x`.
- `firstOrder_gap` [proved]: `E_{1/2} - E_0 = 3x/2`.
- `secondOrder_ground` [proved]: `δE_0 = -8/(3x³)`; `secondOrder_first` [proved]: `δE_{1/2} = 16/(15x³)`.
- `secondOrder_gap` [proved]: `Δ⁽²⁾ = 56/(15 x³) = 56/(15 g⁶)`, for `x ≠ 0`.
The fourth-order term (−20144/3375) and the Rayleigh-Schrodinger formula itself are not formalized.

## Verdict
C4 abstract steps A1/A2 (finite-dimensional, two fibers) and the C5 E0 lower-bound lemma (matrix and form versions) are formally proved. The C5 k-th eigenvalue version is partially proved (form bound on a codimension-k subspace; the min-max step is missing). C2's 56/15 is checked given the cited inputs.

## Open issues
- Prove the easy half of min-max from `Matrix.IsHermitian.eigenvectorBasis` to close the C5 eigenvalue-count version.
- C4 statements cover finitely many fibers in finite dimension; the direct-integral version (Reed-Simon XIII.16) is not formalized.
- Mathlib emits "consider adding `module`" warnings for every import; harmless.

## Build log excerpt (`run.log`)

    error: MassGap/BlochObstruction.lean:65:70: unsolved goals      (build 2, fixed)
    ⚠ [4698/4702] Built MassGap.SchurBound (78s)
    ⚠ [4698/4702] Built MassGap.BlochObstruction (158s)
    ✔ [4700/4702] Built MassGap (22s)
    Build completed successfully (4702 jobs).
    'MassGap.Schur.feshbach_form_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
    'MassGap.StrongCoupling.secondOrder_gap' depends on axioms: [propext, Classical.choice, Quot.sound]
