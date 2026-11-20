import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Analysis.InnerProductSpace.Basic

/-!
# C5: Schur-complement (Feshbach) lower bound, finite dimension

Setting: `𝕜 = ℝ` or `ℂ` (any `RCLike`).

* `schur_lower_bound_iff` (matrix form, exact Schur complement, finite dimension): for
  `H = [[P, B], [Bᴴ, Q]]` with `Q - E > 0`,
  `H - E ≥ 0  ↔  P - E - B (Q - E)⁻¹ Bᴴ ≥ 0`. In particular (`schur_lower_bound`) the Schur
  complement condition implies `H` has no eigenvalue below `E` (the E0 statement).
* `feshbach_form_bound` (form version used in work 006, with floor `d` and coupling bound `X`):
  if `⟪v, D v⟫ ≥ d ‖v‖²`, `‖⟪u, B v⟫‖² ≤ ⟪u, X u⟫ ‖v‖²` (i.e. `X ≥ B Bᴴ`), `E < d`, and
  `A - X/(d - E) ≥ E` on a subspace `V` of the P-space, then the quadratic form of
  `H = [[A, B], [Bᴴ, D]]` is `≥ E ‖·‖²` on `V × Q`. If `V` has codimension `k`, so does `V × Q`;
  turning this into "at most `k` eigenvalues below `E`" needs the min-max principle,
  which mathlib v4.36.0 does not have (see work 021 README).
-/

open Matrix
open scoped InnerProductSpace ComplexOrder

namespace MassGap.Schur

section MatrixForm

variable {𝕜 : Type*} [RCLike 𝕜] {m n : Type*} [Fintype m] [Fintype n]
  [DecidableEq m] [DecidableEq n]

theorem fromBlocks_sub_smul_one (P : Matrix m m 𝕜) (B : Matrix m n 𝕜) (Q : Matrix n n 𝕜)
    (c : 𝕜) :
    fromBlocks P B Bᴴ Q - c • (1 : Matrix (m ⊕ n) (m ⊕ n) 𝕜) =
      fromBlocks (P - c • 1) B Bᴴ (Q - c • 1) := by
  rw [← fromBlocks_one, fromBlocks_smul]
  ext i j
  rcases i with i | i <;> rcases j with j | j <;> simp

/-- Exact Schur complement criterion at energy `E`. -/
theorem schur_lower_bound_iff (P : Matrix m m 𝕜) (B : Matrix m n 𝕜) (Q : Matrix n n 𝕜) (E : ℝ)
    (hQ : (Q - (E : 𝕜) • 1).PosDef) :
    (fromBlocks P B Bᴴ Q - (E : 𝕜) • 1).PosSemidef ↔
      (P - (E : 𝕜) • 1 - B * (Q - (E : 𝕜) • 1)⁻¹ * Bᴴ).PosSemidef := by
  have : Invertible (Q - (E : 𝕜) • 1) := hQ.isUnit.invertible
  rw [fromBlocks_sub_smul_one]
  exact PosDef.fromBlocks₂₂ _ _ hQ

/-- E0 version: `Q > E` and Schur complement `≥ E` imply `H ≥ E` (no eigenvalue below `E`). -/
theorem schur_lower_bound (P : Matrix m m 𝕜) (B : Matrix m n 𝕜) (Q : Matrix n n 𝕜) (E : ℝ)
    (hQ : (Q - (E : 𝕜) • 1).PosDef)
    (hS : (P - (E : 𝕜) • 1 - B * (Q - (E : 𝕜) • 1)⁻¹ * Bᴴ).PosSemidef) :
    (fromBlocks P B Bᴴ Q - (E : 𝕜) • 1).PosSemidef :=
  (schur_lower_bound_iff P B Q E hQ).2 hS

end MatrixForm

section FormVersion

variable {𝕜 : Type*} [RCLike 𝕜] {EP EQ : Type*}
  [NormedAddCommGroup EP] [InnerProductSpace 𝕜 EP]
  [NormedAddCommGroup EQ] [InnerProductSpace 𝕜 EQ]

/-- Scalar core of the bound: `|2 c t| ≤ x + t² a²` when `c² ≤ x a²`, `x ≥ 0`. -/
theorem two_abs_mul_le {c x a t : ℝ} (hx : 0 ≤ x) (hc : c ^ 2 ≤ x * a ^ 2) :
    |2 * c * t| ≤ x + t ^ 2 * a ^ 2 := by
  apply abs_le_of_sq_le_sq _ (by positivity)
  have h1 : 0 ≤ (x - t ^ 2 * a ^ 2) ^ 2 := sq_nonneg _
  have h2 : c ^ 2 * t ^ 2 ≤ x * a ^ 2 * t ^ 2 :=
    mul_le_mul_of_nonneg_right hc (sq_nonneg t)
  nlinarith [h1, h2]

/-- Feshbach-Schur form bound (Lemma of work 006; no finiteness used, linear maps need not be
bounded). The quadratic form of `H = [[A, B], [Bᴴ, D]]` at `(u, v)` is
`re⟪u, A u⟫ + 2 re⟪u, B v⟫ + re⟪v, D v⟫`. -/
theorem feshbach_form_bound (A X : EP →ₗ[𝕜] EP) (B : EQ →ₗ[𝕜] EP) (D : EQ →ₗ[𝕜] EQ)
    (V : Submodule 𝕜 EP) {d E : ℝ}
    (hD : ∀ v, d * ‖v‖ ^ 2 ≤ RCLike.re ⟪v, D v⟫_𝕜)
    (hX : ∀ u v, ‖⟪u, B v⟫_𝕜‖ ^ 2 ≤ RCLike.re ⟪u, X u⟫_𝕜 * ‖v‖ ^ 2)
    (hX0 : ∀ u, 0 ≤ RCLike.re ⟪u, X u⟫_𝕜)
    (hE : E < d)
    (hS : ∀ u ∈ V, E * ‖u‖ ^ 2 ≤ RCLike.re ⟪u, A u⟫_𝕜 - RCLike.re ⟪u, X u⟫_𝕜 / (d - E)) :
    ∀ u ∈ V, ∀ v, E * (‖u‖ ^ 2 + ‖v‖ ^ 2) ≤
      RCLike.re ⟪u, A u⟫_𝕜 + 2 * RCLike.re ⟪u, B v⟫_𝕜 + RCLike.re ⟪v, D v⟫_𝕜 := by
  intro u hu v
  have htpos : 0 < d - E := by linarith
  have hc : (RCLike.re ⟪u, B v⟫_𝕜) ^ 2 ≤ RCLike.re ⟪u, X u⟫_𝕜 * ‖v‖ ^ 2 := by
    have h1 : |RCLike.re ⟪u, B v⟫_𝕜| ≤ ‖⟪u, B v⟫_𝕜‖ := RCLike.abs_re_le_norm _
    have h2 : (RCLike.re ⟪u, B v⟫_𝕜) ^ 2 ≤ ‖⟪u, B v⟫_𝕜‖ ^ 2 := by
      rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) h1 2
    exact h2.trans (hX u v)
  set t := d - E with ht
  set x := RCLike.re ⟪u, X u⟫_𝕜 with hxdef
  set c := RCLike.re ⟪u, B v⟫_𝕜 with hcdef
  set a := ‖v‖ with hadef
  have hkey : 0 ≤ x / t + t * a ^ 2 + 2 * c := by
    have h := two_abs_mul_le (t := t) (hX0 u) hc
    have h' : -(x + t ^ 2 * a ^ 2) ≤ 2 * c * t := by
      have := neg_abs_le (2 * c * t); linarith
    have hmul : 0 ≤ (x / t + t * a ^ 2 + 2 * c) * t := by
      have e1 : (x / t + t * a ^ 2 + 2 * c) * t = x + t ^ 2 * a ^ 2 + 2 * c * t := by
        field_simp
      rw [e1]; linarith
    exact nonneg_of_mul_nonneg_left hmul htpos
  have hSu := hS u hu
  have hDv := hD v
  have hd : d = t + E := by rw [ht]; ring
  rw [hd] at hDv
  nlinarith [hSu, hDv, hkey]

end FormVersion

end MassGap.Schur
