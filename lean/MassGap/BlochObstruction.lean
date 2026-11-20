import Mathlib.LinearAlgebra.Eigenspace.Matrix
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Analysis.RCLike.Basic

/-!
# C4, abstract step: band data do not bound the in-sector gap

Finite-dimensional version of steps A1/A2 of work 008 (claim C4).
Setting: matrices over a field (spectrum lemma) or real symmetric matrices (gap lemmas),
finite dimension. A block-diagonal matrix `fromBlocks A 0 0 D` is the finite-dimensional
stand-in for a Bloch decomposition `H = ⊕_θ H_θ` with two fibers.

* `spectrum_fromBlocks_diag`: spec (A ⊕ D) = spec A ∪ spec D (A1).
* `band_data_do_not_bound_in_sector_gap`: no positive function of the blocks' ground
  energies bounds the in-block gap from below (A1, "Bloch gives no gap bound").
* `twoBlock_counterexample`: two 2x2 blocks with identical ground energy `e`,
  block 2 has gap 1, block 1 has gap `δ`, any `δ > 0`.
* `in_block_gap_of_band_gap`: a uniform gap above the lowest band implies an in-block gap
  of at least the band gap (A2, positive direction).
* `overlapping_bands_example`: the converse fails (both in-block gaps 1, bands overlap).
-/

open Matrix

namespace MassGap.Bloch

section General

variable {K : Type*} [Field K] {m n : Type*} [Fintype m] [Fintype n]
  [DecidableEq m] [DecidableEq n]

theorem algebraMap_sub_fromBlocks (a : K) (A : Matrix m m K) (D : Matrix n n K) :
    algebraMap K (Matrix (m ⊕ n) (m ⊕ n) K) a - fromBlocks A 0 0 D =
      fromBlocks (algebraMap K (Matrix m m K) a - A) 0 0
        (algebraMap K (Matrix n n K) a - D) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    simp [algebraMap_eq_diagonal, diagonal_apply]

/-- A1 (finite dimension): the spectrum of a block-diagonal matrix is the union of the
spectra of its blocks. -/
theorem spectrum_fromBlocks_diag (A : Matrix m m K) (D : Matrix n n K) :
    spectrum K (fromBlocks A 0 0 D) = spectrum K A ∪ spectrum K D := by
  ext a
  simp only [spectrum.mem_iff, Set.mem_union]
  rw [algebraMap_sub_fromBlocks, isUnit_iff_isUnit_det, det_fromBlocks_zero₂₁,
    IsUnit.mul_iff, ← isUnit_iff_isUnit_det, ← isUnit_iff_isUnit_det]
  tauto

end General

/-- `IsGroundGap S e Δ`: `e` is the lowest point of the (finite) spectrum `S`, and the next
distinct level is `e + Δ`, with `Δ > 0`. Levels are counted without multiplicity. -/
def IsGroundGap (S : Set ℝ) (e Δ : ℝ) : Prop :=
  0 < Δ ∧ e ∈ S ∧ e + Δ ∈ S ∧ ∀ x ∈ S, x = e ∨ e + Δ ≤ x

/-- A 2x2 real symmetric (diagonal) block with ground energy `e` and gap `δ`. -/
def blk (e δ : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := diagonal ![e, e + δ]

theorem blk_isHermitian (e δ : ℝ) : (blk e δ).IsHermitian := by
  unfold blk
  exact isHermitian_diagonal _

theorem spectrum_blk (e δ : ℝ) : spectrum ℝ (blk e δ) = {e, e + δ} := by
  unfold blk
  rw [spectrum_diagonal]
  ext x
  simp only [Set.mem_range, Fin.exists_fin_two, Matrix.cons_val_zero, Matrix.cons_val_one,
    Set.mem_insert_iff, Set.mem_singleton_iff]
  constructor <;> rintro (h | h) <;> simp_all

theorem isGroundGap_blk {e δ : ℝ} (hδ : 0 < δ) : IsGroundGap (spectrum ℝ (blk e δ)) e δ := by
  rw [spectrum_blk]
  refine ⟨hδ, by simp, by simp, ?_⟩
  intro x hx
  rcases hx with rfl | rfl
  · exact Or.inl rfl
  · exact Or.inr le_rfl

/-- The two-fiber Hamiltonian `H(δ) = blk e δ ⊕ blk e 1`. -/
def twoBlock (e δ : ℝ) : Matrix (Fin 2 ⊕ Fin 2) (Fin 2 ⊕ Fin 2) ℝ :=
  fromBlocks (blk e δ) 0 0 (blk e 1)

theorem twoBlock_isHermitian (e δ : ℝ) : (twoBlock e δ).IsHermitian := by
  unfold twoBlock
  exact IsHermitian.fromBlocks (blk_isHermitian e δ) (by simp) (blk_isHermitian e 1)

theorem spectrum_twoBlock (e δ : ℝ) :
    spectrum ℝ (twoBlock e δ) = {e, e + δ} ∪ {e, e + 1} := by
  unfold twoBlock
  rw [spectrum_fromBlocks_diag, spectrum_blk, spectrum_blk]

/-- Counterexample: both blocks have the same ground energy `e` (identical band data), block 2
has gap 1, block 1 has gap `δ`, and `δ > 0` is arbitrary. -/
theorem twoBlock_counterexample (e : ℝ) {δ : ℝ} (hδ : 0 < δ) :
    (twoBlock e δ).IsHermitian ∧
    IsGroundGap (spectrum ℝ (blk e δ)) e δ ∧ IsGroundGap (spectrum ℝ (blk e 1)) e 1 :=
  ⟨twoBlock_isHermitian e δ, isGroundGap_blk hδ, isGroundGap_blk one_pos⟩

/-- C4 "band data do not bound the in-sector gap": for any candidate lower bound `f` on the
in-block gap that is a positive function of the two blocks' ground energies, and any band data
`(e₁, e₂)`, there are real symmetric 2x2 blocks with exactly these ground energies whose
block-1 gap is below `f e₁ e₂`. -/
theorem band_data_do_not_bound_in_sector_gap (f : ℝ → ℝ → ℝ) (hf : ∀ e₁ e₂, 0 < f e₁ e₂)
    (e₁ e₂ : ℝ) :
    ∃ H₁ H₂ : Matrix (Fin 2) (Fin 2) ℝ, H₁.IsHermitian ∧ H₂.IsHermitian ∧
      ∃ Δ₁ Δ₂ : ℝ, IsGroundGap (spectrum ℝ H₁) e₁ Δ₁ ∧ IsGroundGap (spectrum ℝ H₂) e₂ Δ₂ ∧
        Δ₁ < f e₁ e₂ := by
  have h := hf e₁ e₂
  exact ⟨blk e₁ (f e₁ e₂ / 2), blk e₂ 1, blk_isHermitian _ _, blk_isHermitian _ _,
    f e₁ e₂ / 2, 1, isGroundGap_blk (by linarith), isGroundGap_blk one_pos, by linarith⟩

/-- A2, positive direction: if the spectrum of `A ⊕ D` has a gap `(b₀, b₁)` above the lowest
band, and block `A` has exactly one level `e` in the lowest band, then every other level of `A`
lies at least `b₁ - b₀` above `e`. -/
theorem in_block_gap_of_band_gap {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m]
    [DecidableEq n] {A : Matrix m m ℝ} {D : Matrix n n ℝ} {b₀ b₁ e : ℝ}
    (hband : ∀ x ∈ spectrum ℝ (fromBlocks A 0 0 D), x ≤ b₀ ∨ b₁ ≤ x)
    (he : e ≤ b₀) (hlow : ∀ x ∈ spectrum ℝ A, x ≤ b₀ → x = e) :
    ∀ x ∈ spectrum ℝ A, x ≠ e → b₁ - b₀ ≤ x - e := by
  intro x hx hne
  have hxH : x ∈ spectrum ℝ (fromBlocks A 0 0 D) := by
    rw [spectrum_fromBlocks_diag]; exact Or.inl hx
  rcases hband x hxH with h | h
  · exact absurd (hlow x hx h) hne
  · linarith

/-- The converse of A2 fails: both blocks have in-block gap 1, yet the bands overlap
(the higher ground energy 2 exceeds the lower first excited level 1). -/
theorem overlapping_bands_example :
    IsGroundGap (spectrum ℝ (blk 0 1)) 0 1 ∧ IsGroundGap (spectrum ℝ (blk 2 1)) 2 1 ∧
      spectrum ℝ (fromBlocks (blk 0 1) 0 0 (blk 2 1)) = {0, 1} ∪ {2, 3} := by
  refine ⟨isGroundGap_blk one_pos, isGroundGap_blk one_pos, ?_⟩
  rw [spectrum_fromBlocks_diag, spectrum_blk, spectrum_blk]
  norm_num

end MassGap.Bloch
