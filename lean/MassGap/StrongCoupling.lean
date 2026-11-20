import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Basic.Real.Basic

/-!
# C2: second-order strong-coupling coefficient 56/15 (SU(2) KS, single plaquette)

Setting: SU(2) Kogut-Susskind Hamiltonian on one plaquette, gauge-invariant sector, lattice,
finite volume: `H = (g²/2) Σ_links E² - (1/g²)(Tr U_p + Tr U_p†)`. The states are `|j⟩`,
`j = 0, 1/2, 1, ...`, with unperturbed energy `2 g² j (j+1)` (four links, Casimir `j(j+1)`),
and the perturbation `V = -(2/g²) Tr U_p` has `⟨j ± 1/2 | V | j⟩ = -2/g²` (Clebsch-Gordan
`χ_{1/2} χ_j = χ_{j-1/2} + χ_{j+1/2}`). These inputs are taken from work 003/005 and are not
formalized here.

What is checked: given these inputs, the Rayleigh-Schrodinger second-order shifts
`δE_j = Σ_{k ≠ j} |V_{kj}|² / (E_j - E_k)` of the ground state (`j = 0`) and of the first
excited state (`j = 1/2`) give `Δ⁽²⁾ = δE_{1/2} - δE_0 = 56 / (15 g⁶)`. Here `x = g²`.
-/

namespace MassGap.StrongCoupling

/-- Unperturbed single-plaquette energy at spin `j`, `x = g²`. -/
noncomputable def e0 (x j : ℝ) : ℝ := 2 * x * j * (j + 1)

/-- Off-diagonal matrix element `⟨j ± 1/2|V|j⟩`. -/
noncomputable def vEl (x : ℝ) : ℝ := -2 / x

/-- Leading order: `E_{1/2} - E_0 = 3 g² / 2`. -/
theorem firstOrder_gap (x : ℝ) : e0 x (1/2) - e0 x 0 = 3 * x / 2 := by
  unfold e0; ring

theorem secondOrder_ground (x : ℝ) (hx : x ≠ 0) :
    vEl x ^ 2 / (e0 x 0 - e0 x (1/2)) = -8 / (3 * x ^ 3) := by
  unfold vEl e0
  field_simp
  ring

theorem secondOrder_first (x : ℝ) (hx : x ≠ 0) :
    vEl x ^ 2 / (e0 x (1/2) - e0 x 0) + vEl x ^ 2 / (e0 x (1/2) - e0 x 1) =
      16 / (15 * x ^ 3) := by
  unfold vEl e0
  field_simp
  ring

/-- C2 second-order coefficient: `Δ⁽²⁾ = 56 / (15 g⁶)`. -/
theorem secondOrder_gap (x : ℝ) (hx : x ≠ 0) :
    (vEl x ^ 2 / (e0 x (1/2) - e0 x 0) + vEl x ^ 2 / (e0 x (1/2) - e0 x 1)) -
      vEl x ^ 2 / (e0 x 0 - e0 x (1/2)) = 56 / (15 * x ^ 3) := by
  rw [secondOrder_first x hx, secondOrder_ground x hx]
  field_simp
  ring

end MassGap.StrongCoupling
