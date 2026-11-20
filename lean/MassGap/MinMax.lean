import Mathlib.Analysis.InnerProductSpace.Spectrum
import Mathlib.Analysis.InnerProductSpace.Projection.FiniteDimensional
import Mathlib.Order.Interval.Finset.Fin
import MassGap.SchurBound

/-!
# Easy half of Courant-Fischer, and the k-eigenvalue Schur/Feshbach count (C5)

Setting: finite-dimensional inner product space `M` over `𝕜 = ℝ` or `ℂ` (any `RCLike`),
`n = finrank M`, `T` symmetric (self-adjoint). Eigenvalues are mathlib's
`LinearMap.IsSymmetric.eigenvalues`, a list `Fin n → ℝ` with multiplicity, sorted decreasing.

* `card_lt_add_finrank_le` (min-max, easy half): if `⟪x, T x⟫ ≥ E ‖x‖²` on a subspace `V`, then
  `#{i | λᵢ < E} + dim V ≤ n`.
* `card_eigenvalues_lt_le`: hence if `dim V ≥ n - k`, `T` has at most `k` eigenvalues `< E`
  (counted with multiplicity).
* `le_eigenvalue_of_card_lt_le`: in sorted form, the `(k+1)`-th lowest eigenvalue is `≥ E`.
* `schur_count` (C5, k-eigenvalue version): `H` symmetric on `M`, Q-space `K`, P-space `Kᗮ`,
  `V ≤ Kᗮ` of codimension `k` in `Kᗮ`. Floor `d` on `K`, `E < d`, coupling bound `X ≥ B Bᴴ`
  on `V`, and `A - X/(d - E) ≥ E` on `V`. Then `H` has at most `k` eigenvalues `< E`.
  Built on `MassGap.Schur.feshbach_form_bound`.
* `schur_count_sorted`: same, as `E ≤ λ_(k+1)` (k = 1: `E` is a lower bound on `E₁`).
-/

open Module RCLike
open scoped InnerProductSpace

namespace MassGap.MinMax

variable {𝕜 : Type*} [RCLike 𝕜] {M : Type*} [NormedAddCommGroup M] [InnerProductSpace 𝕜 M]
  [FiniteDimensional 𝕜 M]

section MinMax

variable {T : M →ₗ[𝕜] M} (hT : T.IsSymmetric) {n : ℕ} (hn : finrank 𝕜 M = n)

/-- The span of the eigenvectors with eigenvalue `< E`. -/
noncomputable def lowSpan (E : ℝ) : Submodule 𝕜 M :=
  Submodule.span 𝕜
    (Set.range fun i : {i : Fin n // hT.eigenvalues hn i < E} => hT.eigenvectorBasis hn i.1)

theorem finrank_lowSpan (E : ℝ) :
    finrank 𝕜 (lowSpan hT hn E) = Fintype.card {i : Fin n // hT.eigenvalues hn i < E} := by
  have hon : Orthonormal 𝕜
      (fun i : {i : Fin n // hT.eigenvalues hn i < E} => hT.eigenvectorBasis hn i.1) :=
    (hT.eigenvectorBasis hn).orthonormal.comp _ Subtype.val_injective
  exact finrank_span_eq_card hon.linearIndependent

/-- On `lowSpan`, the form `⟪x, T x⟫ ≥ E ‖x‖²` forces `x = 0`. -/
theorem eq_zero_of_mem_lowSpan (E : ℝ) {x : M} (hx : x ∈ lowSpan hT hn E)
    (hle : E * ‖x‖ ^ 2 ≤ re ⟪x, T x⟫_𝕜) : x = 0 := by
  set S := {i : Fin n // hT.eigenvalues hn i < E}
  set b := hT.eigenvectorBasis hn
  set μ := hT.eigenvalues hn
  have hon : Orthonormal 𝕜 (fun i : S => b i.1) := b.orthonormal.comp _ Subtype.val_injective
  obtain ⟨c, rfl⟩ := (Submodule.mem_span_range_iff_exists_fun 𝕜).1 hx
  have hTx : T (∑ i, c i • b i.1) = ∑ i : S, (c i * (μ i.1 : 𝕜)) • b i.1 := by
    simp only [map_sum, map_smul, b, μ, LinearMap.IsSymmetric.apply_eigenvectorBasis, smul_smul]
    rfl
  have hnorm : ‖∑ i, c i • b i.1‖ ^ 2 = ∑ i : S, ‖c i‖ ^ 2 := by
    rw [← inner_self_eq_norm_sq (𝕜 := 𝕜), hon.inner_sum c c Finset.univ, map_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [RCLike.conj_mul]; norm_cast
  have hform : re ⟪∑ i, c i • b i.1, T (∑ i, c i • b i.1)⟫_𝕜 = ∑ i : S, μ i.1 * ‖c i‖ ^ 2 := by
    rw [hTx, hon.inner_sum, map_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← mul_assoc, RCLike.conj_mul]
    norm_cast
    ring
  rw [hnorm, hform, Finset.mul_sum] at hle
  have hsum : ∑ i : S, (E - μ i.1) * ‖c i‖ ^ 2 ≤ 0 := by
    have : ∑ i : S, (E - μ i.1) * ‖c i‖ ^ 2 =
        ∑ i : S, E * ‖c i‖ ^ 2 - ∑ i : S, μ i.1 * ‖c i‖ ^ 2 := by
      rw [← Finset.sum_sub_distrib]; exact Finset.sum_congr rfl fun i _ => by ring
    linarith
  have hnn : ∀ i ∈ (Finset.univ : Finset S), 0 ≤ (E - μ i.1) * ‖c i‖ ^ 2 := fun i _ =>
    mul_nonneg (by linarith [i.2]) (sq_nonneg _)
  have hzero := (Finset.sum_eq_zero_iff_of_nonneg hnn).1 (le_antisymm hsum
    (Finset.sum_nonneg hnn))
  have hc : ∀ i, c i = 0 := by
    intro i
    have h := hzero i (Finset.mem_univ _)
    have hpos : 0 < E - μ i.1 := by linarith [i.2]
    have : ‖c i‖ ^ 2 = 0 := by
      rcases mul_eq_zero.1 h with h | h
      · linarith
      · exact h
    simpa using this
  simp [hc]

/-- **Courant-Fischer, easy half.** If the form of `T` is `≥ E` on `V`, the number of
eigenvalues `< E` (with multiplicity) plus `dim V` is at most `n`. -/
theorem card_lt_add_finrank_le (V : Submodule 𝕜 M) (E : ℝ)
    (hV : ∀ x ∈ V, E * ‖x‖ ^ 2 ≤ re ⟪x, T x⟫_𝕜) :
    Fintype.card {i : Fin n // hT.eigenvalues hn i < E} + finrank 𝕜 V ≤ n := by
  have hdisj : Disjoint (lowSpan hT hn E) V := Submodule.disjoint_def.2 fun x hxW hxV =>
    eq_zero_of_mem_lowSpan hT hn E hxW (hV x hxV)
  have h := Submodule.finrank_add_finrank_le_of_disjoint hdisj
  rwa [finrank_lowSpan, hn] at h

/-- Codimension-`k` form: form `≥ E` on `V` with `dim V ≥ n - k` implies at most `k`
eigenvalues `< E`. -/
theorem card_eigenvalues_lt_le (V : Submodule 𝕜 M) (E : ℝ) (k : ℕ)
    (hcodim : n ≤ finrank 𝕜 V + k) (hV : ∀ x ∈ V, E * ‖x‖ ^ 2 ≤ re ⟪x, T x⟫_𝕜) :
    (Finset.univ.filter fun i => hT.eigenvalues hn i < E).card ≤ k := by
  have h := card_lt_add_finrank_le hT hn V E hV
  rw [Fintype.card_subtype] at h
  omega

/-- Sorted form: if at most `k < n` eigenvalues are `< E`, then the `(k+1)`-th lowest
eigenvalue, `eigenvalues ⟨n - 1 - k⟩` in mathlib's decreasing order, is `≥ E`. -/
theorem le_eigenvalue_of_card_lt_le (E : ℝ) (k : ℕ) (hk : k < n)
    (hcard : (Finset.univ.filter fun i => hT.eigenvalues hn i < E).card ≤ k) :
    E ≤ hT.eigenvalues hn ⟨n - 1 - k, by omega⟩ := by
  by_contra hlt
  push Not at hlt
  set m : Fin n := ⟨n - 1 - k, by omega⟩
  have hsub : Finset.Ici m ⊆ Finset.univ.filter fun i => hT.eigenvalues hn i < E := by
    intro j hj
    simp only [Finset.mem_Ici] at hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact lt_of_le_of_lt (hT.eigenvalues_antitone hn hj) hlt
  have h1 := Finset.card_le_card hsub
  rw [Fin.card_Ici] at h1
  simp only [m] at h1
  omega

end MinMax

section SchurCount

variable {H X : M →ₗ[𝕜] M}

/-- The form `⟪w, H w⟫ ≥ E ‖w‖²` on `V ⊔ K`, from `feshbach_form_bound`. -/
theorem form_ge_on_sup (hH : H.IsSymmetric) (K V : Submodule 𝕜 M) (hVK : V ≤ Kᗮ) {d E : ℝ}
    (hD : ∀ v ∈ K, d * ‖v‖ ^ 2 ≤ re ⟪v, H v⟫_𝕜)
    (hX : ∀ u ∈ V, ∀ v ∈ K, ‖⟪u, H v⟫_𝕜‖ ^ 2 ≤ re ⟪u, X u⟫_𝕜 * ‖v‖ ^ 2)
    (hX0 : ∀ u ∈ V, 0 ≤ re ⟪u, X u⟫_𝕜)
    (hE : E < d)
    (hS : ∀ u ∈ V, E * ‖u‖ ^ 2 ≤ re ⟪u, H u⟫_𝕜 - re ⟪u, X u⟫_𝕜 / (d - E)) :
    ∀ w ∈ V ⊔ K, E * ‖w‖ ^ 2 ≤ re ⟪w, H w⟫_𝕜 := by
  -- Blocks of `H` (and `X`) as maps between the subspaces `V` (P side) and `K` (Q side).
  let A : V →ₗ[𝕜] V := V.orthogonalProjectionOnto.toLinearMap ∘ₗ H ∘ₗ V.subtype
  let X' : V →ₗ[𝕜] V := V.orthogonalProjectionOnto.toLinearMap ∘ₗ X ∘ₗ V.subtype
  let B : K →ₗ[𝕜] V := V.orthogonalProjectionOnto.toLinearMap ∘ₗ H ∘ₗ K.subtype
  let D : K →ₗ[𝕜] K := K.orthogonalProjectionOnto.toLinearMap ∘ₗ H ∘ₗ K.subtype
  have eA : ∀ u : V, ⟪u, A u⟫_𝕜 = ⟪(u : M), H u⟫_𝕜 := fun u => by
    simp only [A, LinearMap.coe_comp, Function.comp_apply, Submodule.coe_subtype,
      ContinuousLinearMap.coe_coe, Submodule.inner_orthogonalProjectionOnto_eq_of_mem_left]
  have eX : ∀ u : V, ⟪u, X' u⟫_𝕜 = ⟪(u : M), X u⟫_𝕜 := fun u => by
    simp only [X', LinearMap.coe_comp, Function.comp_apply, Submodule.coe_subtype,
      ContinuousLinearMap.coe_coe, Submodule.inner_orthogonalProjectionOnto_eq_of_mem_left]
  have eB : ∀ (u : V) (v : K), ⟪u, B v⟫_𝕜 = ⟪(u : M), H v⟫_𝕜 := fun u v => by
    simp only [B, LinearMap.coe_comp, Function.comp_apply, Submodule.coe_subtype,
      ContinuousLinearMap.coe_coe, Submodule.inner_orthogonalProjectionOnto_eq_of_mem_left]
  have eD : ∀ v : K, ⟪v, D v⟫_𝕜 = ⟪(v : M), H v⟫_𝕜 := fun v => by
    simp only [D, LinearMap.coe_comp, Function.comp_apply, Submodule.coe_subtype,
      ContinuousLinearMap.coe_coe, Submodule.inner_orthogonalProjectionOnto_eq_of_mem_left]
  have key := MassGap.Schur.feshbach_form_bound A X' B D ⊤ (d := d) (E := E)
    (fun v => by rw [eD]; simpa using hD v v.2)
    (fun u v => by rw [eB, eX]; simpa using hX u u.2 v v.2)
    (fun u => by rw [eX]; exact hX0 u u.2)
    hE
    (fun u _ => by rw [eA, eX]; simpa using hS u u.2)
  intro w hw
  obtain ⟨u, hu, v, hv, rfl⟩ := Submodule.mem_sup.1 hw
  have h := key ⟨u, hu⟩ trivial ⟨v, hv⟩
  rw [eA, eB, eD] at h
  have horth : ⟪u, v⟫_𝕜 = 0 := Submodule.inner_left_of_mem_orthogonal hv (hVK hu)
  have hnorm : ‖u + v‖ ^ 2 = ‖u‖ ^ 2 + ‖v‖ ^ 2 := by
    rw [← inner_self_eq_norm_sq (𝕜 := 𝕜), ← inner_self_eq_norm_sq (𝕜 := 𝕜),
      ← inner_self_eq_norm_sq (𝕜 := 𝕜), inner_add_left, inner_add_right, inner_add_right,
      horth, ← inner_conj_symm v u, horth]
    simp
  have hcross : re ⟪v, H u⟫_𝕜 = re ⟪u, H v⟫_𝕜 := by
    rw [← hH v u, ← inner_conj_symm, RCLike.conj_re]
  have hexp : re ⟪u + v, H (u + v)⟫_𝕜 =
      re ⟪u, H u⟫_𝕜 + 2 * re ⟪u, H v⟫_𝕜 + re ⟪v, H v⟫_𝕜 := by
    rw [map_add, inner_add_left, inner_add_right, inner_add_right]
    simp only [map_add]
    rw [hcross]; ring
  rw [hnorm, hexp]
  exact h

/-- **C5, k-eigenvalue Schur/Feshbach count.** Hamiltonian `H` symmetric on `M` (`dim M = n`),
Q-space `K`, P-space `Kᗮ`, `V ≤ Kᗮ` with `dim Kᗮ ≤ dim V + k` (codimension `≤ k` in the
P-space). Floor `d` on `K`, `E < d`, coupling bound `X ≥ B Bᴴ` on `V` (`B` = P-Q block of `H`),
Schur condition `A - X/(d - E) ≥ E` on `V`. Then `H` has at most `k` eigenvalues `< E`,
counted with multiplicity. -/
theorem schur_count (hH : H.IsSymmetric) {n : ℕ} (hn : finrank 𝕜 M = n) (K V : Submodule 𝕜 M) (hVK : V ≤ Kᗮ)
    (k : ℕ) (hcodim : finrank 𝕜 Kᗮ ≤ finrank 𝕜 V + k) {d E : ℝ}
    (hD : ∀ v ∈ K, d * ‖v‖ ^ 2 ≤ re ⟪v, H v⟫_𝕜)
    (hX : ∀ u ∈ V, ∀ v ∈ K, ‖⟪u, H v⟫_𝕜‖ ^ 2 ≤ re ⟪u, X u⟫_𝕜 * ‖v‖ ^ 2)
    (hX0 : ∀ u ∈ V, 0 ≤ re ⟪u, X u⟫_𝕜)
    (hE : E < d)
    (hS : ∀ u ∈ V, E * ‖u‖ ^ 2 ≤ re ⟪u, H u⟫_𝕜 - re ⟪u, X u⟫_𝕜 / (d - E)) :
    (Finset.univ.filter fun i => hH.eigenvalues hn i < E).card ≤ k := by
  have hform := form_ge_on_sup hH K V hVK hD hX hX0 hE hS
  have hdisj : Disjoint V K :=
    Disjoint.mono_left hVK (Submodule.orthogonal_disjoint K).symm
  have hdim : finrank 𝕜 ↥(V ⊔ K) = finrank 𝕜 V + finrank 𝕜 K := by
    have h := Submodule.finrank_sup_add_finrank_inf_eq V K
    rw [hdisj.eq_bot, finrank_bot, add_zero] at h
    exact h
  have hKK := Submodule.finrank_add_finrank_orthogonal K
  refine card_eigenvalues_lt_le hH hn (V ⊔ K) E k ?_ hform
  omega

/-- Sorted form of `schur_count`: for `k < n`, `E ≤ λ_(k+1)` (the `(k+1)`-th lowest eigenvalue).
With `k = 1` this is the E1 lower bound used in C5. -/
theorem schur_count_sorted (hH : H.IsSymmetric) {n : ℕ} (hn : finrank 𝕜 M = n) (K V : Submodule 𝕜 M) (hVK : V ≤ Kᗮ)
    (k : ℕ) (hk : k < n) (hcodim : finrank 𝕜 Kᗮ ≤ finrank 𝕜 V + k) {d E : ℝ}
    (hD : ∀ v ∈ K, d * ‖v‖ ^ 2 ≤ re ⟪v, H v⟫_𝕜)
    (hX : ∀ u ∈ V, ∀ v ∈ K, ‖⟪u, H v⟫_𝕜‖ ^ 2 ≤ re ⟪u, X u⟫_𝕜 * ‖v‖ ^ 2)
    (hX0 : ∀ u ∈ V, 0 ≤ re ⟪u, X u⟫_𝕜)
    (hE : E < d)
    (hS : ∀ u ∈ V, E * ‖u‖ ^ 2 ≤ re ⟪u, H u⟫_𝕜 - re ⟪u, X u⟫_𝕜 / (d - E)) :
    E ≤ hH.eigenvalues hn ⟨n - 1 - k, by omega⟩ :=
  le_eigenvalue_of_card_lt_le hH hn E k hk
    (schur_count hH hn K V hVK k hcodim hD hX hX0 hE hS)

end SchurCount

end MassGap.MinMax
