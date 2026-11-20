# 20261120-009: skeptic audit of C4 (theta/Bloch) and C7 (Payne-Weinberger on the Gribov region)

Audits `work/20261120-008-theta-obstruction-gribov-scope/` (README, `theta_bands.py`, `gribov_scope.py`). Tags: **[P]** proved here, **[C]** cited, **[N]** numerical observation, **[J]** conjectured.

Files:
- `rerun/` unmodified copies of the 008 scripts (copied so the rerun does not overwrite 008 outputs). `rerun/theta_summary.txt` is the rerun output.
- `skeptic_checks.py` -> `skeptic_summary.txt`, `run.log` (400 s, 1 thread, < 300 MB). Imports the model code from `rerun/`.
- Rerun: `python tools/pin.py .venv/Scripts/python.exe work/20261120-009-skeptic-theta-gribov/skeptic_checks.py`.

## Question
Do the C4 and C7 obstructions hold as stated, is every "proved" step proved, and does either claim overstate what is blocked?

## Method
Reread every step of 008 against its setting; reran `theta_bands.py` (output byte-identical to 008's `theta_summary.txt`); built explicit counterexamples or checks for each target; extended the d = 3 Gribov diameter scan from L = 3, 4 to L = 3..6 with more restarts; recomputed the PW/m_0++ ratios from measured diameters only, with m/sqrt(sigma) read from `data/` (AT2021 SU(2) 0++ = 3.781, same row 008 used).

---

## C4 (theta/Bloch): survives, with two wording fixes

**A1 (no in-sector bound from Bloch decomposition).** Correct: the direct-integral decomposition is a block diagonalization and spec H = closure of union of spec H_θ [C, Reed-Simon IV XIII.16]. [P] holds.

**A2 (Δ_ind ≤ Δ(θ), converse fails).** Correct and elementary: E_1(θ) − E_0(θ) ≥ min E_1 − max E_0. [P] holds.

**A3 (pendulum numbers).** Rerun identical. Spot checks: instanton action 8 sqrt(V0) = ∫ sqrt(2V) dx is right; prefactor s^{3/2} matches ω sqrt(S); χ = W_0/2 is the tight-binding relation; Bragg gap 2|V_1| = V0 is right (off-diagonal −V0/2 between m = 0 and m = −1 at θ = π). `chi_exact` second-order formula is right (H'' = 1/4π² times identity).

**A4 (N-well product model): wording wrong, math right.** 008 says it imposes T_i ψ = e^{iθ} ψ for every copy and calls this "θ couples to the total winding". These are different models:
- With every T_i fixed (008's computation): E_0 = N e_0(θ), Δ^{(N)} = Δ, bandwidth N W_0. Correct [P].
- With only T_total = T_1...T_N fixed (what "θ couples to total winding" means, and the only large-gauge symmetry YM has on R^3, π_3 = Z): the fiber still carries the commuting T_i, so it decomposes over a continuous angle φ and its spectrum is ∪_φ {e_a(φ) + e_b(θ − φ)}. The lowest band is an interval, no isolated ground state, in-sector gap **0** already at N = 2 [P; N in `skeptic_summary.txt` C4-a: V0 = 1, θ = 0, band [0.92987, 0.93520]].
So the fixed-all-T_i model is the right stand-in for "θ as a local Lagrangian parameter, every local tunneling event carries e^{iθ}" (dilute gas), and the extensivity W ∝ N is correct there. It is not a Bloch decomposition along a single symmetry. The T_total version actually strengthens A1 (Bloch along one generator can leave the fiber gapless), but the README must say which model it uses.

**A5 (2D examples).** U(1) on a circle: E_n(θ) = (g²L/2)(n + θ/2π)², χ = g²/4π², crossing at θ = π: all correct [P]. Massless Schwinger: χ = E_0'' = 0 by chiral rotation, gap e/sqrt(π): correct [C: Schwinger, Phys. Rev. 128, 2425 (1962); Coleman-Jackiw-Susskind, Ann. Phys. 93, 267 (1975); Coleman, Ann. Phys. 101, 239 (1976); citations match my knowledge, not re-fetched]. Limitation: both examples are gapped above the vacuum (U(1) trivially, with gap g²L/2 → ∞), so they show "χ = 0 does not preclude a gap" and "χ > 0 does not produce a particle", but not "χ > 0 with a gapless spectrum". The independence claim should be worded as the first two. (U(1) at θ = π does give Δ(π) = 0 by vacuum degeneracy, which 008 mentions.)

**Could a cleverer Bloch-type argument bound the in-sector gap?** I checked the two obvious candidates:
1. θ-curvature sum rule [P, new]: second-order perturbation theory in θ gives E_0'' = ⟨H''⟩ − 2 Σ_n |⟨n|H'|0⟩|²/(E_n − E_0) ≥ ⟨H''⟩ − 2 Var_0(H')/Δ, hence **Δ(θ) ≤ 2 Var_0(H')/(⟨H''⟩ − E_0''(θ))**. This is an upper bound on the gap, never a lower one (Cauchy-Schwarz on the same sum also yields only lower bounds on the inverse moment). Checked on the pendulum (C4-b): V0 = 4, Δ = 1.86557 ≤ 1.86830; all 8 cases hold, near-saturated (single-mode dominance).
2. 't Hooft electric-flux / center sectors and twisted boundary conditions: these are again a block diagonalization ([C] 't Hooft NPB 153, 141 (1979)); they control sector ground energies E(e) ~ σL (confinement data), and twist/flux-insertion arguments of LSM/Oshikawa type are variational, giving upper bounds on gaps. The 0++ gap lives inside the e = 0, m = 0 sector, where none of these says anything. [P for the structure; no counterexample to the obstruction found.]
So "blocked" is fair for Bloch decomposition plus band-gap/θ-data arguments. It should not be stated as "every symmetry-sector method is blocked": a gap theorem could still use sector structure as one ingredient.

**Verdict C4: survives (weakened wording only).**

---

## C7 (Payne-Weinberger on the Gribov region): weakened

**(i) Convexity.** q_A(ω) affine in A, so Ω (closure) is an intersection of half-spaces: convex [P, correct]. Continuum: [C] Zwanziger NPB 209, 336 (1982); Dell'Antonio-Zwanziger CMP 138, 291 (1991). Λ convex: ||A^g||² − ||A||² is affine in A for fixed g, so Λ is an intersection of half-spaces [P, one line, can replace the citation]. Compact Wilson lattice: Ω lives in a product of 3-spheres, Euclidean convexity meaningless; 008 says so correctly.

**Boundedness.** K = M_0^{-1/2} M_1 M_0^{-1/2} on constants-orthogonal ghosts is traceless (colour blocks antisymmetric), so every ray is bounded unless K = 0 [P]; in finite dimension a closed convex set with every ray bounded is bounded [P, recession cone]. The missing piece is still "K(A) ≠ 0 for every nonzero transverse A, L ≥ 3" (L = 2 fails, as 008 found). Open, as 008 states.

**Normalization.** 008 calls the coefficient of the A term in q_A "a convention". It is not: rescaling it rescales Ω and the PW quantity by its square. The coefficient 1 used is the one consistent with Vmag's F = dA + A × A and kinetic term (g²/2) Σ E^a E^a (checked: lattice discretisation of ∂ω · (A × ω)). Drop the word "convention".

**(ii) log det M concave on Ω°.** f'' = −tr[(M^{-1/2} M_1 M^{-1/2})²] ≤ 0: correct [P].

**(iii) Min-max lemma on Λ.** Correct for a flat weighted Laplacian: the θ = 0 quotient form domain (H¹ functions with matching traces on identified boundary faces) is a subspace of H¹(Λ, μ) with the same form and inner product, so λ_k(quotient) ≥ λ_k^Neumann; constants are admissible in both, so λ_1 ≥ λ_1^N ≥ π²/diam(Λ)² [C: FNT 1207.0680 for log-concave μ]. Boundary identifications enter only through the trace constraint, so the lemma does not need J to agree on identified faces. θ ≠ 0 caveat is right. [P holds, conditional on the flat-Laplacian hypothesis, which 008 states is not the actual operator.]

**Direction error: Λ vs Ω.** The lemma uses diam(Λ); all lattice numbers are diam(Ω) estimates. Since Λ ⊂ Ω, diam(Λ) ≤ diam(Ω) and π²/diam(Λ)² ≥ π²/diam(Ω)², so the lattice ratios are **not** upper bounds on the Λ quantity the lemma produces. On the noncompact lattice Λ is not even defined (008 says so), so the lattice numbers describe a Gribov-Zwanziger-style restriction to Ω, not the quotient problem. The only correctly directed statement for Λ is 008's continuum T³ chord: diam(Λ) ≥ 2π sqrt(3L) [P, normalisation checked], giving Δ_PW ≤ g²/(24L), i.e. Δ_PW/m_0++ ≤ g²/(24 m L) = 3.7e-3 g² at L sqrt(σ) = 3. This is weaker than "7e-4" but is an actual bound.

**(iv)1 "Gap not monotone in added V ≥ 0".** The statement H ≥ K does not imply gap(H) ≥ gap(K) is right, and the conclusion (V cannot be dropped) is right. The example given is wrong: a free particle on R has no gap to close. Correct example (C7-a): Neumann −d²/dx² on [−1, 1], gap π²/4 = 2.467; add V = 2000·1{|x| < 0.05} ≥ 0, gap 0.0055 (lower); add V = 50x² ≥ 0, gap 13.69 (higher). The ground-state-transform route and its log-concavity hypothesis are stated correctly.

**(iv)2 Magnetic non-convexity: real, and stronger than claimed.** 008's evidence is at L = 2 (the artifact lattice) at |A| = 9.6-16.4. A cleaner [P] statement for every L, d ≥ 2, at arbitrarily small |A|: constant fields are transverse and the origin is interior to Ω; on constant fields Vmag = (V/2)|A_1 × A_2|²; for A_1 = x e_1, A_2 = x e_2 the restricted Hessian has eigenvalues 3x², −x² (per-link component x), and a negative eigenvalue of a principal restriction is one of the full Hessian. Checked numerically on the full transverse Hessian at L = 3, d = 3 (C7-b): −3.014e-2 at |A| = 1.27 (analytic −3.000e-2), with the chord in that direction extending to ±25.5, so the point is deep inside Ω. Replace the L = 2 numbers with this.

**(i)/(iv)3 Diameter scaling: claim wrong; conclusion survives.** 008 claims D² ∝ V in d = 2 and "between L² and L³" in d = 3, and extrapolates with D² = 163 V. Its own d = 2 table has D²/V falling (112 → 41, local exponent of D ≈ 0.8 at L = 5 → 8), and its d = 3 pair gives D ∝ L^{0.83}, i.e. D² ∝ L^{1.66}, outside the stated range. Extended d = 3 scan (C7-c, gradient ascent, 12/8/4/3 restarts):

| L | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| D_lower | 72.8 | 91.1 | 120.4 | 150.4 |
| D_lower/L | 24.3 | 22.8 | 24.1 | 25.1 |
| D²/V | 196 | 130 | 116 | 105 |

Fit D_lower ~ L^{1.05}: D ∝ L, D² ∝ L², not V [N, 4 sizes, lower bounds]. Consequence: 008's N_s = 8-32 ratios used D² = 163 N_s³, overstating D and understating the ratio by ~10x at N_s = 32 (2.6e-6 at g² = 1 vs ~1.8e-5 with D ∝ L). These were extrapolations, not "lattices tried", and not bounds.

**"≤ 7e-4 m_0++ on all lattices tried".** False as stated. With measured D (C7-d, L sqrt(σ) = 3): N_s = 3, g² = 4: 9.8e-4; N_s = 4, g² = 4: 8.4e-4 (008 used the averaged c = 163 at N_s = 4, which is not the measured L = 4 value). Units check fine: a·Δ_PW = (g²/2)π²/D² with lattice A, a·m = (m/sqrt σ)(L sqrt σ)/N_s. But the g² values are not tied to the lattice spacing (g² = 4 at a sqrt(σ) = 0.75 is not on the SU(2) scaling trajectory, roughly g² ≈ 2 there), and the comparison mixes a noncompact Hamiltonian lattice with a compact Euclidean one. Treat these as orders of magnitude.

**Limits.** With D ∝ L: Δ_PW/m_0++ ∝ g²(a)/N_s at fixed physical L → 0 as a → 0; at fixed a, ∝ 1/L → 0 as L → ∞ (relative to fixed a·m). Continuum Λ bound g²/(24L) → 0 relative to m as L → ∞ and, at fixed L, as g² → 0. So "→ 0 in both limits" survives, as an [N] statement on Ω plus the [P] continuum chord for Λ.

**Verdict C7: weakened.** The obstruction stands (failure point = unproved log-concavity of J ψ_0², plus kinetic bound tending to 0), but the numerics are misstated and partly misdirected.

---

## Required wording changes to CLAIMS.md

C4: replace "band width grows with volume while the in-sector gap does not (N-well product model)" by "in an N-well model where every local tunneling carries e^{iθ} (dilute-gas structure), bandwidth ∝ N while the fiber gap is N-independent; Bloch decomposition along the total winding alone leaves a gapless fiber already at N = 2". Replace "show band data and mass gap are independent" by "show χ = 0 is compatible with a gap and χ > 0 with no particle states". Add: "θ-curvature sum rules give only upper bounds Δ ≤ 2Var(H')/(⟨H''⟩ − E_0'')". Status: proved (abstract + toy models), skeptic passed with rewording. Scope: blocks Bloch decomposition plus band-gap/θ-data arguments, not every use of sector structure.

C7: replace "magnetic potential is non-convex inside the region (Hessian down to -11.3)" by "magnetic potential is non-convex at arbitrarily small |A| inside Ω on every lattice (constant-field Hessian eigenvalue −x², proved)". Replace "Bound <= 7e-4 m_0++ on all lattices tried and -> 0" by "Ω-diameter estimate D ∝ L (d = 3, L = 3..6, lower bounds), kinetic PW quantity ≤ 1e-3 m_0++ on measured lattices (orders of magnitude; Ω not Λ, noncompact discretisation); continuum T³: Δ_PW(Λ) ≤ g²/(24L); → 0 as a → 0 and L → ∞". Status: proved (convexity, concavity, quotient lemma for flat weighted Laplacian, non-convexity), numerical (diameters, 4 sizes in d = 3), skeptic passed with rewording.

008 README errata (for the overseer): duplicated sentence in (i); "free particle" example in (iv)1; "convention" for the FP coefficient; D² ∝ V claims.

## Open issues
- Boundedness of Ω for all L ≥ 3 (K(A) ≠ 0 for transverse A ≠ 0) still unproved.
- No lower bound on diam(Λ) for any lattice discretisation; compact Wilson lattice FMR (Zwanziger NPB 412, 657 (1994)) still unread.
- D_lower at L = 6 used only 3 restarts; the D ∝ L fit could shift upward (only strengthens the verdict).
- Example of χ > 0 with gapless spectrum not supplied (not needed for the obstruction).
- Citations in A5/A6 checked against my own knowledge only, not re-fetched.
