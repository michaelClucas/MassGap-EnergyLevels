# 20261120-007 Frustration-free non-Abelian gauge Hamiltonian and finite-size gap criteria

## Question

Is line of attack 2 (AKLT / Knabe / martingale / Lemm-Mozgunov spectral-gap methods) alive for non-Abelian lattice gauge Hamiltonians? If so, what is the frustration-free (FF) model, what finite-size criterion would prove its gap, and what certified numerical input does it need? What blocks the route from that model to Kogut-Susskind (KS) and to the continuum?

## Method

- Literature: INSPIRE citation lists of arXiv:2409.10612, arXiv:2204.12737 and Dell'Antonio-Zwanziger 1991; web searches for FF / RK / string-net non-Abelian gauge models and Glauber dynamics for lattice gauge theory; Yarotsky's hypotheses read from the paper (Theorem 1, pp. 2-3). Every arXiv ID was checked against the arXiv API. Details in `literature-notes.md`.
- Construction and proofs by hand (below). Each step is tagged **[proved]**, **[cited]**, **[numerical]** or **[conjectured]**.
- `checks.py` (log: `checks.log`, about 5 s, 1 thread): single-link gap of the FF term, and a random-matrix test of the two-projector inequality used in section 4.

Conventions (match `src/ks_su2.py`). SU(2), U = cos θ + i sin θ n·σ, identified with the unit sphere S³ ⊂ R⁴. E² = -Δ_{S³}/4 per link (eigenvalue j(j+1) on spin-j functions). Wilson action of the spatial lattice S(U) = Σ_p (1 - ½ tr U_p). KS in the repo normalization:

    H_KS = (g²/2) [ Σ_ℓ E_ℓ² + x Σ_p (1 - ½ tr U_p) ] + const,   x = 8/g⁴.

Setting throughout unless stated: Hamiltonian lattice, d spatial dimensions (d = 2 or 3), finite box or torus Λ_L ⊂ Z^d, gauge group SU(2), fixed lattice spacing. All bounds are uniform in L where stated.

## Result

### 1. The FF Hamiltonian (RK / stochastic-quantization type)

For a link ℓ let ∇_ℓ be the gradient on the S³ factor of ℓ. Define ψ_β = exp(-β S/2) and

    Q_ℓ = ½ (∇_ℓ + (β/2) ∇_ℓ S) = ½ ψ_β ∇_ℓ ψ_β⁻¹,      h_ℓ = Q_ℓ* Q_ℓ ≥ 0,
    H_FF(β) = Σ_ℓ h_ℓ = Σ_ℓ [ E_ℓ² + (β²/16) |∇_ℓ S|² - (β/8) Δ_ℓ S ].

Facts:

1. h_ℓ ψ_β = 0 for every ℓ, so H_FF is frustration-free with ground energy 0. **[proved]** (Q_ℓ ψ_β = ½ ψ_β ∇_ℓ 1 = 0.)
2. The ground state is unique on the full tensor-product space L²(SU(2)^links). **[proved]** If H_FF φ = 0 then Σ_ℓ ∫ |∇_ℓ(φ/ψ_β)|² ψ_β² = 0, so φ/ψ_β is constant on the connected manifold SU(2)^links.
3. H_FF commutes with every local gauge transformation, so it acts on the Gauss-law subspace with no energetic Gauss term. **[proved]** S is gauge invariant; E_ℓ², |∇_ℓ S|² and Δ_ℓ S are built from bi-invariant objects on each link.
4. Explicit form. On S³, Δ cos θ = -3 cos θ and ∇(u·a)·∇(u·b) = a·b - (u·a)(u·b). Each plaquette has 4 links, so **[proved]**

        H_FF(β) = Σ_ℓ E_ℓ² + (3β/2) Σ_p (1 - ½ tr U_p)
                  + (β²/16) Σ_ℓ Σ_{p,p'∋ℓ} [ ½ tr W_{pp'} - (½ tr U_p)(½ tr U_{p'}) ] - (3β/2) N_p,

    where W_{pp'} is the loop around p ∪ p' with ℓ removed (W_{pp} = 1). The first two terms are KS with x = 3β/2.
5. Relation to KS. With β = 2x/3 = 16/(3 g⁴),

        H_KS/(g²/2) = H_FF(β) - (β²/16) Σ_ℓ |∇_ℓ S|² + const.      **[proved]**

    The two ground states agree to first order in 1/g⁴: KS first-order PT gives ψ ≈ 1 + (x/3) Σ_p ½ tr U_p (one-plaquette loop, excitation energy 4·3/4 = 3), and ψ_β ≈ 1 + (β/2) Σ_p ½ tr U_p. **[proved]** The difference W = (β²/16) Σ_ℓ |∇_ℓ S|² is bounded, gauge invariant, extensive, with local norm ≤ (β²/16)(2(d-1))² per link.
6. H_FF is unitarily equivalent to (1/4) × the generator of Langevin dynamics for the d-dimensional Wilson measure π_β ∝ e^{-βS} dHaar: ⟨fψ_β, H_FF fψ_β⟩ = ¼ ∫ |∇f|² dπ_β. So gap(H_FF) = ¼ × Poincaré constant of π_β for the Dirichlet form ∫|∇f|². **[proved]** This is the standard RK / stochastic-matrix-form correspondence (Henley; Castelnovo-Chamon-Mudry-Pujol; Dijkgraaf-Orlando-Reffert) applied to the Wilson action. **[cited for the general mechanism]**

So the ground state of H_FF in d spatial dimensions is the square root of the d-dimensional Euclidean Wilson measure: in 3+1D it is a 3D lattice gauge theory, in 2+1D a 2D one. At β = 0 it is the KS strong-coupling vacuum (all links j = 0).

### 2. Gap of H_FF at strong coupling (Bakry-Émery) **[proved]**

On unit S³, Ric = 2. For one plaquette, every block of the Hessian of ½ tr U_p in unit tangent directions (geodesics U e^{sX}, ‖X‖_op = 1) has norm ≤ 1, since |½ tr(A X B Y)| ≤ ‖X‖‖Y‖ for unitary A, B. So |Hess(½ tr U_p)(v,v)| ≤ (Σ_{ℓ∈p} |v_ℓ|)² ≤ 4 Σ_{ℓ∈p} |v_ℓ|². Each link is in 2(d-1) plaquettes, so ‖Hess S‖ ≤ 8(d-1), and Ric + β Hess S ≥ 2 - 8(d-1)β. Bakry-Émery then gives the Poincaré inequality with constant K = 2 - 8(d-1)β [cited: Bakry-Émery, as used by Shen-Zhu-Zhu 2204.12737], and by fact 6:

    gap(H_FF(β)) ≥ ½ - 2(d-1)β    for β < 1/(4(d-1)), uniformly in L, d = 2 or 3, SU(2).

The bound holds on the full space and therefore also in the gauge-invariant sector. In KS units the condition reads g⁴ > 64(d-1)/3, i.e. g² > 4.62 (d=2) or g² > 6.53 (d=3), which is strong coupling. This is the Shen-Zhu-Zhu argument transcribed to the Hamiltonian; the content is not new.

### 3. Yarotsky covers KS with Gauss law at strong coupling **[proved, hypotheses checked against Theorem 1 of math-ph/0412040]**

Yarotsky's setting: a translation-invariant lattice, a tensor product of possibly infinite-dimensional site spaces, a "classical" H_0 = Σ h_x that is a function of an orthogonal partition of unity containing the projector onto a preferred vector Ω_x, with Ω non-degenerate and h_x ≥ 1 on its complement, plus a perturbation relatively form-bounded with constants α, β small enough (Theorem 1), or α < 1 with a small bounded part (Theorem 2). Conclusions: unique ground state, gap uniform in volume, exponential clustering, infinite-volume limit.

KS on the unconstrained space L²(SU(2))^links: take a site to be the d links leaving a vertex. H_0 = (scaled) Σ E_ℓ², which is a function of the isotypic projectors {Π_j} per link. Ω = constant function (all j = 0) is non-degenerate, with gap 3/4 per link. The plaquette term is bounded (α = 0, β ∝ 1/g⁴). For g large, Theorem 1 gives a unique gapped ground state on the full space, uniformly in L. The gauge group G^V commutes with H and is connected with no nontrivial characters (SU(2) is perfect), so the unique ground state is gauge invariant and the gap on the Gauss-law subspace is at least the full-space gap.

So the strong-coupling KS gap including the Gauss law is a direct corollary of Yarotsky. The constants are not explicit in Theorem 1. This answers the STATE.md open question.

### 4. Finite-size criteria, stated for this model

Reduction to projectors. Let P_ℓ = M_ψ E_ℓ M_ψ⁻¹, where E_ℓ is the heat-bath conditional expectation of π_β (integrate U_ℓ against its conditional law given the staples) and M_ψ is multiplication by ψ_β. P_ℓ is an orthogonal projection on L²(Haar), range P_ℓ = ker h_ℓ, and P_ℓ commutes with gauge transformations. **[proved]** P_ℓ and P_ℓ' commute unless ℓ and ℓ' share a plaquette (Fubini). **[proved]** Fiberwise in the other links, h_ℓ is the single-link operator

    h(κ) = E² + (κ²/16) sin²θ' - (3κ/8) cos θ',     κ = β |a_ℓ| ∈ [0, 2(d-1)β],

where a_ℓ ∈ R⁴ is the staple sum and θ' is the angle from a_ℓ/|a_ℓ|. Its gap γ(κ) is ≥ 0.75 for every κ on the grid 0..32, with the minimum 3/4 at κ = 0 and γ ≈ κ/4 at large κ (`checks.log`). **[numerical, grid, not certified]** Hence

    H_FF ≥ γ_* H_proj,   H_proj = Σ_ℓ (1 - P_ℓ),   γ_* = inf_κ γ(κ) (= 3/4 numerically).

H_proj is the generator of random-scan heat-bath Glauber dynamics of π_β, conjugated to L²(Haar). It is FF with the same unique ground state.

**Criterion K (Knabe, cube patches) [proved; standard argument, Knabe 1988, Gosset-Mozgunov 1512.00088]:** Let H = Σ_e Π_e with Π_e orthogonal projections on any Hilbert space (finite or infinite dimension), FF. For a family of patches A with H_A = Σ_{e∈A} Π_e, suppose each e is in m_1 patches, each non-commuting pair is in exactly m_2 patches, each commuting pair in at most m_2 patches, and H_A² ≥ γ H_A for all A. Then

    H² ≥ (1 + (m_1/m_2)(γ - 1)) H.

Proof: Σ_A H_A² = m_1 H + Σ_{e≠e'} n(e,e') Π_eΠ_e' ≤ m_1 H + m_2 (H² - H), using Π_eΠ_e' + Π_e'Π_e = 2Π_eΠ_e' ≥ 0 for commuting pairs; the left side is ≥ γ m_1 H. No dimension or finiteness of the local space is used.

Applied to H_proj on a torus Λ_L (L ≥ k+2), patches = all translates of the cube B_k of side k (in plaquettes) with all its links. Every non-commuting pair lies in a common plaquette, whose bounding box is a unit square, so m_2 = k²(k+1)^{d-2} for every such pair. Every other pair has a bounding box at least as large, so it is counted at most m_2 times. m_1 = k(k+1)^{d-1}. **[proved]** Therefore

    gap(H_proj) ≥ 1 + ((k+1)/k) (γ_k - 1),    threshold  γ_k > 1/(k+1),

    gap(H_FF) ≥ γ_* · [1 + ((k+1)/k)(γ_k - 1)],

where γ_k = ess-inf over exterior link configurations b of the smallest nonzero eigenvalue of H_{B_k}(b) = Σ_{ℓ∈B_k}(1 - P_ℓ). This is the gap of heat-bath dynamics of the Wilson measure in the box with frozen exterior links. The threshold 1/(k+1) is the same in d = 2 and d = 3. At β = 0, γ_k = 1 and the bound gives gap ≥ γ_* = 3/4, which is exact.

**Criterion D (two-link angles, Dobrushin-like) [proved here]:** For projections with kernels P, Q (A = 1-P, B = 1-Q), AB + BA ≥ -c(A + B) with c = ‖PQ - P∧Q‖. Proof: by Jordan's lemma reduce to 2×2 blocks; on a block, A and B are rank-one projections onto x, y with ⟨x,y⟩ = t, |t| ≤ c, and AB + BA + c(A+B) has eigenvalues (t+c)(1+t) and (c-t)(1-t), both ≥ 0. Random-matrix check: minimum eigenvalue -4e-15 over 2000 cases (`checks.log`). Summing gives

    gap(H_proj) ≥ 1 - max_ℓ Σ_{ℓ'~ℓ} c_{ℓℓ'},

where ℓ' ~ ℓ means sharing a plaquette (6(d-1) neighbours), and c_{ℓℓ'} = sup over the other links of the maximal correlation of (U_ℓ, U_ℓ') under the conditional Wilson measure (the Friedrichs-angle cosine between range E_ℓ and range E_ℓ'). Criterion: c < 1/(6(d-1)), i.e. 1/12 in d = 3. A crude Doeblin bound c ≤ 1 - e^{-4β} only gives β < 0.022 (d=3), which is worse than Bakry-Émery. The exact c(β) is a 6-dimensional singular-value problem and has not been computed.

**Certified numerical input needed (for the computational agent):**
- (a) γ_* = inf_κ gap h(κ): a 1D radial Sturm-Liouville problem per sector L. A certified lower bound (Temple / interval FD) on κ ∈ [0, 2(d-1)β] is cheap.
- (b) c(β) for criterion D: the second singular value of the two-link conditional kernel on SU(2)×SU(2), sup over staple data (a, b ∈ R⁴ with |a|, |b| ≤ 2(d-1) - 1 plus the shared-plaquette holonomy M ∈ SU(2)). Peter-Weyl truncation with an analytic tail bound (the kernel is entire in β) plus a Lipschitz net over the compact parameter set.
- (c) γ_k for criterion K: the heat-bath gap on L²(SU(2)^{|B_k|}) (12-dim for k = 1, d = 2; 36-dim for d = 3), uniform over continuous boundary data. Not certifiable with current tools for continuous SU(2). It becomes a finite but large problem (120^{links}) with the binary icosahedral subgroup 2I as a regulator.

None of this needs `work/20261120-006` (certified KS bounds). Those bounds are for KS, which is not FF.

### 5. Is KS itself frustration-free?

- Natural decomposition (terms E_ℓ² and plaquette terms, with any constant shifts): frustrated. **[proved]** A common ground state of all E_ℓ² is the all-j=0 state, and ½ tr U_p applied to it is a nonzero one-plaquette loop, so it is not an eigenvector of any plaquette term.
- Exact quasi-local FF form. In finite volume the KS ground state ψ_KS is strictly positive and smooth: -Δ + V on a compact connected manifold has a positivity-improving semigroup [cited: standard Perron-Frobenius for Schrödinger operators, Reed-Simon IV §XIII.12]. With Φ* = -2 log ψ_KS and Q_ℓ* = ½(∇_ℓ + ½∇_ℓΦ*),

        H_KS/(g²/2) - E_0 = Σ_ℓ Q_ℓ** Q_ℓ*,  each term ≥ 0 and annihilating ψ_KS.   **[proved]**

    So KS is FF with these terms, but each term involves ∇_ℓ Φ*, which depends on all links. The terms are local iff Φ* is a finite-range sum.
- **[conjectured]** For every finite g and all large enough L, Φ* is not a finite-range sum (only quasi-local, with range growing as g → 0). Evidence: the strong-coupling expansion of log ψ_KS has connected n-plaquette terms at order (1/g⁴)^n for all n, and adjacent-plaquette second-order amplitudes differ from products of first-order ones (shared-link energy denominators). Hence KS is not FF with finite-range terms, and Knabe-type criteria do not apply to KS directly.

### 6. Obstructions

- O1. Infinite-dimensional link spaces. These do **not** block criteria K, D or the martingale method, which are dimension-free operator inequalities. **[proved for K, D]** They do block (i) certified patch numerics: the patch problem is a fiber family over continuous boundary data on a high-dimensional manifold; and (ii) the stability theorems (Bravyi-Hastings-Michalakis 1001.0344, Michalakis-Pytel 1109.1588), which assume finite local dimension. **[cited assumptions]** Lieb-Robinson bounds with unbounded on-site terms exist (Nachtergaele-Ogata-Sims math-ph/0603064; Nachtergaele-Raz-Schlein-Sims 0712.3820), but no FF stability theorem was found in that setting.
- O2. Non-commuting terms. Not an obstruction: K and D handle them.
- O3. H_FF → KS. H_KS ∝ H_FF - W with local ‖W‖ = O(β²(d-1)²). Any perturbative transfer of the H_FF gap to KS needs local ‖W‖ small compared with gap(H_FF), i.e. small β, which is the strong-coupling region already covered by Yarotsky (section 3). The FF family does not extend the KS coupling range. **[proved for perturbative transfer; conjectured that no other transfer works]**
- O4. Continuum. H_FF's spectrum is the Langevin relaxation spectrum of the d-dimensional Wilson measure. This is RK dynamics: not Lorentz invariant, with diffusive (z = 2) scaling at criticality [cited: Henley cond-mat/0311345, Isakov et al. 1012.3806]. Heuristic for 2D U(1) and Gaussian modes: gap(H_FF) ~ 1/β → 0 as β → ∞. **[conjectured]** In KS units β ∝ 1/g⁴ (not 1/g²), so H_FF's ground state matches KS only at strong coupling. Its continuum limit, if any, is not 3+1D Yang-Mills.
- O5. String-net / quantum-double FF points (Levin-Wen SU(2)_k, Kitaev D(G) for G ⊂ SU(2) finite) are commuting-projector and exactly gapped, but they sit in deconfined topological phases. Reaching the confined (KS strong-coupling) phase from them requires a gap closing (topological order cannot be created by finite-depth local evolution: Bravyi-Hastings-Verstraete quant-ph/0603121). **[cited]** Wrong phase for YM.

## Verdict

**Alive with conditions, as a stepping-stone program. Dead as a route to the KS or continuum Yang-Mills gap beyond strong coupling.**

- Alive: H_FF(β) is a gauge-invariant, frustration-free SU(2) lattice gauge Hamiltonian on the true (infinite-dimensional) link space. Its unique ground state is √(Wilson measure) and its gap is proved for β < 1/(4(d-1)) (section 2). Criteria K and D apply verbatim in infinite dimension and can push the proved range further with certified two-link or patch numerics. No non-Abelian frustration-free gauge Hamiltonian with a proved gap was found in the literature. INSPIRE shows no non-Abelian follow-up to 2409.10612, and Gorantla-Huang's FF point is a different, topological construction.
- Conditions and limits: (i) everything H_FF proves at small β, Yarotsky already proves for KS itself; (ii) the H_FF → KS map is only first-order in 1/g⁴; (iii) H_FF has RK (z = 2) dynamics, so no continuum statement about Yang-Mills follows. A result here is a Result-A candidate only on a "non-Abelian RK gauge model" stepping stone. It is partial and says nothing about the Millennium problem.
- Result-C content: O3 + O4 + section 5 together form a precise obstruction for line 2. KS is FF only with quasi-local terms, local FF deformations agree with KS only at O(1/g⁴), and FF (RK) dynamics are diffusive.

## Next theorem to attempt

**Target T (2+1D, SU(2), all couplings, candidate).** Let d = 2, Λ a finite box with free boundary, and H_FF(β) as above. Then there is c > 0, independent of Λ and β, with gap(H_FF(β)) ≥ c · min(1, 1/β) for all β ≥ 0.

Why tractable: with free boundary in d = 2 the plaquette holonomies (based at a root by a maximal tree) are independent under π_β [cited: standard for 2D lattice YM, e.g. Chatterjee 1803.01950]. So the ground state is an exact product state in plaquette variables. The task is a Poincaré inequality for a product measure with a non-product Dirichlet form: each interior link moves the two adjacent plaquette holonomies. Plan: (1) check on 1×n and 2×2 boxes numerically, and on the U(1) analogue analytically, that gap ~ 1/β is right; (2) prove a block-Poincaré lemma for two plaquettes sharing a link, then tensorize along a comb tree (criterion K with k = 1 strips, or a martingale over rows). If true, T is the first all-coupling gap for a non-Abelian frustration-free gauge Hamiltonian. The 1/β decay is expected (O4) and must be stated as such.

**Short-term certified step (feeds T and d = 3):** compute and certify γ_* and c(β) (items 4a, 4b) and report the largest β where criterion D beats Bakry-Émery's β < 1/(4(d-1)).

## Open issues

- γ_* ≥ 3/4 is from a finite-difference grid (κ ≤ 32, N = 4000); it needs a certified bound, or a proof that γ(κ) is non-decreasing.
- The "unit-norm block ≤ 1" Hessian estimate in section 2 should be re-derived by a second person; constants depend on the tangent-norm convention.
- Section 5's conjecture (no finite-range FF decomposition of KS) has evidence only. A check would compute log ψ_KS on a 3-plaquette open chain with the existing ED and test additivity.
- The novelty of H_FF as a written-down SU(2) model is not established: the RK construction for "any field theory with an explicit real action" is known (Fendley-Fradkin cond-mat/0502071 make this remark), and lattice stochastic quantization of gauge fields is old. Treat H_FF as known in principle; only the gap theorems (T, and criterion D numbers) would be new.
- The infinite-volume statement (uniform finite-volume gap ⇒ gap of the GNS Hamiltonian) was not written out.
- Aizenman-Holley 1987 and L. Wu 2006 (Dobrushin ⇒ Glauber gap/Poincaré) are cited from memory, not verified; criterion D is proved here independently.
