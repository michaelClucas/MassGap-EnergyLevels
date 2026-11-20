# 20261120-011 The 2+1D H_FF gap theorem (target T of work 007)

## Question

Target T (from `work/20261120-007-*`): d = 2 spatial lattice, free boundary, SU(2), H_FF(β) = Σ_ℓ Q_ℓ*Q_ℓ with ground state ∝ exp(-β S_W/2). Is gap(H_FF(β)) ≥ c·min(1, 1/β) for all β > 0, uniformly in volume? Is it known, true, false, or trivial?

## Verdict (read this first)

**Known qualitatively; the stated rate is open and the proposed proof route is blocked.**

1. **Known: volume-uniform gap at every β, on the torus.** Piazza-Shen, "Log-Sobolev inequality for 2D SU(N) lattice gauge theory at all couplings", arXiv:2610.35286 (math.PR, 26 Oct 2026) [arXiv-API; Theorem 1.1, sections 5-7 read from the PDF]. Theorem 1.1: for every N ≥ 2 and β ∈ R there is C_{N,β}, independent of L, with Ent(F²) ≤ 2C Σ_e ∫|∇_e F|² dμ for **every smooth F on G^{E_L}** (not only gauge-invariant F) on the torus (Z/LZ)². This is exactly the link-gradient (Langevin) Dirichlet form of 007 fact 6, up to constants: their inner product is Tr(X†Y) (twice our unit-S³ metric) and their weight is exp(Nβ_PS Σ Re Tr U_p), so β_ours = 4β_PS for SU(2). LSI implies Poincaré with the same constant, so via 007 fact 6, **gap(H_FF(β)) > 0 uniformly in L for every β, in 2D on the torus, on the full Hilbert space.** [cited] Their constant degrades fast: K = e^{4tn} with t = N²|β_PS| and block area n chosen so that M r^{n-1} ≤ 1/3, r = 1 - ½e^{-2t}, so n ~ e^{2t} and the Poincaré constant they give is doubly exponential in β (my reading of sections 5-6, not stated by them as a formula). Their section 7 shows the LSI constant must diverge as β → ∞ (≥ cβ^{1/3}) for SU(N) with 4 | N; that argument concerns LSI, not Poincaré, and does not cover SU(2).
   - Free boundary is not stated in the paper. **[conjectured]** The same two-scale proof works (the block marginal is again a convolution w^{*n}(W_∂B), their eq. 3.1, and free boundary removes the global torus constraint).
   - Consequence for this program: "first all-coupling gap for a non-Abelian FF gauge Hamiltonian" (007 Next theorem) is **not new**. It is an immediate corollary of 2610.35286 plus the ground-state transform. Only explicit β-dependence could be new.
2. **The rate min(1, 1/β) is neither proved nor contradicted.** Our rigorous variational upper bounds (below) show the volume-uniform SU(2) gap is ≤ 16.5 for all β tested (up to 1000), while each fixed box has gap ~ ½βλ_min(M_L) → ∞. Heuristic (section 4): the volume-uniform SU(2) gap is Θ(1) as β → ∞. If so, T is true but not sharp; proving it needs a Poincaré constant polynomial in β, which 2610.35286 does not give.
3. **The U(1) analogue of T is false.** Proved variational bound: inf_L gap ≤ (1/4)·4 Var(g')/Var(g) ≈ β/I_0(β)² ≈ 2πβ² e^{-2β}. The slow mode is the lifted (non-compact) plaquette flux, which interior link moves conserve; it relaxes only by phase slips over the barrier F = π. So any proof of T must use non-Abelian structure (consistent with Piazza-Shen Remark 1.2: their method needs positive Ricci curvature and excludes U(1)).
4. **The 007 proof plan fails.** Comparing the link Dirichlet form to the product (plaquette) Dirichlet form costs the factor λ_min(M_L) = 4 - 4cos(π/(L+1)) ≈ 2π²/L², so "product measure + comparison + tensorization" cannot give a volume-uniform bound. [proved, section 2]
5. **The large-β scaling guessed in 007 (gap ~ 1/β, O4) is wrong in all three settings checked:** single plaquette gap ≈ 2β (grows); U(1) volume-uniform gap ≤ ~e^{-2β}; SU(2) volume-uniform gap ≤ 16.5 and heuristically Θ(1).

Nothing here bears on the 3+1D or continuum Yang-Mills gap.

## Setting and conventions

Hamiltonian lattice, d = 2 spatial dimensions, L×L box of plaquettes with free boundary unless stated, fixed lattice spacing. SU(2): unit S³ metric per link, S = Σ_p (1 - ½ tr U_p), π_β ∝ e^{-βS}. U(1): unit circle per link, S = Σ_p (1 - cos F_p). 007 fact 6: gap(H_FF) = ¼ × (Poincaré constant λ of π_β for E(f) = Σ_ℓ ∫|∇_ℓ f|² dπ_β). All statements below are for the gauge-invariant sector; upper bounds from gauge-invariant test functions also bound the full-space gap.

## Method and results

`hff2d.py` (sections A-D), logs `run_AC.log`, `run_D.log`, `run_B.log`. Run from repo root: `python tools/pin.py .venv/Scripts/python.exe -u work/20261120-011-hff-2d-theorem/hff2d.py ACD` (about 2 min, 1 core; section B is slower). No external data used.

### 1. Single plaquette (L = 1) [reduction proved, values numerical]

Gauge-invariant functions are class functions of U_p, each of the 4 links moves U_p isometrically, so E(f) = 4∫|∇f|² and gap(H_FF) = λ_1(β), the first nonzero eigenvalue of -(w⁻¹)(w f')' on (0, π), w = sin²θ e^{β cos θ}. Ritz (Legendre, deg 40) and finite differences (n = 4000) agree to 1e-6:

| β | 0 | 1 | 2 | 5 | 10 | 50 | 200 |
|---|---|---|---|---|---|---|---|
| gap | 3.0000 | 3.1997 | 3.7945 | 7.7417 | 17.8008 | 97.973 | 397.994 |

gap → 2β - O(1). This matches the Ornstein-Uhlenbeck limit (radial degree-2 mode, eigenvalue 2β·4/4). Not certified (Ritz gives an upper bound; FD agrees).

### 2. Link vs plaquette Dirichlet form [proved]

- U(1), exact: with free boundary in 2D, d: links → plaquettes is onto, gauge classes ↔ F ∈ U(1)^{L²}, and π_β is a product of von Mises laws. For gauge-invariant f(F): E(f) = Σ_ℓ (Σ_p d_ℓp ∂_p f)², i.e. the plaquette Dirichlet form with matrix M = d dᵀ = 4I - A_grid (Dirichlet dual Laplacian), λ_min(M_L) = 4 - 4cos(π/(L+1)).
- SU(2): at a flat configuration (all plaquettes = 1, gauge-equivalent to all links = 1) every link move inserts the same Lie-algebra element into its adjacent plaquette holonomies and tree-link conjugations act trivially, so the metric on plaquette gradients is M_L ⊗ I_3. Near-flat configurations have positive measure, so the pointwise comparison constant is ≤ λ_min(M_L) + o(1).
- 1×2 strip, explicit form: with root-based holonomies W_1, W_2, right/left derivatives r_p, l_p, and conjugation invariance l_1 + l_2 = r_1 + r_2, E = 3|r_1|² + 3|r_2|² + |l_1 - r_2|². Checked against direct differentiation of the 7 link variables (section B; result in `run_B.log`).

### 3. U(1), L×L box: the analogue of T is false [proved bound, quadrature numerical]

Test function f = Σ_p v_p g(F_p), v the bottom eigenvector of M_L, g periodic with g' = 1 - λ/ρ (ρ the von Mises density, λ = 2π/∫ρ⁻¹; the Kramers choice). Exactly, gap ≤ ¼[(E g')² λ_min(M_L) + 4 Var g']/Var g. Quadrature: 2·10⁵-point trapezoid on a smooth periodic integrand.

| β | 2 | 4 | 6 | 8 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| bound, L → ∞ | 0.32 | 2.7e-2 | 1.2e-3 | 4.1e-5 | 1.2e-6 | 1.3e-10 | 1.0e-14 |
| 1/β | 0.5 | 0.25 | 0.17 | 0.13 | 0.1 | 0.067 | 0.05 |

At L = 64 the λ_min term dominates (bound 1.1e-2 at β = 10). So inf_L gap_U(1)(β) ≤ ~2πβ²e^{-2β}, and no bound c·min(1, 1/β) holds.

### 4. SU(2), large box: volume-uniform upper bound [proved bound, quadrature numerical]

Test function f = Σ_k a_k cos θ_k, θ_k the angle of the Wilson loop around the concentric square R_k of odd side s_k = 1, 3, ..., s_max (boundaries disjoint). Exact ingredients:
- Only links on ∂R_k move W_{∂R_k}, each isometrically, so E = Σ_k a_k² · 4s_k · E sin²θ_k. [proved]
- The law of W_{∂R} is the |R|-fold convolution of the plaquette law [cited: 2610.35286 eq. 3.1; standard 2D lasso independence, cf. Chatterjee 1803.01950]. With a = E cos θ, c_1 = E χ_1/3 for one plaquette: E cos θ_k = a^{n_k}, E cos²θ_k = (1 + 3c_1^{n_k})/4, Cov(cos θ_k, cos θ_m) = a^{n_m-n_k}(1+3c_1^{n_k})/4 - a^{n_k+n_m} for R_k ⊂ R_m (class-function convolution identity). [proved]
- gap ≤ ¼ min_a (aᵀDa)/(aᵀCa), valid for every free-boundary box of side ≥ s_max.

| β | 2 | 5 | 10 | 20 | 50 | 100 | 200 | 500 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| s_max = 1 | 3.99 | 8.47 | 18.2 | 38.1 | 98.0 | 198 | 398 | 998 | 1998 |
| best single loop | 3.99 | 8.47 | 10.5 | 16.0 | 24.2 | 33.8 | 47.9 | 75.9 | 107 |
| nested, s_max = 201 | 3.99 | 8.37 | 10.3 | 12.9 | 14.9 | 15.7 | 16.1 | 16.4 | 16.5 |

The best single loop scales as ≈ 3.4√β (optimal side ~ √β); the nested family saturates near 16.5. Heuristic **[conjectured]**: in the Gaussian regime (region area ≪ β) nested loops give ≈ 15β/s², the linearized (Abelian) conservation law protects modes only up to the scale s ~ √β where non-Abelian holonomies decorrelate, so the volume-uniform SU(2) gap is Θ(1) as β → ∞. This is the mechanism that U(1) lacks (section 3).

### 5. 1×2 strip (section B)

Reduced Dirichlet form vs direct 7-link differentiation: max relative deviation 5.6e-10 (and conjugation-invariance residual), so E = 3|r_1|² + 3|r_2|² + |l_1 - r_2|² is confirmed. Ritz values for the gauge-invariant gap (invariant polynomials of degree D in 1 - x_0, 1 - y_0, x⃗·y⃗; Gauss quadrature 30×30×16 in (θ_1, θ_2, φ)):

| β | 0 | 1 | 5 | 20 | 50 |
|---|---|---|---|---|---|
| D = 4 | 3.0000 | 3.1498 | 6.877 | 29.58 | 74.61 |
| D = 6 | 3.0000 | 3.1498 | 6.873 | 29.66 | 74.88 |
| 1.5β | 0 | 1.5 | 7.5 | 30 | 75 |

Gaussian-limit prediction gap ≈ ½βλ_min(M) = 1.5β (λ_min = 3) is matched at β = 20, 50. **[numerical, not converged]**: D = 6 exceeds D = 4 at β ≥ 20, which a Ritz sequence on nested spaces cannot do, so quadrature or the conditioning cutoff is limiting accuracy at the 0.3% level there. The strip has no slow modes (most links are on the boundary) and does not test volume uniformity.

## Literature check (novelty)

- arXiv:2610.35286 Piazza-Shen [arXiv-API, PDF read]: settles qualitative T on the torus (above). Already listed in `literature/02-rigorous-lattice.md`; 007 missed it.
- arXiv:2204.12737 Shen-Zhu-Zhu, 2401.13299 Shen-Zhu-Zhu [arXiv-API]: Bakry-Émery LSI/Poincaré at strong coupling, any d.
- arXiv:2302.12160 Chevyrev-Shen [arXiv-API]: invariance of the continuum 2D YM measure under the Langevin dynamic; no gap.
- Kawasaki / conservative dynamics gap ~ 1/L² (Lu-Yau, CMP 156 (1993), cited from search; Chafaï math/0211108 [arXiv-API]): the mechanism behind sections 2-3.
- Searches (web; INSPIRE title queries "Langevin"+"lattice gauge"+"spectral gap", "stochastic quantization"+lattice+2D+gap, "Yang-Mills"+"log-Sobolev"/"Poincaré"): no explicit-in-β Poincaré constant for 2D lattice YM found. Absence of evidence only.

## Issues found in work 007

- "No non-Abelian FF gauge Hamiltonian with a proved gap" and "T would be the first all-coupling gap": wrong. 2610.35286 (in our own `literature/02`) gives it in 2D on the torus at every β.
- O4 heuristic "gap(H_FF) ~ 1/β" for 2D U(1) and Gaussian modes: wrong. U(1) volume-uniform gap ≤ ~e^{-2β}; Gaussian modes at fixed L have gap ~ βλ_min(M_L), growing in β.
- Target T plan "(2) block-Poincaré + tensorize along a comb tree": blocked by the λ_min(M_L) comparison loss (section 2). A working route must block at scale ≳ e^{β} (Piazza-Shen) or ~ √β (heuristic here).

## Open issues

- Polynomial-in-β Poincaré constant for 2D SU(2) (would prove T, and the Θ(1) conjecture would sharpen it). Most promising: redo Piazza-Shen blocking with blocks of area ~ β, using the explicit convolution w^{*n} contraction rate at large β (heat-kernel regime) instead of the crude r = 1 - ½e^{-2t}.
- Free-boundary version of 2610.35286 not written out.
- U(1): a matching lower bound gap ≳ e^{-2β}·poly uniformly in L (expected, not attempted).
- Numbers are floating-point quadrature/Ritz, not interval-certified. Section 4 and 3 bounds are exact formulas; only the 1D integrals are numerical.
- Section 4's Θ(1) conjecture rests on one test family plus a heuristic; better test functions (non-class functions of several loops) could lower the bound.
