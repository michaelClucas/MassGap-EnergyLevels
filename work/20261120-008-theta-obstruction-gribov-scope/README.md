# 20261120-008: theta/Bloch obstruction (line 4, claim C4) and Payne-Weinberger/Gribov scoping (line 3)

Status tags: **[P]** proved here (argument given), **[C]** cited, **[N]** numerical observation, **[J]** conjectured.
Every statement carries its setting. Nothing here is a statement about the Millennium problem being solved.

Files:
- `theta_bands.py` -> `theta_bands.csv`, `theta_summary.txt` (Task A, pendulum bands; < 1 s).
- `gribov_scope.py` -> `gribov_widths.csv`, `gribov_summary.txt`, `run.log` (Task B; ~10 min on 1 thread).
- Rerun: `python tools/pin.py .venv/Scripts/python.exe work/20261120-008-theta-obstruction-gribov-scope/<script>`.

---

## Task A. What the theta/Bloch band structure controls, and why it does not give the gap

### Question
The theta-vacuum picture (Jackiw-Rebbi 1976; Callan-Dashen-Gross 1976): large gauge transformations act on configuration space as a Z translation T, [H, T] = 0, so H = ∫⊕ H_θ dθ/2π. Does an electron-style "Bloch theorem + band gap" argument bound the YM mass gap?

### Definitions (any setting)
Let H be self-adjoint, T unitary, [H,T] = 0, T generating a free Z action on the configuration space.
- Fiber H_θ: H restricted to ψ with Tψ = e^{iθ}ψ. E_0(θ) ≤ E_1(θ) ≤ ... its levels.
- **Direct (fiber) gap** Δ(θ) = E_1(θ) − E_0(θ). This is the physical mass gap of the theory at vacuum angle θ.
- **Band (indirect) gap** Δ_ind = min_θ E_1(θ) − max_θ E_0(θ). This is what the electron-level "band gap" argument controls: a gap in spec(H) above the lowest band.
- **Bandwidth** W_0 = max E_0 − min E_0; **curvature** χ = E_0''(0) (in field theory E_0(θ) = V ε(θ) and ε''(0) is the topological susceptibility, [C] e.g. Witten 1998 below).

### Steps
A1 **[P], abstract.** H ≅ ∫⊕ H_θ and spec H = closure of ∪_θ spec H_θ (standard direct-integral decomposition, [C] Reed-Simon IV §XIII.16). Bloch's theorem is a block diagonalization: it gives no lower bound on any Δ(θ). Every gap statement must still be proved fiber by fiber.

A2 **[P], abstract.** Δ_ind ≤ Δ(θ) for every θ, so Δ_ind > 0 implies Δ(θ) > 0 for all θ. The converse fails: Δ(θ) > 0 for all θ is compatible with Δ_ind < 0 (overlapping bands). So the band-gap argument can only ever deliver a direct gap through an indirect one.

A3 **[N], QM, quantum pendulum** H = p²/2 + V0(1 − cos x) on R, T: x → x+2π (I = 1, ω = √V0), `theta_summary.txt`:

| V0 | ω | W_0 | χ = E_0''(0) | Δ(0) | Δ(π) | Δ_ind |
|---|---|---|---|---|---|---|
| 0.05 | 0.224 | 1.02e-1 | 2.48e-2 | 0.502 | 0.0500 | 0.0500 |
| 1 | 1.000 | 2.67e-3 | 1.32e-3 | 0.878 | 0.822 | 0.822 |
| 4 | 2.000 | 2.71e-6 | 1.36e-6 | 1.866 | 1.865 | 1.865 |
| 8 | 2.828 | 6.14e-9 | 3.07e-9 | 2.697 | 2.697 | 2.697 |

- Bands control vacuum-energy data: W_0 ∝ s^{3/2} e^{−8s} with s = √V0 (8s is the instanton action ∫√(2V)dx; log W_0 + 8s − 1.5 log s = 2.07, 2.12, 2.14, 2.16 for V0 = 1..8) and χ = W_0/2 to 3 digits for V0 ≥ 2 (tight-binding E_0 = c − (W_0/2) cos θ). [N] (V0 = 16 row is at double-precision limit, not used.)
- The direct gap Δ(θ) ≈ ω − 1/8 is the single-well oscillator spacing (curvature of the well plus first anharmonic correction), independent of θ to the shown digits once tunneling is small. It is not produced by the Bloch structure; the bands only add an O(W) θ-dependence. [N]
- Weak potential: Bragg gap at the zone edge Δ(π) = V0 (2|V_1|) to 6e-4 relative, V0 ≤ 0.05. [N]

A4 **[P], product model (finite N wells, common θ).** Take N copies with T_i acting on copy i and impose T_iψ = e^{iθ}ψ for all i (θ couples to the total winding, exactly the structure of the dilute instanton gas). Then H_θ = ⊕-sum of single-well H_θ's: E_0^{(N)}(θ) = N e_0(θ), Δ^{(N)}(θ) = Δ(θ) (intensive), W_0^{(N)} = N W_0 (extensive). Hence Δ_ind^{(N)} < 0 for N > Δ/W_0 while every fiber keeps its gap (V0 = 4: bands overlap from N ≈ 6.9e5). [P] (proof: tensor-product spectrum), numbers [N].
In 4D YM the same holds with N → volume: ε(θ) is a nonconstant function of θ, E_0(θ) = V ε(θ), so W_0 ∝ V and Δ_ind → −∞ as V → ∞ [C for ε nonconstant: dilute-gas formula of Callan-Dashen-Gross; Witten 1998; nonzero lattice χ_top]. In infinite volume the θ's are superselection sectors and only Δ(θ) has meaning.

A5 **[P/C], two exact 2D counterexamples showing band data and the mass gap are logically independent.**
- Pure U(1) gauge theory on a spatial circle of length L (continuum, 1+1D, finite volume): only the zero mode survives Gauss's law, H = (g²L/2)Π², T = large gauge shift of the holonomy, E_n(θ) = (g²L/2)(n + θ/2π)² exactly [P]. Vacuum energy density ε(θ) = (g²/2)(θ/2π)² on |θ| < π, χ = g²/4π² > 0, level crossing at θ = π (Bragg gap exactly zero, there is no periodic potential). The fiber gap at θ = 0 is g²L/2 → ∞: an electric-flux (string) energy σL with σ = g²/2, not a particle; there are no particles at all. Nonzero band curvature, no glueball.
- Massless Schwinger model (QED_2, continuum, infinite volume): θ is removed by a chiral rotation so E_0(θ) is θ-independent (χ = 0), yet the spectrum has a gap e/√π [C: Coleman-Jackiw-Susskind, Ann. Phys. 93, 267 (1975); Coleman, Ann. Phys. 101, 239 (1976), web-verified via OSTI; Schwinger 1962]. Zero band curvature, nonzero gap.

A6 **4D SU(2)/SU(N), continuum Hamiltonian on T³ of size L.**
- Small L **[C]** (Lüscher, NPB 219, 233 (1983); van Baal review hep-ph/0008206 Sec. 3, web-read): the 0++ mass is m_0 = ε_1 g^{2/3}(L)/L, ε_1 ≈ 4.1167, set by the quartic effective potential in the constant-abelian (toron) valley. This is the "single-well spacing" of A3: a bound-state (line 1/3) result, θ-independent to all orders of perturbation theory.
- **Where the band picture does give something [C]** (Koller-van Baal NPB 302, 1 (1988); van Baal hep-ph/0008206 Sec. 3): (a) tunneling between gauge-equivalent vacua through the sphaleron (barrier 72.605/(L g²)) is implemented as Bloch boundary conditions Ψ → e^{iθν}Ψ on the fundamental domain; van Baal states explicitly "the relevant 'Bloch momentum' is in this case the θ parameter". Instanton effects (θ-dependence of low levels, onset of χ_top) switch on at L ≈ 0.75 fm. (b) Electric-flux energies ('t Hooft, NPB 153, 141 (1979)) are Bloch momenta with respect to center-twisted gauge transformations; their energies are tunneling splittings between toron vacua, ∝ exp(−S/g) with S ≈ 12.46 (the "pinchon"), and they switch on as the wave functional spreads over the vacuum valley. At large L the flux energy grows like σL (confinement). So the band structure computes θ-dependence of levels in intermediate volume and the string tension, not the 0++ gap.
- Large L / infinite volume: A4 applies (W_0 ∝ V). The dilute instanton gas is IR divergent in the instanton size [C: 't Hooft PRD 14, 3432 (1976); CDG], so even the tight-binding form ε(θ) ∝ −cos θ is unjustified. At large N, ε(θ) = min_k N² h((θ + 2πk)/N): multiple branches crossing at θ = π [C: Witten, PRL 81, 2862 (1998)]. At θ = π the mixed center/time-reversal anomaly forbids a trivially gapped unique vacuum [C: Gaiotto-Kapustin-Komargodski-Seiberg 1703.00501]: the analog of the Bragg gap at the zone edge is expected to vanish (bands 0 and 1 cross), the opposite of the electron mechanism A3.

### Result A (obstruction, precise)
Blocked technique: **Bloch's theorem along large gauge transformations + band-gap argument (gap of H above the lowest band, with the band gap estimated by tight-binding single-well spacing or nearly-free Bragg splitting).**
Blocking structural properties:
1. **[P] One collective coordinate in an infinite-dimensional space.** T translates only the Chern-Simons coordinate; each fiber H_θ is the full interacting theory at angle θ. Bloch reduction removes one degree of freedom and leaves the gap question untouched (A1). The quantity the band structure controls, E_0(θ) (ε(θ), χ_top, flux energies), is logically independent of the fiber gap (A5, two exact counterexamples).
2. **[P in product model, C for YM] Extensivity.** Bandwidth ∝ volume, fiber gap intensive, so in infinite volume the band gap Δ_ind = −∞ and the only meaningful statement is Δ(θ) > 0, i.e. the original problem (A2, A4).
3. **[C/J] Tight-binding input fails.** The "single well" is the perturbative vacuum, whose spacing is 2π/L for free gluons and ε_1 g^{2/3}/L for the toron valley, → 0 as L → ∞ at fixed coupling. The well picture yields no L-independent gap; the real gap needs the nonperturbative running coupling, i.e. the problem itself.
4. **[C] Bragg input fails at θ = π** (Witten large-N branches; GKKS anomaly).

What survives: in small/intermediate volume on T³ the Bloch boundary conditions are the correct way to include tunneling (Koller-van Baal), and the band widths give θ-dependence and electric-flux energies. This is prior art, not new. The finite-volume glueball gap there comes from the well curvature (Lüscher), a bound-state technique.

**Verdict A: Result-C obstruction for line 4, C4 can move from conjectured to "obstruction, argued" (steps 1-2 proved in the model/abstract setting, 3-4 cited).** No band-structure technique gives Δ(θ); line 4 is dead as a gap-proof route. It can still give the electric-flux (string tension) energies in finite volume, which belongs to a confinement program, not the gap.

---

## Task B. Payne-Weinberger bound on the Gribov region: scoping

### Question
Does λ_1 ≥ π²/diam² (Payne-Weinberger, ARMA 5, 286 (1960), corrected for n ≥ 3 by Bebendorf, Z. Anal. Anwend. 22, 751 (2003); log-concave weights: Ferone-Nitsch-Trombetti arXiv:1207.0680, all web-verified) for the (weighted) Laplacian on the Gribov region give a usable gap bound for YM on a finite lattice?

### Setting used for numerics
**Noncompact (linear-A) SU(2) lattice Coulomb gauge**, spatial Z_L^d periodic, d = 2, 3, a = 1, finite volume. FP form q_A(ω) = Σ|ω(x+μ) − ω(x)|² + Σ A_μ(x)·(ω(x) × ω(x+μ)) on ghosts orthogonal to constants (the linear-A limit of the compact lattice Landau FP form). Ω = {A transverse : q_A ≥ 0}. Caveat **[P]**: for non-Abelian groups the noncompact lattice action is not invariant under lattice gauge transformations, so the fundamental modular region Λ is not defined in this discretization. It is a stand-in for the continuum theory with a UV cutoff, where M = −∂·D(A) is affine in A. For the compact Wilson lattice, the configuration space is SU(2)^E; M(U) is affine in the R⁴ embedding coordinates but the domain is a product of 3-spheres, so Euclidean convexity is meaningless and geodesic convexity of Ω is unknown **[J/open]**.

### (i) Convexity, boundedness, diameter
- **[P]** Ω is convex: q_A(ω) is affine in A for each ω, so Ω is an intersection of half-spaces. Same argument in the continuum ([C] Zwanziger NPB 209, 336 (1982); Dell'Antonio-Zwanziger CMP 138, 291 (1991)).
- **[C + N]** Bounded in every direction: the A-linear part is traceless on the constants-orthogonal space, so it has a negative eigenvalue unless it vanishes. Numerically, every sampled direction had finite extent for L ≥ 3. **[N] Exception: L = 2** (d = 3) gave a width ~5e6, i.e. Ω is effectively unbounded along directions where the two links joining the same pair of sites cancel. L = 2 is a lattice artifact, excluded from fits.
- **[N] Diameter** (lower bound on diam Ω from gradient ascent of the width t_+ + t_- over unit directions in the L² norm Σ_{x,μ,a}(A_μ^a)²; `gribov_widths.csv`):

| d | L | V | dim Ω | D_lower | median width | softest-mode width |
|---|---|---|---|---|---|---|
| 2 | 3 | 9 | 30 | 31.8 | 20.9 | 22.9 |
| 2 | 4 | 16 | 51 | 32.8 | 20.9 | 19.6 |
| 2 | 5 | 25 | 78 | 35.8 | 23.4 | 18.9 |
| 2 | 6 | 36 | 111 | 41.0 | 26.8 | 18.4 |
| 2 | 8 | 64 | 195 | 51.4 | 31.4 | 17.7 |
| 3 | 3 | 27 | 165 | 72.4 | 43.6 | 39.8 |
| 3 | 4 | 64 | 387 | 91.9 | 61.9 | 39.2 |

  D_lower grows roughly linearly in L in both d = 2 (L = 5..8: slope ~5) and d = 3 (L = 3 → 4: 72 → 92), i.e. D² ∝ V in d = 2 and between L² and L³ in d = 3 (two points, D²/V = 194, 132). The width along the softest (lowest-|p|) mode is L-independent in d = 2 (23 → 18) and in d = 3 (39.8, 39.2), and the median random-direction width also grows. The scaling is an observation on tiny lattices **[N]**. Gradient ascent gives only lower bounds on D, which is the direction needed: the PW bound is at most π²/D_lower². Gradient ascent gives only lower bounds on D, which is the direction needed: the PW bound is at most π²/D_lower².
- **Continuum T³ analytic bound [C + P]:** along constant abelian fields A_j = C_j σ_3/(2L) the fundamental domain is |C_j| ≤ π and the horizon is at |C_j| = 2π ([C] van Baal NPB 369, 259 (1992) and hep-ph/0008206 Sec. 3). With ||A||² = Σ_a ∫(A_j^a)², the diagonal C = (c,c,c) segment in Λ has length 2π√(3L) [P], so D(Λ) ≥ 2π√(3L).

### (ii) Log-concavity of the FP determinant
- **[P]** On Ω° (M > 0), f(t) = log det M(A + tB) has f'' = −tr[(M^{−1/2} M_1(B) M^{−1/2})²] ≤ 0, since M is affine in A. So J = det M is log-concave on Ω.
- **[N]** Max second difference of log det along 10 random chords in Ω (d = 3, L = 3): −0.070 (all ≤ 0, concave as proved).

### (iii) Coulomb-gauge Hamiltonian
[C] Christ-Lee, PRD 22, 939 (1980) (web-verified abstract): in Coulomb gauge
H = (g²/2) ∫ J^{−1} Π_tr · J Π_tr + (g²/2) ∫ J^{−1} ρ^a [M^{−1}(−∂²)M^{−1}]_{ab} J ρ^b + (1/2g²) ∫ B² + V_1 + V_2 (order-ħ² Christ-Lee terms), with ρ^a = f^{abc} A^b·Π^c the gluon colour charge and J = det M.
- The first term is −(g²/2) J^{−1}∇·(J∇): a weighted Laplacian with log-concave weight J (by ii) on convex Ω. **[P given Christ-Lee]**
- The Coulomb term is a second-order operator in Π with A-dependent coefficients, nonnegative on Ω° and singular at the horizon (M^{−1}). Together with the first term it is a Laplace-Beltrami-type operator for a non-flat metric G(A) on Ω (the orbit-space metric, [C] Babelon-Viallet CMP 81, 515 (1981)); Euclidean convexity of Ω is not geodesic convexity for G. **[J]**
- Domain: the physical space is the fundamental modular region Λ ⊂ Ω (convex [C] van Baal hep-th/9711070; Dell'Antonio-Zwanziger) with boundary points identified by gauge transformations (including large ones, carrying the θ phase). **[P] Lemma:** for a flat weighted Laplacian with weight μ on Λ, the θ = 0 quotient problem has form domain contained in H¹(Λ, μ) with the same quadratic form, so by min-max its eigenvalues are ≥ the Neumann eigenvalues on Λ, both with λ_0 = 0 (constants). Hence λ_1(quotient) ≥ λ_1^N(Λ, μ) ≥ π²/diam(Λ)² for log-concave μ on convex Λ (FNT 2012). For θ ≠ 0 constants are not admissible and the gap comparison fails.
- So the **only provable statement** is: in the continuum-with-cutoff setting, at θ = 0, the flat part of the kinetic operator alone has gap ≥ (g²/2)π²/diam(Λ)². **[P modulo cited convexity of Λ]**

### (iv) Failure point and scaling
1. **[P] The gap is not monotone in added nonnegative terms.** H ≥ K does not imply gap(H) ≥ gap(K) (a double well added to a free particle closes the gap). The magnetic term and the Coulomb term cannot simply be dropped. The correct route is the ground-state transform: H − E_0 ≅ weighted Laplacian with weight J ψ_0², which needs ψ_0² log-concave on Λ (Bakry-Emery CD(0,∞)). That is exactly the unproved hypothesis of Moncrief-Marini-Maitra (1809.06318) and Mondal (2301.06996).
2. **[N] Brascamp-Lieb cannot supply it**: the noncompact lattice magnetic energy (d = 3, L = 2, transverse slice) has Hessian min eigenvalue 0 at A = 0 and 0.040 at a point with |A| = 2.8 inside Ω, but −7.0, −3.3, −11.3 at points with |A| = 13.2, 9.6, 16.4 inside Ω (random rays, each at a fraction 0.5-0.95 of the distance to the horizon). The potential is not convex on Ω, so log-concavity of ψ_0 does not follow from Brascamp-Lieb. (L = 2 is the artifact lattice of (i), but the points used are inside Ω on finite rays; the non-convexity is generic: V = αt² + βt³ + γt⁴ along a line is non-convex when 3β² > 8αγ.)
3. **[N/P] Scaling.** Even the kinetic-only bound vanishes. Taking D² ≈ c V with c = 163 (mean D²/V of d = 3, L = 3, 4), a·Δ_PW ≤ (g²/2)π²/(c N_s³). Against the SU(2) 4D 0++ mass (m/√σ = 3.781, AT2021, read from `data/`) at L√σ = 3: Δ_PW/m_0++ ≤ 1.7e-4 (g² = 1, N_s = 4) down to 2.6e-6 (g² = 1, N_s = 32); ≤ 6.7e-4 at g² = 4, N_s = 4 (full table in `gribov_summary.txt`). At fixed physical L this falls as g²(a)·a² (D² ∝ V) or at best g²(a)·a (if D ∝ L, the slower growth allowed by the two d = 3 points); either way → 0 in the continuum limit, and → 0 as L → ∞ at fixed a. In the continuum small-volume T³ setting [P from (i)]: Δ_PW ≤ (g²/2)π²/(12π²L) = g²/(24L), versus Lüscher's m_0 = 4.1167 g^{2/3}/L: ratio ≤ g^{4/3}/98.8 → 0 as L → 0 (asymptotic freedom). The PW route is parametrically too weak even where the gap is known.

**Verdict B: Result-C obstruction; not worth a full proof attempt.** Precise failure point: the gap of the reduced Hamiltonian equals the PW-type quantity only after a ground-state transform, which requires log-concavity of J ψ_0² on Λ (equivalently Bakry-Emery CD(0,∞) for the orbit-space measure), and that is the open Moncrief/Mondal hypothesis, not supplied by convexity of the magnetic potential (which fails, [N]). Even granting it, the diameter bound in physical units scales as g² a²/L³ if D² ∝ volume (g² a/L² if D ∝ L; tiny-lattice data cannot separate these), vanishing in both limits and ≤ 7e-4 of m_0++ on every lattice tried, and is parametrically below the known small-volume gap by g^{4/3}. The salvageable piece is a trivial finite-dimensional theorem about the kinetic operator alone (θ = 0, cutoff continuum, flat metric), which says nothing about YM.

---

## Open issues
- Boundedness of Ω on the noncompact lattice is cited from the continuum argument plus sampling; a direct proof that the traceless A-part is nonzero on constants-orthogonal ghosts for every nonzero transverse A, L ≥ 3, was not written. L = 2 is a genuine exception.
- D_lower comes from 3-6 gradient-ascent restarts; true diameters may be larger (which only strengthens the verdict). d = 3 has only L = 3, 4, so the D² ∝ V fit in d = 3 is two points.
- Compact Wilson lattice: geodesic convexity of the lattice Gribov region in SU(2)^E and the analog of the Λ lemma were not examined. Zwanziger NPB 412, 657 (1994) (web-verified bibliographic data) treats the lattice FMR; not read.
- Christ-Lee operator form is quoted from the abstract-level check plus standard secondary use (Andrasiunas-Taylor 1204.1722); the exact ordering terms V_1, V_2 were not reread.
- Lüscher's ε_1 = 4.1167, the sphaleron barrier 72.605/(Lg²), the pinchon action 12.46, the 0.75 fm onset and |C_j| ≤ π vs 2π are taken from the van Baal review (ar5iv text read via web), not from the original papers.
- A5 Schwinger statements are textbook; only Coleman 1976 was web-checked.
