# 20261120-015 Line 6: calibrating the gap-proof toolkit on the O(3) rotor chain and spin-1 chains

## Question

On stepping-stone models whose gap is known, which electron-level gap techniques reproduce the gap, over which coupling range, and exactly where does each one fail on the way to the continuum limit? Do the failure reasons match those found for Yang-Mills (C4-C9)?

Models:
1. O(3) quantum rotor chain, 1+1D Hamiltonian lattice (continuous time, spatial spacing a = 1): H = (g/2) Σ L_i² - (1/g) Σ n_i·n_{i+1}, n_i ∈ S². Its Euclidean action is (1/2g) ∫ [ṅ² + (n_{i+1}-n_i)²], so β = 1/g in the standard β/2 (∂n)² normalization; the continuum limit is g → 0.
2. Spin-1 chain h(θ) = cos θ S·S + sin θ (S·S)²: Heisenberg θ = 0, AKLT tan θ = 1/3 (bond term ∝ P_2 + const, P_2 = ½S·S + ⅙(S·S)² + ⅓ the spin-2 projector).

## Known answers (what we calibrate against)

- O(3) continuum mass gap. Two-loop running: m a = C · (1/g) · exp(-2π/g) (1 + O(g)) **[derived here from the standard O(N) β-function coefficients b0 = (N-2)/2π, b1 = (N-2)/(2π)², for which b1/b0² = 1 at N = 3; standard]**. The exact ratio m/Λ_MS̄ = 8/e (Hasenfratz-Maggiore-Niedermayer, PLB 245 (1990) 522; listed in `literature/07`) **[cited, not re-verified this task: the 1990 paper is paywalled and a public-source check failed]**. It rests on the Zamolodchikov S-matrix plus the TBA, i.e. on integrability assumptions; there is no rigorous proof of a gap of the 2D O(3) model at weak coupling (`literature/07`). The lattice-Hamiltonian constant C (Λ_H/Λ_MS̄) was not found; nothing below depends on it.
- Haldane gap of the S = 1 Heisenberg chain: 0.41050(2) (White-Huse, PRB 48, 3844 (1993)), 0.41048(6) (Todo-Kato, PRL 87, 047203 (2001)) **[cited, abstracts checked]**. Numerical, not proved.
- AKLT: gap > 0 in the thermodynamic limit, proved (AKLT, CMP 115 (1988) 477) **[cited]**; Knabe's criterion (J. Stat. Phys. 52 (1988) 627) and Gosset-Mozgunov (arXiv:1512.00088, J. Math. Phys. 57, 091901 (2016)) **[criterion forms checked from GM's statement]**: ring gap ε°_m ≥ ((n-1)/(n-2))(ε_n - 1/(n-1)) (Knabe, n > 2, m > n) and ≥ (5/6)((n²+n)/(n²-4))(ε_n - 6/(n(n+1))) (GM, m > 2n), with ε_n the gap of the OPEN n-site chain Σ_{i<n} h_{i,i+1}, h projectors.

## Results table

Every entry: setting = finite chain, lattice spacing 1, unless "infinite volume" is stated.

| technique | model | proves a gap where | fails where | failure reason (structural) | status |
|---|---|---|---|---|---|
| Certified ED (min-max + Schur/Feshbach, angular-momentum truncation) | rotor chain | Finite chains only. Smallest g certified: L = 2 open g = 0.25 (gap [0.1318, 0.1318]); L = 3 g = 0.35 (open [0.1311, 0.1320], ring [0.1245, 0.1282]); L = 4 open g = 0.5 (l_max 3); L = 4 ring g = 0.7 ([0.232, 0.240], l_max 3); L = 5, 6 ring g = 1 (L = 6: [0.413, 0.457]); L = 8, 10 g = 2 only | L = 4 g ≤ 0.35-0.5, L = 5, 6 g ≤ 0.7, L = 8, 10 g ≤ 1 (l_max 3, 2, 1) | (i) Truncation floor is extensive: Q-space floor d needs α(g/2)(l_max+1)(l_max+2) to beat E(H_α) - E0 ~ N·O(1), so l_max ~ √(N/g), and dim grows as (l_max+1)^{2N}. (ii) Even when certified, the value is the finite-volume rotor tower ≈ g/L·O(1), not m(g): at g = 0.35, certified gap ≈ 0.13 while C·e^{-2π/g}/g ≈ 4.6e-8·C. A finite-L certificate is never a statement about m(g). | numerical-certified (floating point, Lanczos completeness assumed) |
| same | spin-1 Heisenberg | exact finite-chain gaps (finite local dimension, no truncation needed) | infinite volume | no finite-size criterion applies (not frustration-free, below) | numerical |
| Knabe / Gosset-Mozgunov | AKLT | Infinite chain: GM at n = 12 gives gap ≥ 0.3026, Knabe at n = 12 gives ≥ 0.3008 (open ε_12 = 0.36436); first n with ε_n > 1/(n-1): n = 4. Ring ED → 0.3501 (n = 12). | — (works) | — | proved (cited criterion + exact ED of a 3^12-dim matrix; floating point) |
| Knabe / GM | spin-1 Heisenberg | nowhere | at the hypothesis | Not frustration-free: E0/bond = -1.4058 (N = 12 ring) vs bond minimum -2; frustration 0.594/bond. Shifting h ≥ 0 does not help: H²≥ΔH-type local inequality needs h_b ψ0 = 0. Open chains also have edge-state splittings ε_n = 0.31, 0.20, 0.14, 0.097 (n = 6..12) → 0, which kills any open-chain criterion. | proved (structural), numerical values |
| FF deformation (path θ: Heisenberg → AKLT) | spin-1 | AKLT endpoint only | all θ < arctan(1/3) | Frustration per bond falls linearly to 0 only at AKLT (0.594, 0.488, 0.376, 0.257, 0.132, 0 along 6 equally spaced θ). Ring gap (N = 12) stays 0.48-0.66 along the path (no closing seen), but stability theorems (Michalakis-Pytel) need a perturbation small compared with an unspecified constant; the distance is O(1) per bond. | numerical (gap), structural (no transfer) |
| Knabe / GM | rotor chain | nowhere | at the hypothesis | Not frustration-free (kinetic ground state ⊗|l=0⟩ is not a ground state of -n·n), local Hilbert space L²(S²) infinite-dimensional. Same as KS (literature/03 structural facts). | proved (structural) |
| RK parent H_FF + Bakry-Émery | rotor chain (H_FF) | H_FF: gap ≥ 1 - 4β for β < 1/4, uniformly in N (ring or open) | β ≥ 1/4, i.e. g ≤ 2 when matched to H (β = 1/g²) | Hess(βS) ≥ -4β is sharp (antiparallel neighbours, equal tangent displacements); Ric(S²) = 1 is fixed. Same as C6 (Ric 2 vs Hess -8(d-1)β): curvature wins only at strong coupling. And H_FF = (2/g)H + (β²/4) Σ_i |∇_i S|² exactly at β = 1/g²: the correction is ≥ 0 and extensive, so no gap transfer to H (as C6). | proved (Hessian bound and identity derived here; random search reproduces -2 per bond to 4e-3) |
| same | H_FF at large β | true gap of H_FF > 0 for every β is expected (1D classical chain); upper bound gap ≤ (2/3)/χ = 2(1-u)/(1+u) ≈ 1/β (u = coth β - 1/β) | — | H_FF's ground state |ψ0|² = e^{βΣn·n} is the 1D classical Heisenberg chain (ξ ≈ β), dynamics are Langevin (z = 2). Its g → 0 limit is not the O(3) sigma model. Same as C8 ("H_FF dynamics has no YM continuum limit"). | upper bound proved (exact 1D correlations); all-β positivity conjectured |
| Strong-coupling series in x = 1/g² | rotor chain, infinite chain | gap/g = 1 - (2/3)x + (2/27)x² + (7/405)x³ + 0.0026194 x⁴ + O(x⁵); converges for small x (finite chain-independent coefficients; convergence radius not computed) | the continuum needs x → ∞ | Target m a ∝ (1/g) e^{-2π/g} = √x e^{-2π√x} has an essential singularity at x = ∞; no finite-radius series, truncation or Padé tail reaches it rigorously. Same structure as C2 (KS plaquette) and Osterwalder-Seiler/Yarotsky strong-coupling gaps. | coefficients numerical-exact (float RSPT, checked across L and l_max); obstruction proved (standard) |

Notes on the table:
- The task text asked for an RK ground state exp(-βS/2) with S the classical **2D** O(3) action. For a 1D spatial chain the RK state built from local terms has S = the 1D classical chain energy (above). The ground state of the actual rotor chain has |ψ0|² = equal-time marginal of the 2D Euclidean O(3) model, which is not a local Gibbs measure on the chain; so no local FF parent with that ground state was found. A 2D spatial lattice of rotors with S = 2D O(3) action gives the 2D Langevin model (Bakry-Émery: gap ≥ 1 - 8β by the same Hessian count with 4 bonds per site **[derived, not checked numerically]**), which again is a z = 2 model, not the 2+1D O(3) QFT.
- AKLT GM bound uses m > 2n = 24; Knabe uses m > n.

## Method

### Certified ED (`rotor.py`, `run_rotor.py`, `summary.py`)

Basis |l, m⟩ per site, l ≤ l_max, total M = 0 sector (contains a member of every J multiplet, so λ_0, λ_1 there are those of H). n_z, n_± matrix elements from the Condon-Shortley recursion; validated by x²+y²+z² = 1 on l ≤ l_max - 1 to 2e-16.

Upper bounds: min-max with P = {all l_i ≤ l_max}. Lower bound: the Schur-complement lemma of `work/20261120-006` (proved there): λ_k(H) ≥ λ_k(A - X/(d - U_k)) if U_k = λ_k(A) < d, with X ≥ PHQHP and D = QHQ ≥ d. New inputs here (proved):
- **X.** V = -(1/g)Σ_b h_b, h_b = n_i·n_j changes l only on the two sites of b by ±1. Split Q by the exact set R of sites above l_max; from P, h_b only reaches R ⊂ b, each R is reached by at most 2 bonds, so ‖QVPψ‖² ≤ (2/g²) Σ_b ⟨ψ, P h_b (1 - P_b) h_b P ψ⟩, and K_b = P h_b (1 - P_b) h_b P is a two-site operator computed exactly with one extra shell. X = (2/g²) Σ_b K_b.
- **d.** On Q some site has l ≥ l_max + 1, so Σ L² ≥ c = (l_max+1)(l_max+2). For α ∈ {0.3, 0.6, 0.85}: H = α(g/2)ΣL² + H_α, H_α the chain with kinetic (1-α)g/2, so D ≥ α g c/2 + LB(E0(H_α)). LB(E0(H_α)) = max of (a) Anderson bond bound Σ_b e0(two-rotor bond Hamiltonian), each certified with the same lemma at l = 10 (and ≥ -1/g trivially), and (b) one level of the same certified ED for H_α.
- Rounding: eigsh residual + 10 n ε ‖A‖₁. Lanczos completeness assumed (not proved), so: numerical-certified, finite-volume, fixed spacing.

Two runs: `run_rotor.py` (g ∈ {5, 3, 2, 1.5, 1, 0.7, 0.5, 0.35, 0.25}, L = 2-4 and part of the list; stopped for time during L = 6 open (l_max 2); rows present are in `rotor_certified.csv`) and `run_rotor.py large` (rings L = 6 (l_max 2), 8, 10 (l_max 1) and open L = 10, g ∈ {2, 1, 0.7, 0.5}; `rotor_certified_large.csv`).

### Spin-1 (`spin1.py`)

Exact ED in S_z = 0. AKLT open-chain ε_n, n = 3..12 (first eigenvalue > 1e-8 above the 4-fold edge-state ground space), Knabe and GM bounds, ring gaps. Heisenberg → AKLT path: frustration per bond (E0 - N e_min)/N and ring gap, N = 8, 10, 12. Output `spin1_aklt_knabe.csv`, `spin1_path.csv`, `run_spin1.log`.

### Bakry-Émery (`bakry_emery.py`)

S = -Σ_b n_i·n_j on (S²)^N. On S², Hess(a·n) = -(a·n) g, so for one bond Hess(-n_1·n_2)(u,u) = (n_1·n_2)(|u_1|²+|u_2|²) - 2u_1·u_2 ≥ -2(|u_1|²+|u_2|²), equality at n_2 = -n_1, u_2 = u_1 **[proved]**; per site (2 bonds) Hess(βS) ≥ -4β, Ric = 1, so gap(H_FF) ≥ 1 - 4β (Bakry-Émery) **[proved]**. Ground-state transform: Q_i*Q_i = L_i² + (β²/4)|∇_i S|² - (β/2)Δ_i S and Δ_i(a·n_i) = -2a·n_i give H_FF = Σ L_i² - 2β Σ n·n + (β²/4)Σ|∇_i S|² **[proved]**. Output `bakry_emery.csv`.

### Strong-coupling series (`sc_series.py`)

Nondegenerate Rayleigh-Schrödinger PT in the momentum-0 sectors (M = 0 for the ground state, M = 1 for the magnon) of rings. The x^n gap coefficient equals the infinite-chain one for n ≤ L - 1 and n ≤ 2 l_max - 1 (the magnon starts at l = 1); checked: l_max = 1 vs 2 differ at order 2, l_max = 2 vs 3 differ at order 4, L = 4, 5, 6 agree through their common exact order. Through x⁴ from (L, l_max) = (5, 3). Output `sc_series.csv`.

## Results (detail)

Rotor chain certified intervals: `rotor_certified.csv`, `rotor_certified_large.csv`; per-configuration summary from `summary.py` (smallest g certified and the first failure):

```
L bc lmax dim | smallest certified g: gap in [lo, hi] | next g (fails): floor vs E1_hi | g/L, exp(-2pi/g)/g at smallest g
2 open 4 main 85 | g=0.25: [0.073224, 0.14608] | - | 0.125, 4.86e-11
2 open 8 main 489 | g=0.25: [0.1318, 0.1318] | - | 0.125, 4.86e-11
3 open 5 main 4332 | g=0.35: [0.13109, 0.13204] | g=0.25: floor -3.377 vs E1_hi -5.404 | 0.117, 4.57e-08
3 per 3 main 580 | g=0.5: [0.10625, 0.24573] | g=0.35: floor -4.474 vs E1_hi -5.205 | 0.167, 6.97e-06
3 per 5 main 4332 | g=0.35: [0.12446, 0.12824] | g=0.25: floor -7.028 vs E1_hi -8.672 | 0.117, 4.57e-08
4 open 1 main 70 | g=1.5: [1.0251, 1.2198] | g=1.0: floor 0.7307 vs E1_hi 0.2378 | 0.375, 1.01e-02
4 open 2 main 1107 | g=0.7: [0.10692, 0.37501] | g=0.5: floor -2.19 vs E1_hi -2.363 | 0.175, 1.81e-04
4 open 3 main 8092 | g=0.5: [0.10156, 0.2016] | g=0.35: floor -4.314 vs E1_hi -4.833 | 0.125, 6.97e-06
4 per 1 main 70 | g=1.5: [0.85948, 1.1464] | g=1.0: floor -0.05513 vs E1_hi -0.02489 | 0.375, 1.01e-02
4 per 2 main 1107 | g=1.0: [0.46689, 0.49105] | g=0.7: floor -0.7021 vs E1_hi -1.596 | 0.25, 1.87e-03
4 per 3 main 8092 | g=0.7: [0.23213, 0.2395] | g=0.5: floor -2.277 vs E1_hi -3.691 | 0.175, 1.81e-04
5 open 2 main 8953 | g=1.0: [0.53796, 0.55467] | g=0.7: floor -0.5551 vs E1_hi -1.405 | 0.2, 1.87e-03
5 per 2 main 8953 | g=1.0: [0.43584, 0.46729] | g=0.7: floor -1.77 vs E1_hi -2.019 | 0.2, 1.87e-03
6 open 2 main 73789 | g=3.0: [2.8036, 2.8036] | - | 0.5, 4.10e-02
6 per 2 large 73789 | g=1.0: [0.41289, 0.45673] | g=0.7: floor -2.691 vs E1_hi -2.443 | 0.167, 1.87e-03
8 per 1 large 12870 | g=2.0: [1.5797, 1.6992] | g=1.0: floor -1.93 vs E1_hi -0.6286 | 0.25, 2.16e-02
10 open 1 large 184756 | g=2.0: [1.59, 1.7145] | g=1.0: floor -2.114 vs E1_hi -0.7436 | 0.2, 2.16e-02
10 per 1 large 184756 | g=2.0: [1.5536, 1.7008] | g=1.0: floor -2.637 vs E1_hi -0.925 | 0.2, 2.16e-02
```

Readings:
1. At fixed L the certificate reaches g ≈ 0.25-0.35 for L ≤ 3, 0.5-0.7 for L = 4, 1 for L = 5-6 and 2 for L = 8-10 at the l_max affordable in 3 GB. The failing ingredient is the Q-space floor d: either d < U_1 outright (L = 6 ring g = 0.7, L = 8, 10 at g = 1) or d - U_1 is too small for the coupling term X/(d - U_1), which then pushes the lower bound on λ_1 below U_0 (all L ≤ 5 failures). This is the same floor failure as KS (C5, C9: "floor d drops below E0").
2. The certified numbers are finite-volume rotor-tower gaps. For a ring of L rotors at small g the collective rotor has moment of inertia L/g, so the lowest J = 1 level sits at ~ g/L above the ground state; the certified gaps track this (e.g. 0.13 at L = 3, g = 0.35 vs g/L = 0.12) and are 10⁶ times larger than e^{-2π/g}/g there. To see m(g) one needs L ≫ ξ ~ e^{2π/g}: at g = 0.5, ξ/a ~ 10⁵/C. **Observation, not proof.**
3. AKLT: Knabe/GM turn n ≤ 12 ED into an infinite-volume lower bound 0.3026 vs the true ≈ 0.350 (ring ED). This is the only technique-model pair here that proves an infinite-volume gap, and it uses frustration-freeness and finite local dimension.

## Verdict

**Which failure reasons are the same as for Yang-Mills (C4-C9).**
- Truncation floor extensive / exponential dimension (certified ED): same as C5 and C9 for KS (floor below E0 at g² ≤ 1 on the 2×2 torus).
- Not frustration-free, infinite local Hilbert space (Knabe/GM/martingale): same as KS (literature/03, line 2). The spin-1 Heisenberg chain shows the obstruction is not specific to gauge theory: even a finite-dimensional, gapped (numerically), same-phase model cannot be handled once frustration-freeness is lost, and no stability theorem reaches O(1) deformations.
- Bakry-Émery: curvature (Ric = 1 on S², Ric = 2 on S³ per link for SU(2)) beats Hess(βS) only at β < 1/4 here, β < 1/(4(d-1)) in C6. Same structure. The RK parent has an extensive non-negative correction (C6) and z = 2 dynamics with no relativistic continuum limit (C8). Same.
- Strong-coupling expansion: finite radius in x = 1/g² vs a continuum gap with an essential singularity at x = ∞. Same as C2, Osterwalder-Seiler, Yarotsky.
- Not tested here: Bloch/θ (C4) and Gribov/PW (C7). The O(3) model has θ (Haldane's θ = 2πS), so C4 has a direct analogue: θ = π for half-integer spin gives gaplessness (LSM, cited in literature/03), θ = 0 gives no Bloch-type gap mechanism; not computed.

**Does any technique reach the O(3) continuum gap?** No. None of the four techniques gives a lower bound that survives g → 0 at fixed physical length, on any chain size: certified ED is finite-volume and dominated by the g/L rotor tower; Knabe/GM need frustration-freeness; Bakry-Émery stops at g = 2 and concerns H_FF, not H; the strong-coupling series has a finite radius. **[cited, literature/07]** The known continuum mass gap (m/Λ_MS̄ = 8/e) comes from the Zamolodchikov S-matrix plus TBA (integrability, a bootstrap assumption) matched to perturbation theory, with lattice Monte Carlo confirmation; there is no rigorous proof of a mass gap of the 2D O(3) model at weak coupling. The 2D O(3) model is therefore an exact analogue of the Yang-Mills situation for this toolkit: every electron-level technique tested proves its gap only where Yang-Mills already has one (strong coupling, finite volume, or a frustration-free model in the same phase).

This calibrates the toolkit as follows: a technique that cannot reach weak coupling on the O(3) chain will not reach it in 4D YM. A candidate method for line 1-3 should first be run on this chain as a test.

## Open issues

- Λ_H/Λ_MS̄ for the Hamiltonian O(3) chain not found or computed; it would let the finite-L data be compared with m(g) quantitatively (needs L ≫ ξ, so DMRG/iTEBD, not ED).
- HMN 8/e and the AKLT 0.350 value not re-verified from the primary source (paywalled); AKLT ring value 0.3501 here is our ED.
- Radius of convergence of the strong-coupling series not estimated: 5 exact coefficients only (x⁰..x⁴). Longer series need a linked-cluster code. Hamer-Kogut-Susskind PRD 19 (1979) 3091 has longer O(3) series **[not checked]**.
- All-β positivity of the 1D H_FF gap (classical 1D chain, expected) not proved here; the 2D-lattice H_FF (S = 2D O(3) action) gap at large β is the O(3) analogue of C8 and is open.
- Lanczos completeness and floating-point rounding: certificates are not interval arithmetic.
- `run_rotor.py` main run stopped for time at L = 6 open, g = 3 (L = 6 open, 8, 10 open/periodic at the fine g grid not run); the large-L run covers rings L = 6, 8, 10.

## Reproduce

From repo root (1 thread each; total ≈ 30 min):

    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/run_rotor.py        > .../run.log
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/run_rotor.py large  > .../run_large.log
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/spin1.py
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/sc_series.py
    python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/bakry_emery.py
    .venv/Scripts/python.exe work/20261120-015-line6-rotor-chain/summary.py

No data from `data/` is used (all inputs are model definitions).
