# 20261120-014: Line 5, interaction-driven gaps (BCS, Mott, Peierls, excitonic) as a route to the YM mass gap

## Question

Line of attack 5: can the mechanisms that open interaction-driven gaps in electron systems (BCS pairing, Mott localisation, Peierls dimerisation, excitonic condensation) and their published Yang-Mills transcriptions (gluon BCS vacuum, Cornwall dynamical mass, <A^2> condensate, dual superconductor, Savvidy / Nielsen-Olesen, large-N mean field) give a controlled lower bound on the Yang-Mills spectral gap? For each mechanism: why it is controlled in condensed matter, the YM version, where error control is lost (shown by explicit calculation or counterexample), and whether any version has error control.

## Verdict (read this first)

**Result-C obstruction for all four mechanisms as routes to a gap lower bound in 4D (or 3D) continuum SU(N).** Every controlled condensed-matter gap equation gets its control from one of three things, and pure YM has none of them at the scale of the gap:

1. a **Fermi surface** (BCS, Peierls, excitonic): the pairing susceptibility has a 1D log, so the gap is exp(-1/lambda) with a free small lambda. Gluons are bosons; there is no Fermi surface and no log (section B, G1).
2. a **free small/large parameter independent of the gap scale** (Mott: t/U; BCS: lambda, omega/E_F; vector large N: 1/N). In YM, m = kappa Lambda with kappa a fixed pure number, so g^2(m) is fixed (section C). The Mott parameter is 1/g^2, which is exactly the strong-coupling regime already covered (Osterwalder-Seiler, C2) and goes the wrong way in the continuum limit.
3. a **classical / static field** (Peierls with adiabatic phonons, Ginzburg-Landau). The YM background fields (Savvidy) are not classical and their fluctuation operator has a tachyon (section G2).

Two further failures apply to every YM transcription:
- **The quantity solved for is not a spectral gap.** Covariant-gauge "gluon masses" from Hartree/gap equations depend on the gauge parameter (section E); Gribov-type propagators violate reflection positivity (section F); so no Hamiltonian eigenvalue is being bounded.
- **Mean-field is not a bound on E1 - E0.** Variational vacua (BCS-type, Gaussian) give upper bounds on E0 only. The resulting gap can be wrong by any factor in either qualitative direction (section D), and a leading-order gap equation can predict a gap where a theorem says there is none (section A, 2D XY).

**Only controlled versions found:** (i) vector large-N gap equations (exact at N = inf, rigorous 1/N in 2D: Kopper-Magnen-Rivasseau for Gross-Neveu); (ii) Polyakov/Göpfert-Mack 3D compact U(1) (Debye screening of a dilute monopole gas); (iii) strong-coupling lattice ("Mott") expansions. None reaches continuum non-Abelian YM: (i) is a vector model (YM is a matrix model, planar limit not a saddle point), (ii) needs a small fugacity e^{-S_0} that pure YM lacks, (iii) needs g^2 large.

**Open sub-question worth pursuing (conjectured to be open, literature check needed):** a rigorous large-N gap-equation proof for the **2D CP^{N-1} model**, the simplest asymptotically free model with a gauge field where the gap equation is exact at N = inf and the gap is a gauge-invariant spectral quantity (bound state of z-bar z). This would be the "BCS-with-gauge-field" stepping stone; it feeds line 6, not 4D YM. See Open issues.

## Setting conventions

"Lattice d" = Euclidean Z^d; "continuum" = R^d Euclidean or (d-1)+1 Minkowski; all statements are infinite volume unless stated. Status tags: **[proved]** (here or by elementary algebra shown in the script), **[cited]** (published; reference tagged [Crossref] was checked against the Crossref API on 2026-11-20 for authors, title, journal, volume, pages; the theorem content is from memory unless stated), **[numerical]**, **[conjectured]**.

## Method

`interaction_gaps.py` (sections A-G), output in `run.log`. Run from repo root:
`python tools/pin.py .venv/Scripts/python.exe -u work/20261120-014-line5-interaction-gaps/interaction_gaps.py` (about 1 min, 1 core, < 200 MB). No lattice data needed; none hardcoded.

## Explicit calculations

### A. Large-N gap equation vs exact spectrum; N = 2 counterexample [proved + cited]
Setting: classical lattice O(N) model on Z^d (= Euclidean (d-1)+1 quantum field on a lattice), infinite volume.
N = inf (spherical model, Stanley 1968 Phys. Rev. 176, 718 [Crossref]) gap equation: beta = W_d(m), W_d(m) = int d^dk/(2pi)^d 1/(m^2 + k-hat^2). Computed: W_3(0) = 0.252731, W_4(0) = 0.154933 (Watson integrals), so for d >= 3 there is a gapless phase beta > W_d(0); for d = 2, W_2 diverges logarithmically and m ~ sqrt(32) e^{-2 pi beta} > 0 for every beta (table in `run.log`, matches the asymptotic form to 1e-7 at beta = 1.5).
**Counterexample.** Applying the leading-order gap equation to N = 2 (2D XY, beta_sph = beta_XY/2) predicts m > 0 for every beta_XY (e.g. m = 1.06e-2 at beta_XY = 2). Fröhlich-Spencer, Commun. Math. Phys. 81, 527 (1981) [Crossref] prove power-law decay of <s_0 . s_x> at large beta_XY, i.e. m = 0. So "gap equation has a massive solution" does not imply a gap. (The numerical BKT point beta_XY ~ 1.12, Hasenbusch J. Phys. A 38, 5869 (2005) [Crossref; value from memory], is not used.) The leading-order equation does not know that the one-loop beta-function coefficient is proportional to N - 2; only the 1/N corrections see it. The control is the 1/N expansion, not the gap equation.

### B. Fermi-surface log vs bosonic Hartree [proved]
BCS with constant density of states and cutoff omega: 1 = lambda asinh(omega/Delta), Delta = omega / sinh(1/lambda) ~ 2 omega e^{-1/lambda}: a solution exists for every lambda > 0 because the Fermi surface makes the pair susceptibility a 1D integral int dxi/|xi|. Bosonic Hartree with massless bare field, Z^4: m^2 = lambda W_4(m) gives m = sqrt(lambda W_4(0)) (1 + o(1)) in lattice units, a cutoff-scale tadpole (lambda = 0.1: m = 0.124 vs BCS 9.1e-5). This is a mass renormalisation, cancelled by a counterterm; in a gauge theory it is forbidden by BRST and is zero in dimensional regularisation. A bosonic IR log appears only for d = 2 (section A). In 4D YM the only log is the UV running of g, which gives dimensional transmutation but no small parameter at the gap scale (C).

### C. Coupling at the gap scale [proved, one loop]
Continuum 4D SU(N). Any RG-invariant gap is m = kappa Lambda with kappa a pure number. One loop: a(mu) = alpha N/(4 pi) = 3/(11 ln(mu^2/Lambda^2)), so a(m) = 3/(11 ln kappa^2): 0.098 at kappa = 4, 0.050 at kappa = 15. Unlike BCS (lambda free, Delta/omega -> 0 as lambda -> 0 with relative corrections O(lambda)), nothing makes a(m) -> 0. Getting kappa for the 0++ requires Lambda_MS-bar r0, not in `data/` (open issue); this does not change the argument.

### D. Gaussian (bosonic BCS / HFB) variational gap is not a bound [numerical, exact FD to ~1e-4 abs.]
Setting: 1D QM H = p^2/2 + mu x^2/2 + lambda x^4 (zero-dimensional field theory; the Gaussian squeezed state is the bosonic analogue of the BCS vacuum and of the Szczepaniak-Swanson Coulomb-gauge vacuum).

| mu | lambda | gap exact | Gaussian MF gap (symmetric) | broken MF | ratio sym/exact |
|---|---|---|---|---|---|
| 1 | 0.1 | 1.2104 | 1.2212 | no | 1.009 |
| 0 | 1 | 1.7257 | 1.8171 | no | 1.053 |
| -1 | 0.1 | 0.2969 | 0.4855 | no | 1.64 |
| -2 | 0.1 | 1.51e-2 | 0.288 | yes, predicts 0 | 19 |
| -4 | 0.1 | 2.08e-7 | 0.149 | yes, predicts 0 | 7e5 |

E0(MF) >= E0 always (Rayleigh-Ritz) [proved], but the gap ratio is unbounded above, and the best (broken) MF predicts degeneracy where the exact gap is > 0. In field theory the symmetric-MF failure is the dangerous one: an O(1) mean-field gap can sit over a true gap that vanishes in the limit (as in A). Gap = difference of two energies; an upper bound on each gives no bound on the difference [proved, trivial].

### E. Gauge-parameter dependence of the Hartree gluon mass [proved, standard one-loop algebra]
Continuum 4D SU(N), linear covariant gauge xi, cutoff regularisation. The four-gluon-vertex tadpole at p = 0 has Lorentz structure g_mn D_rr - D_mn with D_rs = [(delta - kk/k^2) + xi kk/k^2]/k^2, giving Pi_mn = c g^2 C_A I_Lambda (d-1+xi)(d-1)/d g_mn (verified numerically in `run.log`; overall constant c not tracked). The Hartree "mass" therefore scales as (d-1+xi): Feynman/Landau = 4/3. A genuine pole mass is xi-independent (Nielsen, Nucl. Phys. B 101, 173 (1975) [Crossref]), so this quantity is not one. The same applies to <A^2> itself, = (d-1+xi) x const in the free theory; Gubarev-Stodolsky-Zakharov (hep-ph/0010057) fix this by minimising over the orbit, which selects Landau gauge, and the LCO effective potential (Verschelde et al., hep-th/0105018) is computed in Landau gauge. Pinch-technique masses (Cornwall 1982; Aguilar-Binosi-Papavassiliou 1511.08361) are xi-independent by construction, but are momentum-dependent m(q^2), not poles, and rest on truncated Schwinger-Dyson equations with no error bound [cited, literature/06].

### F. Gribov propagator violates positivity [proved]
D(p) = p^2/(p^4 + gamma^4). Zero-spatial-momentum Schwinger function C(t) = (gamma/2 sqrt2) e^{-a}(cos a - sin a), a = gamma t/sqrt2 (closed form checked by quadrature to 1e-7). C(t) < 0 for a in (pi/4, 5pi/4). A positive Källén-Lehmann measure gives C(t) > 0, so gamma is not an energy gap of any positive-metric Hamiltonian. Lattice Landau-gauge propagators also violate positivity: Cucchieri-Mendes-Taurines PRD 71, 051902 (2005); Bowman et al. PRD 76, 094505 (2007) [both Crossref]. Kugo-Ojima, Prog. Theor. Phys. Suppl. 66, 1 (1979) [Crossref]: gluons are outside the physical (BRST-cohomology) space in the confining scenario. So any gap equation for the gluon propagator is solving for a non-spectral quantity.

### G. Peierls and Landau levels: fermion vs boson [proved]
G1. Dimerised ring, hoppings t(1 +- delta). Half-filled spinless fermions: band-energy gain per site fits a delta^2 ln(1/delta) + b delta^2 with a = 0.3181 (1/pi = 0.3183), which beats any elastic cost K delta^2/2: dimerisation at every K (rigorous version: Kennedy-Lieb, PRL 59, 1309 (1987); Lieb-Nachtergaele, PRB 51, 4777 (1995) [both Crossref]; phonons classical/static). Free bosons condense at the band bottom, energy -(t1 + t2) = -2t per particle independent of delta: gain exactly 0, no Peierls instability at any K.
G2. Landau levels, E^2 = (2n+1)B - 2 B s_z (g = 2). Dirac electron (s_z = +-1/2): min E^2 = 0. Charged spin-1 (off-diagonal gluon in a Cartan background): n = 0, s_z = +1 gives E^2 = k_z^2 - B < 0, the Nielsen-Olesen tachyon (Nucl. Phys. B 144, 376 (1978), literature/06). Fermions in a background field are bounded below; gluons are not, so the Savvidy "condensate" is not a stable mean-field vacuum and its one-loop energy is complex.

## Per-mechanism assessment

### 1. BCS
(a) Condensed matter: Cooper instability from the Fermi-surface log; Delta = 2 omega e^{-1/lambda}. Rigorous for the BCS functional (mean-field model, not the many-body Hamiltonian): Frank-Hainzl-Naboko-Seiringer, J. Geom. Anal. 17, 559 (2007); Hainzl-Seiringer, PRB 77, 184517 (2008) [both Crossref]; GL limit Frank-Hainzl-Seiringer-Solovej, JAMS 25, 667 (2012) [Crossref]. The reduced pairing Hamiltonian is exactly solvable (Richardson, Phys. Lett. 3, 277 (1963) [Crossref]) and BCS is exact for it as the number of levels -> inf [cited, from memory]. Derivation of BCS from a full many-body Hamiltonian with Coulomb repulsion is not rigorous [cited, standard].
(b) YM: Coulomb-gauge BCS vacuum + TDA/RPA glueballs (Szczepaniak-Swanson-Ji-Cotanch hep-ph/9511422; Szczepaniak-Swanson hep-ph/0107078); LCO <A^2> gap equation (hep-th/0105018), structurally Gross-Neveu.
(c) Lost: no Fermi surface (B, G1); no small parameter at the gap (C); gap equation for a gauge-dependent / non-positive object (E, F); variational vacuum bounds E0 only (D). In Coulomb gauge the Christ-Lee Hamiltonian is truncated (Coulomb kernel approximated, FP determinant partially kept) with no error estimate [cited, literature/06 summary].
(d) Controlled version: Gross-Neveu at large N in 2D, Kopper-Magnen-Rivasseau, Commun. Math. Phys. 169, 121 (1995) [Crossref] (rigorous mass generation; vector model, fermions). No matrix-model analogue known. **Obstruction.**

### 2. Mott
(a) Hubbard at large U/t: charge gap ~ U, controlled by t/U via Schrieffer-Wolff with rigorous bounds (Bravyi-DiVincenzo-Loss, Ann. Phys. 326, 2793 (2011) [Crossref]) and quantum Pirogov-Sinai (Datta-Fernández-Fröhlich, J. Stat. Phys. 84, 455 (1996) [Crossref]); 1D gap for all U > 0 from the Bethe ansatz (Lieb-Wu, PRL 20, 1445 (1968) [Crossref]).
(b) YM: the Kogut-Susskind Hamiltonian at large g^2 is the Mott limit: electric (on-link Casimir) term = U, plaquette = hopping. Gap 3g^2/2 + ... (C2), rigorous strong-coupling gap (Osterwalder-Seiler; Yarotsky with Gauss law; literature/02).
(c) Lost: the control parameter 1/g^2 must go to infinity in the continuum limit (asymptotic freedom), the opposite end. On small lattices this program found the strong-coupling series matching ED at g^2 = 3, 5, 10 (C2) and a geometry-dominated harmonic gap at weak coupling (C3), with no continuum content in between. Even if no phase transition separates the regimes (expected for 4D SU(N) [conjectured, standard lore]), the strong-coupling expansion has a finite radius of convergence and does not reach g -> 0.
(d) Controlled only at strong coupling, already known. **Obstruction (same as the strong-coupling limitation in STATE).**

### 3. Peierls
(a) Static classical lattice distortion + free fermions: E(delta) contains -delta^2 ln(1/delta) from 2k_F nesting (G1). Rigorous: Kennedy-Lieb 1987, Lieb-Nachtergaele 1995 (period-2 minimiser for half filling, classical phonons).
(b) YM: translation-symmetry-breaking or background-field vacua, Savvidy constant field and the Nielsen-Olesen "spaghetti" domain picture; nothing breaks translations in the YM vacuum (lattice data show no such order [numerical lore]).
(c) Lost: the control comes from (i) fermions with a nested Fermi surface and (ii) a classical static field. YM has bosons (no log: G1 gives zero gain) and a quantum background whose fluctuation operator is unbounded below (G2).
(d) None. **Obstruction.**

### 4. Excitonic insulator
(a) Jérome-Rice-Kohn, Phys. Rev. 158, 462 (1967) [Crossref]: electron and hole pockets with nested Fermi surfaces pair via Coulomb attraction; BCS-type gap equation, same log. In a semiconductor, excitons are hydrogenic bound states below an existing band gap.
(b) YM: glueballs as bound states of constituent (massive) gluons, TDA/RPA on the BCS vacuum (hep-ph/9511422); two-gluon Schrödinger models (literature/08).
(c) Lost: either needs a Fermi surface (no) or presupposes a one-particle gap (the gluon mass, which is the thing to be proved and is not spectral, E/F). A bound-state binding energy gives an upper bound on the lowest two-body level, never a lower bound on the gap.
(d) None as a proof route; useful as phenomenology (Result B task). **Obstruction.**

### 5. Dual superconductor (cross-cut)
Controlled: 3D compact U(1)/Georgi-Glashow (Polyakov 1977; Göpfert-Mack, literature/02): photon mass^2 ~ e^{-S_0}, a Debye-screening (electron-gas) gap, controlled by the small monopole fugacity, S_0 ~ 4 pi m_W/g^2 large because the Higgs breaks SU(2) -> U(1) at weak coupling. Seiberg-Witten (hep-th/9407087): controlled by SUSY. Pure YM is the m_W -> 0 limit: S_0 -> 0, dilute-gas control lost; abelian projection is gauge dependent ('t Hooft 1981). Ginzburg-Landau control (FHSS 2012) needs a microscopic BCS model, absent. **Obstruction.**

### 6. Large-N mean field (cross-cut)
Vector models: gap equation exact at N = inf (Stanley 1968), 1/N rigorous in some cases (Kupiainen, Commun. Math. Phys. 73, 273 (1980) [Crossref]; Kopper-Magnen-Rivasseau 1995). Even there, applying the N = inf equation at finite N can be qualitatively wrong (A). YM: adjoint (matrix) fields; the N = inf limit is the planar sum, governed by loop equations, not by a saddle point of a finite set of collective fields with N in front of the action. The one exception with an electron reduction is single-matrix QM, whose singlet sector is N free fermions (Brézin-Itzykson-Parisi-Zuber, Commun. Math. Phys. 59, 35 (1978) [Crossref]); YM has many non-commuting matrices, and no such reduction is known [cited/conjectured]. Rigorous large-N YM results are strong-coupling only (literature/02). **Obstruction for 4D; no lower bound obtained.**

## Which version gives a gap lower bound?
None for continuum YM. A variational (BCS/Gaussian/Coulomb-gauge) calculation gives a rigorous upper bound on E0, and a variational upper bound on E1 only in a symmetry sector orthogonal to the ground state; E1 - E0 is bounded by neither (D). Lower bounds come only from: strong-coupling expansions (Mott, known), vector large N (not YM), dilute monopole gas (needs Higgs). [proved for the logical point; cited for the three exceptions]

## Open issues
- **CP^{N-1} sub-question [conjectured open]:** is there a rigorous 1/N proof (Kupiainen / KMR style) of mass generation in 2D CP^{N-1} at weak coupling, all large N? Needs a literature check before any work: D'Adda-Lüscher-Di Vecchia 1978 and Witten 1979 give the N = inf solution (not verified here). If open and tractable, it is the cleanest "gap equation + gauge field + controlled" stepping stone for line 6.
- Lambda_MS-bar r0 for pure SU(3) should be added to `data/` (with source) if section C is to quote kappa for the 0++.
- Section E covers the cutoff Hartree truncation only; the full one-loop gauge dependence of the LCO <A^2> mass in linear covariant gauges (Dudal et al. papers) was not checked.
- The 2D XY counterexample uses Fröhlich-Spencer qualitatively (large beta); no explicit beta_0 is claimed.
- Section D FD accuracy: absolute ~1e-4; the 2e-7 splitting at mu = -4 is order-of-magnitude only.
- Theorem statements of cited math papers (FHNS, HS, KMR, DFF, BDL) are from memory; only bibliographic data were verified.
