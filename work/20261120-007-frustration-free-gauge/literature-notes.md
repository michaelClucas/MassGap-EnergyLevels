# Literature notes for 20261120-007 (for the overseer to merge into `literature/03` and `literature/02`)

Checked 2026-11-20. Tags: [arXiv-API] = ID, title and authors confirmed via export.arxiv.org; [INSPIRE] = record found via inspirehep.net API; [search] = seen in web search only.

## Citation checks (INSPIRE `refersto:recid:N`)

**arXiv:2409.10612** (Gorantla-Huang, PRB 111, 245110 (2025)), INSPIRE recid 2829234: 12 citing records. All are on non-invertible symmetries, duality defects, kinetically constrained / stochastic models, or tensor networks: 2506.01131, 2410.17317, 2605.12601, 2502.20440, 2512.04174, 2411.05004, 2611.22928, 2504.01258, 2406.12978, 2412.18606, plus 2 without eprints. **None builds a non-Abelian frustration-free gauge model or reuses the martingale gap proof.** Abstract re-read: the FF point has 9 degenerate ground states (one product state plus 3+1D toric-code states), so it is a topological / symmetry-breaking construction, not an RK √(Gibbs) state. [arXiv-API, INSPIRE]

**arXiv:2204.12737** (Shen-Zhu-Zhu), recid 2072947: 19 citing records. All are Euclidean / probabilistic: area law (2505.16585, 2509.04688), U(N) (2510.22788), SO(3) deconfinement (2605.16162), master loop (2202.00880, 2512.00570), 2D YM scaling limits (2602.08591, 2607.22236, 2611.03640), heat-kernel master field (2606.28945), free energy (2511.07297), stochastic quantisation (2201.03487, 2401.13299, 2208.02381, 2306.05166), 2406.19321, 2207.05443, 2603.15770 (formalization of QFT), and the unvetted 2503.15539. **None turns the Langevin / Bakry-Émery gap into a Hamiltonian (RK) statement.** [INSPIRE]

**Dell'Antonio-Zwanziger, "Every gauge orbit passes inside the Gribov horizon", Commun. Math. Phys. 138, 291-299 (1991)**, recid 323367, 200 citing records. Keyword filter (gap, Hamiltonian, spectral, lattice, convex, Poincaré, log-concave) gave mostly gauge-fixing / propagator papers. One is relevant to line 3: **Tedesco, "Some Remarks on the Spectral Geometry of the Gribov Horizon", arXiv:2607.13228** [arXiv-API]. It is a local spectral framework for the Landau-gauge horizon (Morse-Bott Hessian, index theorem), with no Payne-Weinberger / log-concavity gap bound per the abstract. No paper found applying a Payne-Weinberger or log-concave Poincaré bound on the Gribov region. Line-3 novelty survives this check, but only to the depth of the keyword filter.

## Frustration-free / RK constructions relevant to gauge theory

- **Castelnovo, Chamon, Mudry, Pujol, "From quantum mechanics to classical statistical physics: generalized Rokhsar-Kivelson Hamiltonians and the 'Stochastic Matrix Form' decomposition", Ann. Phys. 318, 316 (2005), cond-mat/0502068** [arXiv-API]. SMF-decomposable Hamiltonians correspond one-to-one with classical master equations; ground state = √(Gibbs weight); gap = relaxation rate. This is the general mechanism behind H_FF.
- **Henley, "From classical to quantum dynamics at Rokhsar-Kivelson points", cond-mat/0311345** [arXiv-API]. RK spectrum = classical relaxation spectrum; z = 2 at criticality.
- **Ardonne, Fendley, Fradkin, "Topological order and conformal quantum critical points", Ann. Phys. 310, 493 (2004), cond-mat/0311466** [arXiv-API]. Quantum Lifshitz / RK field theories with ground state e^{-S_cl/2}.
- **Isakov, Fendley, Ludwig, Trebst, "Dynamics at and near conformal quantum critical points", PRB 83, 125114 (2011), 1012.3806** [arXiv-API]. z = 2 dynamics at RK points.
- **Dijkgraaf, Orlando, Reffert, "Relating field theories via stochastic quantization", Nucl. Phys. B 824, 365 (2010), 0903.0732** [arXiv-API]. (d+1)-dim Hamiltonian with ground state e^{-S_d/2} from stochastic quantization.
- **Freedman, Nayak, Shtengel, "A line of critical points in 2+1 dimensions: quantum critical loop gases and non-Abelian gauge theory", PRL 94, 147205 (2005), cond-mat/0408257** [arXiv-API]. SU(2)-type loop-gas RK states.
- **Fendley, Fradkin, "Realizing non-Abelian statistics", PRB 72, 024412 (2005), cond-mat/0502071** [arXiv-API]. RK-form Hamiltonian (sum of projectors) whose ground state is a weighted superposition of Wilson loops (O(n) loop-model weights, n = q + q⁻¹ for SU(2)_q). Per the search snippet, they remark that the trick applies to any field theory with an explicit real action. So **H_FF in this work's README is an instance of a known construction**, not a new idea. It is SU(2) KS-type on L²(SU(2)) links with the Wilson weight, a combination not found written down.
- **Greensite, Olejník, "Dimensional reduction and the Yang-Mills vacuum state in 2+1 dimensions", PRD 77, 065003 (2008), 0707.2860** [arXiv-API]. Vacuum wavefunctional ansatz; related to e^{-S_{d}}-type (dimensional-reduction) vacua. Follow-ups 1011.5046, 1102.3941 [search]. Physics context for "ground state ≈ √(lower-dim YM measure)". These are not FF Hamiltonians.

## Exactly solvable FF (commuting-projector) gauge models

- **Kitaev, "Fault-tolerant quantum computation by anyons", Ann. Phys. 303, 2 (2003), quant-ph/9707021** [arXiv-API]. Quantum double D(G) for finite G, including non-Abelian G: commuting projectors, gap exactly 1, deconfined topological phase.
- **Bombin, Martin-Delgado, "A family of non-Abelian Kitaev models on a lattice: topological confinement and condensation", PRB 78, 115421 (2008), 0712.0190** [arXiv-API]. Non-Abelian quantum doubles with extra terms.
- **Levin, Wen, "String-net condensation", PRB 71, 045110 (2005), cond-mat/0404617** [arXiv-API]. SU(2)_k string nets: commuting projectors, exactly gapped, topological.
- **Zache, González-Cuadra, Zoller, PRL 131, 171902 (2023), 2304.02527** and **Hayata, Hidaka, 2305.05950** [arXiv-API]. q-deformed KS = string-net magnetic projectors plus the KS electric term. The electric term breaks the commuting-projector structure. Hayata-Hidaka study scars, not gaps.
- **Bravyi, Hastings, Verstraete, PRL 97, 050401 (2006), quant-ph/0603121** [arXiv-API]. Topological order cannot be generated by local evolution in short time, so the string-net FF points are separated from the trivial confined phase.

## Finite-size criteria (beyond those already in `literature/03`)

- **Lemm, "Finite-size criteria for spectral gaps in D-dimensional quantum spin systems", 1902.07141** [arXiv-API]. D-dim Knabe-type criteria with deformed (weighted) patches.
- **Anshu, "Improved local spectral gap thresholds for lattices of finite dimension", PRB 101, 165104 (2020), 1909.01516** [arXiv-API]. Improved thresholds via detectability-lemma ideas.
- **Aharonov, Arad, Landau, Vazirani, "The detectability lemma and quantum gap amplification", 0811.3412** [arXiv-API].
- These proofs are operator inequalities on projectors. The Knabe argument written in the README uses no finite-dimensionality. The finite-dimensional assumption in these papers enters only through numerics and through the stability theorems.

## Glauber / Langevin dynamics for lattice gauge theories

- **Malacan, "Slow and rapid mixing in Z2 hyperlattice gauge theory", 2611.19051** [arXiv-API]. Heat-bath Glauber for Z2 gauge theory: Dobrushin condition (rapid mixing) at small β, exponentially slow at large β. Abelian finite group only.
- **Yin, Khemani, "Fast mixing and thermal phases from strong decay of correlations", 2611.17433** [arXiv-API]. General: strong correlation decay ⇒ constant Glauber gap (abstract per search). Possibly a route to criterion-K-like statements for π_β; not checked for continuous spins or gauge symmetry.
- No paper found on Glauber or Langevin spectral gaps for SU(N) Wilson measures beyond Shen-Zhu-Zhu's Bakry-Émery regime. Searched 2026-11-20; absence of evidence only.
- Aizenman-Holley (1987) and L. Wu (Ann. Probab. 2006), "Dobrushin ⇒ Glauber gap / Poincaré": recalled from memory, **unverified**.

## Yarotsky (resolves a STATE.md open question)

**Yarotsky, "Ground states in relatively bounded quantum perturbations of classical lattice systems", CMP 261, 799 (2006), math-ph/0412040** [arXiv-API; hypotheses read from the PDF, pp. 2-4].
- Site spaces may be infinite-dimensional; tensor-product structure on Z^ν; translation-invariant interactions of finite range Λ_0.
- H_0 = Σ h_x "classical": a function of an orthogonal partition of unity containing the projector onto Ω_x; Ω non-degenerate ground state with h_x ≥ 1 on its complement (eq. 1).
- Perturbation φ_x form-bounded: |φ_x(v,v)| ≤ α‖h_x^{1/2}v‖² + β‖v‖² (eq. 2). Theorem 1: small α, β ⇒ unique gapped ground state uniformly in volume, a weak* infinite-volume limit, exponential clustering, analyticity. Theorem 2: any α < 1 if the bounded part is small enough.
- KS with Gauss law: apply on the unconstrained tensor-product space (E² diagonal in isotypic projectors, plaquette term bounded). Gauge invariance of the unique ground state follows because SU(2)^V has no nontrivial characters. **The strong-coupling KS gap with Gauss law is covered.** Constants are not explicit.
