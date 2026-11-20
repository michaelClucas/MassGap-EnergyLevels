# 02. Rigorous lattice and constructive results

Cross-cuts lines of attack 1, 2, 3, 6. Tags as in `README.md`.

## Classical results

**Osterwalder, Seiler, "Gauge field theories on a lattice", Ann. Phys. 110, 440 (1978).** DOI 10.1016/0003-4916(78)90039-8 [Crossref]
- Proves reflection (physical) positivity of Wilson's lattice action, a transfer matrix, and at small β convergent cluster/polymer expansion: analyticity, exponential clustering (mass gap of the transfer matrix), and Wilson area law (confinement).
- Setting: Euclidean lattice, any d, compact gauge group, infinite volume, strong coupling only.
- Relevance: the rigorous strong-coupling gap. Any Hamiltonian gap proof at strong coupling reproves a version of this; novelty needs extending the coupling range.

**Balaban, "Renormalization group approach to lattice gauge field theories" I, Commun. Math. Phys. 109, 249 (1987)** DOI 10.1007/BF01215223 [Crossref], **II (cluster expansions), CMP 116, 1 (1988)** DOI 10.1007/BF01239022 [Crossref]. Earlier ultraviolet-stability papers in CMP 1984-1985 [unverified individually].
- Multiscale RG for lattice YM in 3D and 4D at weak coupling: ultraviolet stability bounds uniform in lattice spacing for finite volume.
- Setting: Euclidean lattice, d=3,4, SU(N)/compact group, finite volume, weak coupling.
- Relevance: the only serious control of the continuum limit (UV) for 4D YM. Does not give a gap or infinite volume. The infrared problem is exactly what this program is about.

**Magnen, Rivasseau, Sénéor, "Construction of YM4 with an infrared cutoff", Commun. Math. Phys. 155, 325 (1993).** DOI 10.1007/BF02097397 [Crossref]
- Continuum construction of 4D YM Schwinger functions with an infrared cutoff (finite volume / regulator), in a regularized gauge.
- Setting: continuum, 4D, SU(2)-type, IR-cutoff.
- Relevance: UV is (largely) under control in finite volume; the open problem is the IR (gap, infinite volume).

**Göpfert, Mack, "Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant", Commun. Math. Phys. 82, 545 (1982).** DOI 10.1007/BF01961240 [Crossref]
- Rigorous version of Polyakov's monopole mechanism: confinement and a mass gap at all couplings for 3D Villain U(1).
- Setting: Euclidean lattice, 3D, U(1), infinite volume, all couplings.
- Relevance: the only "all couplings" rigorous gap for a 3D gauge theory. Stepping stone; mechanism is a Debye-screened monopole plasma.

**Yarotsky, "Ground states in relatively bounded quantum perturbations of classical lattice systems", Commun. Math. Phys. 261, 799 (2006).** arXiv:math-ph/0412040 [arXiv-API]
- Uniform spectral gap, unique ground state and clustering for small relatively bounded quantum perturbations of classical lattice Hamiltonians, allowing infinite-dimensional local Hilbert spaces.
- Setting: quantum lattice systems, any d, infinite volume.
- Relevance: the natural tool to prove the Kogut-Susskind gap at strong coupling (electric term = "classical" part, plaquette term = perturbation). Whether it applies directly to KS with the Gauss-law constraint has not been checked by this team. Treat a strong-coupling Hamiltonian gap as essentially known.

## Probabilistic / stochastic-analysis results (2015-2026)

**Chatterjee, "Rigorous solution of strongly coupled SO(N) lattice gauge theory in the large N limit", Commun. Math. Phys. (2019).** arXiv:1502.07719 [arXiv-API; journal ref from memory, unverified]
- Wilson loop expectations at strong coupling and large N via master loop (Makeenko-Migdal type) equations; sum over random surfaces.
- Setting: Euclidean lattice, any d, SO(N), N→∞ at fixed small 't Hooft coupling.

**Chatterjee, "Yang-Mills for probabilists".** arXiv:1803.01950 [arXiv-API]
- Survey of the problem for probabilists; states lattice versions of the gap question.

**Chatterjee, "A probabilistic mechanism for quark confinement", Commun. Math. Phys. (2021).** arXiv:2006.16229 [arXiv-API]
- If a pure lattice gauge theory (compact connected group) has exponential decay of correlations under arbitrary boundary conditions, and unbroken center symmetry, then Wilson loops obey a (weak) confinement law. Mass gap ⇒ confinement.
- Setting: Euclidean lattice, d≥2, compact connected group.
- Relevance: links "gap" to "confinement" rigorously in one direction.

**Adhikari, Cao, "Correlation decay for finite lattice gauge theories at weak coupling".** arXiv:2202.10375 [arXiv-API]
- Exponential decay of correlations (a gap) at weak coupling for finite (possibly non-Abelian) gauge groups.
- Setting: Euclidean lattice, d=4 (and general), finite group, infinite volume, weak coupling.
- Relevance: weak-coupling gap for finite groups. Finite groups have a gapped vacuum for trivial reasons at weak coupling (frozen discrete link variables); does not extend to continuous groups.

**Shen, Zhu, Zhu, "A stochastic analysis approach to lattice Yang-Mills at strong coupling", Commun. Math. Phys. (2023).** arXiv:2204.12737 [arXiv-API; journal ref unverified]
- For |β| < 1/(16(d−1)) (SU(N)) or (N−2)/(32(d−1)N) (SO(N)) with 't Hooft scaling: Langevin dynamics on SU(N)^links, **Bakry-Émery criterion** (Ricci curvature of SU(N) dominates the Hessian of the Wilson action), log-Sobolev and Poincaré inequalities, unique infinite-volume limit, large-N factorization, and a **mass gap**.
- Setting: Euclidean lattice, any d≥2, SU(N)/SO(N), infinite volume, strong coupling (uniform in N).
- Relevance: **prior art for "Bakry-Émery curvature gives the YM gap"** on the lattice. Proven only where the curvature of the group beats the action, i.e. strong coupling.

**Shen, Zhu, Zhu, "Langevin dynamics of lattice Yang-Mills-Higgs and applications".** arXiv:2401.13299 [arXiv-API]
- Extends the above to YM-Higgs.

**Cao, Park, Sheffield, "Random surfaces and lattice Yang-Mills".** arXiv:2307.06790 [arXiv-API]
- Wilson loop expectations as sums over embedded random surfaces for U(N), SO(N), Sp(N); strong-coupling and large-N analysis.

**Cao, Nissim, Sheffield, "Expanded regimes of area law for lattice Yang-Mills theories".** arXiv:2505.16585 [arXiv-API]
- Extends the β range where area law is proved for U(N), especially at large N, improving Osterwalder-Seiler, via an a priori bound for the master loop equation.

**Cao, Nissim, Sheffield, "Dynamical approach to area law for lattice Yang-Mills".** arXiv:2509.04688 [arXiv-API]
- Uses the Shen-Zhu-Zhu mass gap plus Durhuus-Fröhlich (1980) to get area law in the 't Hooft regime for U(N), SU(N), SO(2N).

**Nissim, "U(N) lattice Yang-Mills in the 't Hooft regime".** arXiv:2510.22788 [arXiv-API]
- Mass gap, unique infinite-volume limit, large-N limit for U(N), where Bakry-Émery fails because U(N) Ricci curvature is not uniformly positive (U(1) factor is flat); handled by a random-environment SU(N) model plus cluster expansion.
- Relevance: shows exactly where the curvature method breaks (flat directions). Directly informs a Result-C obstruction for the Bakry-Émery route.

**Chandra, Chevyrev, Hairer, Shen, "Langevin dynamic for the 2D Yang-Mills measure", Publ. Math. IHÉS 136, 1 (2022).** arXiv:2006.04987 [arXiv-API]
**Chandra, Chevyrev, Hairer, Shen, "Stochastic quantisation of Yang-Mills-Higgs in 3D", Invent. Math. 237, 541 (2024).** arXiv:2201.03487 [arXiv-API]
- Continuum stochastic quantization: gauge-covariant state spaces for 2D YM and 3D YMH (local in time).
- Setting: continuum, d=2,3, finite volume (torus).
- Relevance: current frontier of continuum construction; no gap results.

**Chatterjee, "A scaling limit of SU(2) lattice Yang-Mills-Higgs theory", Prob. Math. Phys. 7, 339 (2026).** arXiv:2401.10507 [arXiv-API]
- A nontrivial scaling limit of 4D lattice YMH in a particular regime.

**Piazza, Shen, "Log-Sobolev inequality for 2D SU(N) lattice gauge theory at all couplings".** arXiv:2610.35286 [arXiv-API]
- LSI with volume-uniform constant at any fixed β in 2D. 2D gap itself is elementary.

## Unvetted claims (do not cite as results)

- Jacobsen, "A Constructive Proof of Existence and Mass Gap for Pure SU(3) Yang-Mills in Four-Dimensional Space-Time", arXiv:2506.00284 [arXiv-API], physics.gen-ph. Claims a full solution; no independent acceptance found.
- Gauvin, arXiv:2503.15539 [arXiv-API]; abstract states no continuum contraction is established.
- Mihajlovic, "SU(2) Yang-Mills: Spectral Bounds, Vacuum Concentration, and Criteria for the Continuum Limit", arXiv:2611.11534 [arXiv-API], November 2026. Explicitly partial: finite-cutoff transfer-matrix bounds, fidelity obstructions for hard electric cutoffs. Worth reading for the hard-cutoff obstruction claims, but unrefereed.
- Several Zenodo/ResearchGate/Academia records claiming mass-gap proofs via Bakry-Émery or "mobility gap" arguments. Not peer reviewed; ignore unless a specific lemma is needed.
