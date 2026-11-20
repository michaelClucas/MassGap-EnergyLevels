# 20261120-001 Literature survey

## Question

For each line of attack in `OVERSEER_PROMPT.md`: what prior work exists, what is already known, what is open, and which unexplored angle is most promising? Which proposed tasks would only rediscover published work?

## Method

- Web search (standard and extended) for each topic, 2026-11-20.
- Every arXiv ID checked against the arXiv API (title, authors, journal ref). Every DOI for pre-arXiv papers checked against the Crossref API. Unchecked items are tagged "unverified" in `literature/`.
- Abstracts fetched via the arXiv API for recent papers (2015-2026) to write the summaries. Classic papers summarized from standard knowledge; flagged where a summary was not checked against the source.
- `verify_ids.py` re-checks all cited arXiv IDs. Run: `python tools/pin.py python work/20261120-001-literature-survey/verify_ids.py`. Last run: `verify_ids.log` (all cited IDs returned by the API).
- No compute beyond API calls.

Output: `literature/README.md` (index and prior-art list) and `literature/01-...08-*.md`.

## Result

Status labels: (cited) = stated in a cited source; (conj) = this survey's judgment, not proved.

### Line 1. Hamiltonian lattice YM as a many-body problem

- Known (cited): KS Hamiltonian is a quantum rotor (group-manifold Laplacian) per link plus plaquette hopping. Strong-coupling gap series (Hamer-Irving 1985, 2+1D SU(2)), t-expansion (Horn-Weinstein), **coupled cluster** (Schuette, Wichmann et al. 1998-2002), analytic variational SU(2)-SU(5) (Carlsson-McKellar 2003), ED/TN of truncated SU(2) and SU(3) plaquette chains and small 2D lattices (Klco 2019, Ciavarella 2021, Yao 2023, Ebner 2024, Hayata-Hidaka 2023), gauge-invariant LSH basis (Raychowdhury-Stryker 2020), q-deformed KS (Zache et al. 2023), large-N_c Krylov truncations (Ciavarella-Burbano-Bauer 2025). Qubit-regularized SU(2) plaquette chain mapped to TFIM with a continuum limit and E8 "glueballs", √σ/m₁ = 0.249(1) (Siew-Chandrasekharan-Bhattacharya 2025/26).
- Open: no Hamiltonian calculation with controlled truncation error reaches the 3+1D scaling region; small-lattice ED gaps are far from the continuum.
- Already done, do not claim novelty: ED gap vs coupling on small SU(2) lattices; CC and variational glueball masses; plaquette chain as a gapped spin chain.
- Most promising unexplored angle (conj): **certified two-sided gap bounds** on small SU(2) KS lattices in the LSH basis: Temple/Weinstein lower bounds plus an explicit bound on the j_max truncation error. Found no paper giving certified (not just converged) KS gap intervals. This feeds line 2 (finite-size criteria need certified finite-patch gaps).

### Line 2. Spin-system spectral-gap methods

- Known (cited): AKLT, Knabe, Nachtergaele martingale, Gosset-Mozgunov and Lemm-Mozgunov finite-size thresholds, computer-assisted 2D AKLT gap proofs (Lemm-Sandvik-Wang, Pomata-Wei 2020), SDP hierarchy of certificates (Rai et al. 2026), stability of frustration-free gaps (Michalakis-Pytel, Bravyi-Hastings-Michalakis). **Gorantla-Huang 2025 proved a gap with the martingale method for a frustration-free deformation of 3+1D Z2 lattice gauge theory.** Yarotsky 2006 gives gaps for small quantum perturbations of classical lattice systems with infinite-dimensional local spaces, the natural strong-coupling KS tool.
- Structural obstruction (cited/standard): KS is not frustration-free (electric and magnetic terms have no common ground state); link spaces are infinite-dimensional, so stability theorems need a truncation or q-deformation.
- Open: no application of Knabe/martingale/Lemm-Mozgunov to any non-Abelian continuous-group gauge Hamiltonian was found.
- Most promising unexplored angle (conj): **a frustration-free non-Abelian gauge model in the KS family** (SU(2)_k string-net / q-deformed KS, which is commuting-projector at its deconfined fixed point, or a Rokhsar-Kivelson-type deformation of truncated KS whose ground state is the strong-coupling or a PEPS state), prove its gap with a finite-size criterion verified by certified ED, then use stability to move toward KS. The Z2 case (Gorantla-Huang) is the template. Expected endpoint: a gap only in a neighbourhood of the frustration-free point at fixed truncation, i.e. a Result-A candidate on a stepping stone or a Result-C statement of where it stops.

### Line 3. Configuration-space quantum mechanics

- Known (cited): Feynman 1981 heuristic (bounded reduced configuration space ⇒ gap); KKN 1996-98 gauge-invariant variables in 2+1D with a gap m = e²N/2π from the configuration-space measure and a leading-order string tension; Leigh-Minic-Yelnikov glueball spectrum; Gribov region is convex, bounded, and meets every orbit (Dell'Antonio-Zwanziger 1991); orbit space has non-negative O'Neill curvature (Singer 1981, Babelon-Viallet 1981). **Bakry-Émery for YM already exists**: rigorously on the lattice (Euclidean Wilson measure, strong coupling: Shen-Zhu-Zhu 2022; fails for U(N) because of flat U(1) directions: Nissim 2025) and formally in the continuum Hamiltonian picture (Moncrief-Marini-Maitra 2018, Mondal 2023, conditional on a ground state / existence of the theory).
- Open: no proof that the Bakry-Émery tensor of |Ψ₀|² is positive beyond strong coupling. No Hardy-inequality or Payne-Weinberger argument on the Gribov region found.
- Most promising unexplored angle (conj): **Payne-Weinberger-type bound on the Gribov region.** In Coulomb gauge on a finite lattice the FP operator M(A) is affine in A, so log det M(A) is concave on the Gribov region (where M > 0), and the region is convex. The FP-weighted measure is then log-concave on a convex body, for which λ₁ ≥ π²/diam² (Payne-Weinberger for log-concave measures; Bonnefont-Joulin arXiv:2302.05207, not added to the index as it is pure math, ID from search only). Caveats to check before any claim: (a) the Christ-Lee Hamiltonian has an extra non-negative Coulomb term and the potential, so the gap of the full H is not bounded by the gap of the weighted Laplacian without an extra argument (e.g. ground-state transformation, then Bakry-Émery with the potential); (b) the physical domain is the fundamental modular region with boundary identifications, not the Gribov region; (c) diam of the Gribov region grows with lattice volume, so the bound alone vanishes in infinite volume. Concavity of log det M is elementary and checkable now. A finite-volume, fixed-lattice-spacing gap lower bound from it, with the volume-scaling obstruction stated explicitly, is a plausible Result-A-on-stepping-stone or Result-C item.

### Line 4. Band and Bloch structure

- Known (cited): θ is the Bloch quasi-momentum for the Z of winding-number vacua (Jackiw-Rebbi, Callan-Dashen-Gross 1976). Finite-volume Bloch-type effective Hamiltonians with θ as a boundary phase at sphalerons, Rayleigh-Ritz glueball masses (Lüscher 1983, Koller-van Baal 1988, van den Heuvel-van Baal 1995-97). Large N: E(θ) is multi-branched, not a cosine band (Witten 1998). θ = π has a mixed center/T anomaly forbidding a trivially gapped vacuum (Gaiotto-Kapustin-Komargodski-Seiberg 2017), the YM counterpart of LSM.
- Open: nothing obviously open that bears on the gap at θ = 0.
- Judgment (conj): the θ-band describes the vacuum energy E(θ) across sectors; the glueball gap is an excitation inside a fixed-θ sector, so Bloch band-gap reasoning does not produce it. Bandwidth ~ e^{−8π²/g²} → 0 at weak coupling. Most promising use: a clean **Result-C entry** (computed in the van Baal effective model: bandwidth vs 0++ gap as functions of L). Low priority for new work.

### Line 5. Interaction-driven gaps

- Known (cited): dual superconductor ('t Hooft 1981, Mandelstam 1976), proved analytically only with N=2 SUSY softly broken (Seiberg-Witten 1994) and rigorously only for 3D compact U(1) (Polyakov 1977, Göpfert-Mack 1982). Savvidy vacuum and its Nielsen-Olesen tachyon (Landau levels + g=2 Zeeman term of charged gluons). Gap equations for a gluon mass: Cornwall 1982 (pinch technique), Aguilar-Binosi-Papavassiliou (Schwinger mechanism), LCO ⟨A²⟩ effective potential (Verschelde et al. 2001), Gribov-Zwanziger horizon gap equation. **BCS variational vacuum for gluons in Coulomb gauge with glueball TDA/RPA spectrum: Szczepaniak-Swanson-Ji-Cotanch 1996, Szczepaniak-Swanson 2001.**
- Open: none of these has error control. The BCS/mean-field and gap-equation versions are already published, so a new mean-field gap equation is a rederivation.
- Most promising unexplored angle (conj): the precise Hubbard-model reading of KS at strong coupling (on-link Casimir = on-site U, plaquette = correlated hopping) makes the strong-coupling gap a Mott gap; the useful question is a **rigorous variational upper bound plus Temple lower bound on the Coulomb-gauge BCS glueball masses** on a lattice regulator, turning the Szczepaniak-Swanson model into a bounded statement. Moderate priority; likely ends as a Result-B competitor, not Result A.

### Line 6. Stepping stones

- Known (cited): O(N≥3) 2D sigma model mass gap exact via TBA (Hasenfratz-Maggiore-Niedermayer 1990), relying on integrability, **not rigorously proved** at weak coupling. Haldane conjecture: rigorous only at the AKLT point, not for S=1 Heisenberg. 3D compact U(1): rigorous gap at all couplings (Göpfert-Mack). 2D YM: gap elementary, LSI at all couplings (Piazza-Shen 2026). Large-N lattice YM: gap and area law at strong 't Hooft coupling (Chatterjee 2015, Shen-Zhu-Zhu 2022, Cao-Nissim-Sheffield 2025, Nissim 2025). 2+1D and 3+1D SU(N) glueball spectra from Teper, Athenodorou-Teper, Lucini-Teper(-Wenger), Morningstar-Peardon.
- Open: rigorous gap for the 2D O(3) model at weak coupling; Haldane gap for S=1 Heisenberg.
- Most promising unexplored angle (conj): use the **O(3) quantum rotor chain** (Hamiltonian lattice O(3) model, rotor on S² per site, the 1D analog of KS with S² in place of SU(2)) as the test bed for every method above (certified ED + finite-size criterion; Payne-Weinberger/Bakry-Émery on the configuration space). Any method that fails there will fail for YM; one that works there at couplings beyond the strong-coupling regime would itself be notable.

### Glueball models for Result-B comparison (file 08)

Flux tube (Isgur-Paton 1985), bag (Jaffe-Johnson 1976), Coulomb-gauge many-body (Szczepaniak et al.), constituent gluon with helicity (Mathieu-Buisseret-Semay 2008), holographic (Csáki et al. 1999, Brower-Mathur-Tan 2000; one scale, ratios predicted from a 1D Schrödinger problem), 2+1D Leigh-Minic-Yelnikov (Bessel zeros). The holographic and constituent models already are "bound state in an effective potential" models, so a Result-B model of that type must beat them on the same data with honest parameter counts.

## Verdict

Every line of attack has substantial prior work. The literal electron-level analogies (θ as Bloch momentum, BCS vacuum for gluons, Landau-level instability, KS as tight-binding/Mott, coupled cluster for KS glueballs, Bakry-Émery curvature gap, configuration-space bound states) are all published. Novelty is available only in (i) certified/rigorous versions on stepping stones and (ii) the non-Abelian frustration-free gap-proof route (line 2) and Payne-Weinberger on the Gribov region (line 3), neither of which was found in the literature. Recommended priority: line 2 frustration-free non-Abelian model, line 3 Payne-Weinberger check, both tested first on the O(3) rotor chain / SU(2) plaquette chain. Line 4 is expected to give a Result-C entry.

## Open issues

- Absence claims ("no paper found") rest on web search on 2026-11-20; a follow-up search on Google Scholar / INSPIRE citations of Gorantla-Huang 2409.10612, Shen-Zhu-Zhu 2204.12737 and Dell'Antonio-Zwanziger 1991 should be done before any novelty claim.
- Classic-paper summaries are from standard knowledge; reread the source before using a precise statement (especially Koller-van Baal 1988, Balaban, Yarotsky's hypotheses and whether they cover KS with Gauss law).
- Journal refs tagged "unverified" in `literature/` (e.g. Chatterjee 1502.07719 CMP, Shen-Zhu-Zhu CMP, Mondal JHEP) need checking before citation in a result.
- Not covered: Hamer's later SU(3) Hamiltonian series, Kogut-Shigemitsu string-tension series, Durhuus-Fröhlich 1980 (mass gap ⇒ area law), Fröhlich-Spencer, Brydges-Fröhlich-Seiler, Seiler's LNP 159 monograph, Jaffe-Witten Clay problem description, Gribov-Zwanziger refined (Dudal et al.), DMRG for 2+1D non-Abelian LGT beyond the cited TN papers. Add if a task needs them.
- Unvetted claimed proofs (2506.00284, 2611.11534, 2503.15539, Zenodo records) listed in `literature/02-rigorous-lattice.md`; 2611.11534 has explicit hard-cutoff obstruction claims that a skeptic should read if line 1 truncation bounds are pursued.
