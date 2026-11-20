# 04. Configuration-space quantum mechanics

Line of attack 3. Tags as in `README.md`.

## Feynman and Karabali-Kim-Nair (2+1D)

**Feynman, "The qualitative behavior of Yang-Mills theory in 2+1 dimensions", Nucl. Phys. B 188, 479 (1981).** DOI 10.1016/0550-3213(81)90005-5 [Crossref]
- Argues heuristically that the gauge-reduced configuration space is compact in the relevant sense (large field configurations are gauge equivalent to small ones), so the Schrödinger operator on it has a gap, in the same way a particle in a bounded domain does.
- Setting: continuum Hamiltonian, 2+1D, SU(2), infinite volume, heuristic.
- Relevance: **the original "YM as a particle in a bounded configuration space" argument**. Line 3 starts here.

**Karabali, Nair, "A gauge-invariant Hamiltonian analysis for non-Abelian gauge theories in (2+1) dimensions", Nucl. Phys. B 464, 135 (1996).** arXiv:hep-th/9510157 [arXiv-API]
**Karabali, Kim, Nair, "Planar Yang-Mills theory: Hamiltonian, regulators and mass gap", Nucl. Phys. B 524, 661 (1998).** arXiv:hep-th/9705087 [arXiv-API]
**Karabali, Kim, Nair, "On the vacuum wavefunction and string tension of Yang-Mills theories in (2+1) dimensions", Phys. Lett. B 434, 103 (1998).** arXiv:hep-th/9804132 [arXiv-API]
- Gauge-invariant matrix variables A = −∂M M⁻¹, H = SL(N,C)/SU(N); the kinetic operator acquires a term m·(number operator) with m = e²N/2π from the measure (WZW Jacobian), giving a mass gap at leading order; vacuum wavefunctional and string tension √σ = e² √((N²−1)/8π) at leading order.
- Setting: continuum Hamiltonian, 2+1D, SU(N), infinite volume, regulated, leading-order approximation.
- Relevance: **the gap comes from the volume measure of the reduced configuration space**, exactly a "kinetic term on curved space gives a gap" mechanism. Already done; compare with Teper's lattice numbers (file 07).

**Leigh, Minic, Yelnikov, "On the Glueball Spectrum of Pure Yang-Mills Theory in 2+1 Dimensions", Phys. Rev. D 76, 065018 (2007).** arXiv:hep-th/0604060 [arXiv-API]
- Extends KKN to a vacuum wavefunctional ansatz; glueball masses from zeros of Bessel functions, compared with lattice.
- Relevance: a k-parameter glueball spectrum model in 2+1D (Result-B competitor).

## Gribov ambiguity and orbit-space geometry

**Gribov, "Quantization of non-Abelian gauge theories", Nucl. Phys. B 139, 1 (1978).** DOI 10.1016/0550-3213(78)90175-X [Crossref]
- Coulomb/Landau gauge does not fix the gauge uniquely (Gribov copies); restriction to the Gribov region (FP operator positive) changes the IR: gluon propagator suppressed, ghost enhanced.

**Singer, "Some remarks on the Gribov ambiguity", Commun. Math. Phys. 60, 7 (1978).** DOI 10.1007/BF01609471 [Crossref]
- No continuous global gauge slice exists for non-Abelian YM on S⁴ (or S³): the bundle A → A/G is nontrivial.

**Singer, "The geometry of the orbit space for non-Abelian gauge theories", Phys. Scr. 24, 817 (1981).** DOI 10.1088/0031-8949/24/5/002 [Crossref]
**Babelon, Viallet, "The Riemannian geometry of the configuration space of gauge theories", Commun. Math. Phys. 81, 515 (1981).** DOI 10.1007/BF01208272 [Crossref]
- The orbit space A/G with the L² metric has (O'Neill) non-negative sectional curvature; Singer suggests (regularized) Ricci positivity could give a spectral gap of the Laplacian on it.
- Setting: continuum, d=3 spatial (and general), compact group.
- Relevance: the geometric input for a Lichnerowicz / Bakry-Émery gap argument.

**Zwanziger, "Local and renormalizable action from the Gribov horizon", Nucl. Phys. B 323, 513 (1989).** DOI 10.1016/0550-3213(89)90122-3 [Crossref]
**Zwanziger, "Renormalizability of the critical limit of lattice gauge theory by BRS invariance", Nucl. Phys. B 399, 477 (1993).** DOI 10.1016/0550-3213(93)90506-K [Crossref]
- Local action implementing the restriction to the Gribov region (horizon function); introduces the Gribov mass parameter γ fixed by a gap equation.
- Relevance: a "gap equation" derived from configuration-space boundedness; links lines 3 and 5.

**Dell'Antonio, Zwanziger, "Every gauge orbit passes inside the Gribov horizon", Commun. Math. Phys. 138, 291 (1991).** DOI 10.1007/BF02099494 [Crossref]
- Proves the Gribov region is convex, bounded in every direction, and meets every orbit.
- Setting: continuum (finite volume / lattice-regulated), compact group.
- Relevance: **the rigorous "bounded domain" statement Feynman wanted**. A Dirichlet-type lower bound on the Gribov region is the obvious next step; the boundary (horizon) is where the FP determinant vanishes, which acts as a Hardy-type weight.

**Vandersickel, Zwanziger, "The Gribov problem and QCD dynamics", Phys. Rep. (2012).** arXiv:1202.1491 [arXiv-API; journal unverified]
- Review of Gribov-Zwanziger.

**Greensite, Olejník, Zwanziger, "Center Vortices and the Gribov Horizon", JHEP 0505:070 (2005).** arXiv:hep-lat/0407032 [arXiv-API]
**Greensite, "Faddeev-Popov spectra at the Gribov horizon", Phys. Rev. D 81, 114011 (2010).** arXiv:1001.0784 [arXiv-API]
**Cucchieri, Mendes, "Ghost sector and geometry in minimal Landau gauge", Phys. Rev. D 88, 114501 (2013).** arXiv:1308.1283 [arXiv-API]
- Lattice studies of the FP spectrum near the horizon; rigorous-style lower bounds on the smallest nonzero FP eigenvalue in terms of the Gribov-region geometry (Cucchieri-Mendes).
- Relevance: numbers and inequalities for the geometry of the reduced space on the lattice.

## Curvature / Bakry-Émery approaches (prior art)

**Moncrief, Marini, Maitra, "Orbit Space Curvature as a Source of Mass in Quantum Gauge Theory".** arXiv:1809.06318 [arXiv-API]
- Revisits Singer: including the potential term, the relevant object is the **Bakry-Émery Ricci tensor** of the measure |Ψ₀|² dμ on the orbit space (Ψ₀ the ground state). Positivity of it gives a gap. Conditional on knowing Ψ₀; uses their Euclidean-signature semiclassical approximation.
- Setting: continuum Hamiltonian, 2+1 and 3+1D, compact group, formal.
- Relevance: **prior art for "Bakry-Émery curvature bound gives the YM gap"**. Not a proof.

**Marini, Maitra, Moncrief, "A Euclidean Signature Semi-Classical Program".** arXiv:1901.02380 [arXiv-API]
- Formal asymptotic construction of ground-state wavefunctionals for scalar and YM fields.

**Mondal, "A Geometric Approach to the Yang-Mills Mass Gap", JHEP 12 (2023) 191 [journal ref from search result, unverified].** arXiv:2301.06996 [arXiv-API]
**Mondal, "Mass and infinite dimensional geometry".** arXiv:2308.09304 [arXiv-API]
- Regularized Bakry-Émery Ricci curvature on the orbit space, if positive, gives a mass gap in 2+1 and 3+1D, assuming a quantized YM theory exists. Author calls it heuristic.
- Relevance: same as above. A team "Bakry-Émery for YM" claim must cite these and Shen-Zhu-Zhu (file 02).

**Shen, Zhu, Zhu, arXiv:2204.12737** (file 02): the rigorous lattice version, using Bakry-Émery on SU(N)^E for the Euclidean Wilson measure at strong coupling.

## Other Hamiltonian / Schrödinger-functional work

**Maitra, "Gauge-invariant ground state for canonically quantized Yang-Mills theory".** arXiv:0804.3015 [arXiv-API]
- Hamilton-Jacobi construction of a zero-energy candidate ground state with a nonlinear factor ordering. Formal.

**Orland, "Integrable Models and Confinement in (2+1)-Dimensional Weakly-Coupled Yang-Mills Theory", Phys. Rev. D 74, 085001 (2006).** arXiv:hep-th/0607013 [arXiv-API]
- Anisotropic 2+1D YM with couplings e, e'; in axial gauge H = H₀ + e'² H₁ with H₀ a set of 1+1D principal chiral sigma models (exactly gapped, exact S-matrix); string tension from form factors at small e'.
- Relevance: a 2+1D gap derived from an exactly gapped 1+1D sigma model. Bridges lines 3 and 6.

## Hardy-inequality search

No paper found using a Hardy-type inequality on the Gribov region or orbit space to bound the YM gap (search 2026-11-20). This appears open.
