# 03. Spectral-gap methods from spin systems

Line of attack 2. Tags as in `README.md`.

## Core methods

**Affleck, Kennedy, Lieb, Tasaki, "Rigorous results on valence-bond ground states in antiferromagnets", Phys. Rev. Lett. 59, 799 (1987).** DOI 10.1103/PhysRevLett.59.799 [Crossref]
**AKLT, "Valence bond ground states in isotropic quantum antiferromagnets", Commun. Math. Phys. 115, 477 (1988).** DOI 10.1007/BF01218021 [Crossref]
- Spin-1 chain H = Σ P_2(i,i+1): unique frustration-free VBS ground state (an MPS), exponential decay, and a proven spectral gap in the thermodynamic limit.
- Setting: 1D lattice spin chain, infinite volume.
- Relevance: the model case of a rigorously gapped, frustration-free, SU(2)-symmetric Hamiltonian.

**Knabe, "Energy gaps and elementary excitations for certain VBS-quantum antiferromagnets", J. Stat. Phys. 52, 627 (1988).** DOI 10.1007/BF01019721 [Crossref]
- Finite-size criterion: for frustration-free H = Σ h_i (projectors), if the gap of an n-site open chain exceeds 1/(n−1) (1D), the infinite chain is gapped. Uses H² ≥ Δ₀ H.
- Setting: 1D and some 2D lattices.
- Relevance: converts a finite ED computation into a thermodynamic-limit gap proof. Requires frustration-freeness.

**Nachtergaele, "The spectral gap for some spin chains with discrete symmetry breaking", Commun. Math. Phys. 175, 565 (1996).** DOI 10.1007/BF02099509 [Crossref]
- Martingale method: gap from approximate orthogonality of ground-state projectors of overlapping regions.
- Setting: 1D (generalizable), frustration-free, infinite volume.

**Gosset, Mozgunov, "Local gap threshold for frustration-free spin systems", J. Math. Phys. (2016).** arXiv:1512.00088 [arXiv-API; journal unverified]
- Sharpens Knabe in 1D (threshold ~ 6/(n(n+1))) and 2D square lattices; gapless frustration-free 1D systems have gap O(1/n²).

**Lemm, Mozgunov, "Spectral gaps of frustration-free spin systems with boundary", J. Math. Phys. (2019).** arXiv:1801.08915 [arXiv-API; journal unverified]
- Knabe-type finite-size criteria for 2D/higher-dimensional lattices with open boundaries.

**Lemm, Sandvik, Wang, "Existence of a Spectral Gap in the AKLT Model on the Hexagonal Lattice", PRL 124, 177204 (2020).** arXiv:1910.11810 [arXiv-API]
**Pomata, Wei, "Demonstrating the AKLT spectral gap on 2D degree-3 lattices", PRL 124, 177203 (2020).** arXiv:1911.01410 [arXiv-API]
- Computer-assisted gap proofs: finite-size criterion verified numerically (DMRG/ED) on patches.
- Relevance: **the exact template for "ED on a finite patch + finite-size criterion = rigorous infinite-volume gap"**, but for frustration-free models.

**Kastoryano, Lucia, "Divide and conquer method for proving gaps of frustration free Hamiltonians", J. Stat. Mech. (2018).** arXiv:1705.09491 [arXiv-API]
**Young, "On a bulk gap strategy for quantum lattice models", Rev. Math. Phys. 36, 2460007 (2024).** arXiv:2308.01405 [arXiv-API]
**Rai, Kull, Emonts, Tura, "A Hierarchy of Spectral Gap Certificates for Frustration-Free Spin Systems", Quantum 10, 2065 (2026).** arXiv:2411.03680 [arXiv-API]
- Recursive / SDP-hierarchy generalizations of Knabe and martingale bounds.
- Relevance: the SDP hierarchy is a numerical certificate machine; a candidate tool if a frustration-free gauge model is found.

**Michalakis, Pytel (Zwolak), "Stability of Frustration-Free Hamiltonians", Commun. Math. Phys. 322, 277 (2013).** arXiv:1109.1588 [arXiv-API]
**Bravyi, Hastings, Michalakis, "Topological quantum order: stability under local perturbations", J. Math. Phys. 51, 093512 (2010).** arXiv:1001.0344 [arXiv-API]
- Gap of a frustration-free (or commuting-projector, topologically ordered) Hamiltonian survives small local perturbations, under local topological quantum order.
- Setting: lattice, any d, finite-dimensional local Hilbert spaces.
- Relevance: route from a frustration-free gauge model to a nearby non-frustration-free one (e.g. KS). Needs finite local dimension, so a truncated / q-deformed regulator.

## Lieb-Schultz-Mattis family

**Lieb, Schultz, Mattis, "Two soluble models of an antiferromagnetic chain", Ann. Phys. 16, 407 (1961).** DOI 10.1016/0003-4916(61)90115-4 [Crossref]
**Hastings, "Lieb-Schultz-Mattis in Higher Dimensions", Phys. Rev. B 69, 104431 (2004).** arXiv:cond-mat/0305505 [arXiv-API]
- Half-integer spin per unit cell + U(1) symmetry ⇒ no unique gapped ground state.
- Relevance: an obstruction theorem. Its modern form is an 't Hooft anomaly (next entries, and file 05 for θ=π).

**Cho, Hsieh, Ryu, "Anomaly Manifestation of Lieb-Schultz-Mattis Theorem and Topological Phases", Phys. Rev. B 96, 195105 (2017).** arXiv:1705.03892 [arXiv-API]
- LSM as a mixed anomaly between lattice translation and on-site symmetry.
- Relevance: tells us what the YM analog of LSM is (the mixed center/θ anomaly of Gaiotto-Kapustin-Komargodski-Seiberg, file 05). At θ=0 there is no such anomaly, so no LSM obstruction to a gap.

## Applications to gauge Hamiltonians (prior art check)

**Kitaev toric code / Levin-Wen string nets** (standard; not separately fetched). Z_N and quantum-group gauge theories at their "deconfined" fixed points are commuting-projector, frustration-free, gapped Hamiltonians.

**Gorantla, Huang, "Duality-preserving deformation of 3+1d lattice Z2 gauge theory with exact gapped ground states", Phys. Rev. B 111, 245110 (2025).** arXiv:2409.10612 [arXiv-API]
- A frustration-free point of a deformed 3+1D Z2 gauge theory; **gap proved in the thermodynamic limit with the martingale method**; generalization to Z_N discussed.
- Setting: lattice Hamiltonian, 3+1D, Z2 (Z_N), infinite volume.
- Relevance: **prior art: martingale gap proof already applied to a lattice gauge Hamiltonian** (Abelian, finite group). A non-Abelian continuous-group version was not found.

**Zache et al. arXiv:2304.02527; Hayata-Hidaka arXiv:2305.05950** (see file 01): truncated / q-deformed SU(2) KS written as string-net models. Levin-Wen string-nets for SU(2)_k are frustration-free; the KS electric term is not part of that frustration-free structure.

**Searches found no paper applying Knabe, martingale, or Lemm-Mozgunov criteria to an SU(2)/SU(3) Kogut-Susskind Hamiltonian.** (Web search 2026-11-20; absence of evidence, not proof.)

## Structural facts for this line (cited / standard)

- The KS Hamiltonian is not frustration-free: the electric term's ground state (all links in the trivial irrep) is not annihilated by the plaquette term, and vice versa. Knabe/martingale do not apply directly. [standard, follows from the form of H]
- Local Hilbert space per link is L²(G), infinite-dimensional. Stability theorems (Bravyi-Hastings-Michalakis, Michalakis-Pytel) assume finite local dimension; truncation or q-deformation is needed. [cited assumptions in those papers]
- At strong coupling the gap follows from perturbation theory around the electric ground state (Yarotsky, file 02); a spin-system method adds nothing there unless it reaches further into weak coupling.
