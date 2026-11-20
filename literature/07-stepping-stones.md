# 07. Lower-dimensional stepping stones and large N

Line of attack 6. Tags as in `README.md`.

## 2D O(N) sigma models

**Polyakov, "Interaction of Goldstone particles in two dimensions. Applications to ferromagnets and massive Yang-Mills fields", Phys. Lett. B 59, 79 (1975).** DOI 10.1016/0370-2693(75)90161-6 [Crossref]
- Asymptotic freedom of the 2D O(N≥3) sigma model; dynamically generated mass expected.
- Setting: continuum, 1+1D, O(N).

**Zamolodchikov, Zamolodchikov, "Factorized S-matrices in two dimensions as the exact solutions of certain relativistic quantum field theory models", Ann. Phys. 120, 253 (1979).** DOI 10.1016/0003-4916(79)90391-9 [Crossref]
- Exact S-matrix of the O(N) sigma model (and others), assuming integrability and a massive O(N) vector multiplet.

**Hasenfratz, Maggiore, Niedermayer, "The exact mass gap of the O(3) and O(4) non-linear σ-models in d=2", Phys. Lett. B 245, 522 (1990).** DOI 10.1016/0370-2693(90)90685-Y [Crossref]
**Hasenfratz, Niedermayer, "The exact mass gap of the O(N) σ-model for arbitrary N ≥ 3 in d=2", Phys. Lett. B 245, 529 (1990).** DOI 10.1016/0370-2693(90)90686-Z [Crossref]
- Exact m/Λ_MS ratio from TBA with a chemical potential, matched to perturbation theory. Confirmed by lattice MC.
- Setting: continuum, 1+1D, O(N), infinite volume.
- Relevance: the gap is "exact" but rests on the integrability/S-matrix assumption. **No rigorous proof of a gap for the 2D O(3) model at weak coupling exists** (the analog of the Millennium problem is open here too). Any method that proves the O(3) gap at all couplings would be a major result and a real stepping stone.

## Haldane gap

**Haldane, "Continuum dynamics of the 1-D Heisenberg antiferromagnet: Identification with the O(3) nonlinear sigma model", Phys. Lett. A 93, 464 (1983).** DOI 10.1016/0375-9601(83)90631-X [Crossref]
**Haldane, "Nonlinear Field Theory of Large-Spin Heisenberg Antiferromagnets", Phys. Rev. Lett. 50, 1153 (1983).** DOI 10.1103/PhysRevLett.50.1153 [Crossref]
- Spin-S chain ↔ O(3) sigma model with θ = 2πS; integer S gapped, half-integer gapless.
- Relevance: the exact bridge between "electron/spin energy gap" and a relativistic QFT gap. Rigorous gap known for the AKLT point (file 03), not for the S=1 Heisenberg chain. θ = 2πS is the 2D analog of YM θ (file 05).

## 2+1D Yang-Mills (lattice numbers)

**Teper, "SU(N) gauge theories in 2+1 dimensions", Phys. Rev. D 59, 014512 (1999).** arXiv:hep-lat/9804008 [arXiv-API]
**Meyer, Teper, "High Spin Glueballs from the Lattice", Nucl. Phys. B 658, 113 (2003).** arXiv:hep-lat/0212026 [arXiv-API]
**Athenodorou, Teper, "SU(N) gauge theories in 2+1 dimensions: glueball spectra and k-string tensions", JHEP (2017).** arXiv:1609.03873 [arXiv-API; journal unverified]
- Euclidean lattice continuum-extrapolated glueball spectra and string tensions for SU(N), N=2..∞ (extrapolated).
- Setting: Euclidean lattice, 2+1D, SU(N), continuum extrapolated (numerical).
- Relevance: data for KKN / Leigh-Minic-Yelnikov comparisons and for Result-B fits in 2+1D. Goes into `data/`.

## 3+1D and large N (lattice numbers)

**Morningstar, Peardon, "The glueball spectrum from an anisotropic lattice study", Phys. Rev. D 60, 034509 (1999).** arXiv:hep-lat/9901004 [arXiv-API]
**Lucini, Teper, "SU(N) gauge theories in four dimensions: exploring the approach to N = ∞", JHEP 0106:050 (2001).** arXiv:hep-lat/0103027 [arXiv-API]
**Lucini, Teper, Wenger, "Glueballs and k-strings in SU(N) gauge theories: calculations with improved operators", JHEP 0406:012 (2004).** arXiv:hep-lat/0404008 [arXiv-API]
**Athenodorou, Teper, "The glueball spectrum of SU(3) gauge theory in 3+1 dimension", JHEP (2020).** arXiv:2007.06422 [arXiv-API; journal unverified]
- Standard glueball spectra and N-dependence. Belong in `data/` (phenomenologist task).

**'t Hooft, "A planar diagram theory for strong interactions", Nucl. Phys. B 72, 461 (1974).** DOI 10.1016/0550-3213(74)90154-0 [Crossref]
- Large-N expansion; at N=∞ glueballs are free, stable, with an infinite tower.

Rigorous large-N lattice results (Chatterjee 2015; Shen-Zhu-Zhu 2022; Cao-Nissim-Sheffield 2025; Nissim 2025): see file 02. All require strong coupling.

## Exactly solvable gauge stepping stones

- 2D YM: exactly solvable, gap elementary; LSI at all couplings (Piazza-Shen, arXiv:2610.35286, file 02).
- 3D compact U(1): rigorous gap at all couplings (Göpfert-Mack, file 02).
- Qubit-regularized SU(2) plaquette chain → TFIM → E8 massive spectrum (Siew-Chandrasekharan-Bhattacharya, arXiv:2512.11068, file 01).
- Anisotropic 2+1D YM → coupled principal chiral models (Orland, arXiv:hep-th/0607013, file 04).
