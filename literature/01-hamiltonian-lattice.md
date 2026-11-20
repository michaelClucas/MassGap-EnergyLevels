# 01. Hamiltonian lattice Yang-Mills as a many-body problem

Line of attack 1 in `OVERSEER_PROMPT.md`. Verification tags: [arXiv-API] = ID and title checked against the arXiv API; [Crossref] = DOI and title checked against Crossref; [unverified] = not checked.

## Foundations

**Kogut, Susskind, "Hamiltonian formulation of Wilson's lattice gauge theories", Phys. Rev. D 11, 395 (1975).** DOI 10.1103/PhysRevD.11.395 [Crossref]
- Defines H = (g²/2) Σ_links E² − (1/g²) Σ_plaquettes Re Tr U_p on a spatial lattice, continuous time, temporal gauge, Gauss law as a constraint. The electric term is the Laplacian on the group manifold per link (a quantum rotor), so the strong-coupling spectrum is Casimir eigenvalues.
- Setting: lattice, any spatial dimension, SU(N) or U(1), finite or infinite volume.
- Relevance: this is the "tight-binding model on the group manifold". Every ED / DMRG / gap-bound task in this program starts from this H.

**Kogut, "An introduction to lattice gauge theory and spin systems", Rev. Mod. Phys. 51, 659 (1979).** DOI 10.1103/RevModPhys.51.659 [Crossref]
- Review that explicitly maps lattice gauge theories onto quantum spin systems and duality; includes the strong-coupling Hamiltonian expansion.
- Setting: lattice, various d, Z2/U(1)/SU(N).
- Relevance: the original spin-system/gauge dictionary. Read before claiming any "spin chain analogy" is new.

**Kogut, Sinclair, Susskind, "A quantitative approach to low-energy quantum chromodynamics", Nucl. Phys. B 114, 199 (1976).** DOI 10.1016/0550-3213(76)90586-1 [Crossref]
- Strong-coupling Hamiltonian perturbation theory for the hadron/glueball spectrum.
- Setting: lattice, 3+1D, SU(3), infinite volume, strong coupling.
- Relevance: first strong-coupling glueball-gap estimates in the Hamiltonian picture.

## Strong-coupling series, t-expansion, coupled cluster, variational

**Hamer, Irving, "SU(2) lattice gauge theory in (2+1)D", Z. Phys. C 27, 307 (1985).** DOI 10.1007/BF01556621 [Crossref]
- Cluster-expansion strong-coupling series for vacuum energy, string tension and the scalar (0+) mass gap.
- Setting: lattice Hamiltonian, 2+1D, SU(2), infinite volume.
- Relevance: the gap at strong coupling as a power series in 1/g⁴. Baseline for any ED gap-vs-coupling curve.

**Horn, Weinstein, "The t expansion: A nonperturbative analytic tool for Hamiltonian systems", Phys. Rev. D 30, 1256 (1984).** DOI 10.1103/PhysRevD.30.1256 [Crossref]
- Lanczos/moment-type expansion of <e^{-tH}>; later extended to SU(2) LGT for the 0++ mass and string tension.
- Setting: lattice Hamiltonian.
- Relevance: a Krylov/Lanczos method, which is the same machinery as ED.

**Schuette, Wichmann, Weichmann, "Application of the Coupled Cluster Method to a Hamiltonian Lattice Field Theory".** arXiv:hep-lat/9805011 [arXiv-API]
**Wichmann, Schuette, Metsch, Wethkamp, "The Coupled Cluster Method in Hamiltonian Lattice Field Theory: SU(2) Glueballs", Phys. Rev. D 65, 094511 (2002).** arXiv:hep-lat/0112015 [arXiv-API]
- Apply the quantum-chemistry coupled-cluster method (exponential ansatz) to the Kogut-Susskind Hamiltonian and compute SU(2) glueball masses. Reported convergence to β≈3 where the 28th-order strong-coupling series breaks down near β≈2.1.
- Setting: lattice Hamiltonian, 2+1D, SU(2), infinite volume.
- Relevance: **prior art for "quantum chemistry methods (coupled cluster) on YM"**. Any CC proposal from this team is a rederivation unless it goes further.

**Carlsson, McKellar, "SU(N) Glueball Masses in 2+1 Dimensions", Phys. Rev. D 68, 074502 (2003).** arXiv:hep-lat/0303016 [arXiv-API]
**Carlsson, McKellar, "SU(N) Lattice Gauge Theory on a Single Cube".** arXiv:hep-lat/0303022 [arXiv-API]
**Carlsson, "Improvement and Analytic Techniques in Hamiltonian Lattice Gauge Theory" (thesis).** arXiv:hep-lat/0309138 [arXiv-API]
- Analytic variational (Rayleigh-Ritz) glueball masses for SU(2)-SU(5) in 2+1D with improved Hamiltonians; single-cube exact spectra.
- Setting: lattice Hamiltonian, 2+1D, SU(N≤5), infinite volume (variational) and one cube (exact).
- Relevance: **prior art for variational bound-state methods applied to KS glueballs**, and for N-dependence.

**Loan, Luo, Luo, "Monte Carlo study of glueball masses in the Hamiltonian limit of SU(3) lattice gauge theory", Int. J. Mod. Phys. A 21, 2905 (2006).** arXiv:hep-lat/0503038 [arXiv-API]
- Anisotropic Euclidean MC extrapolated to a_t → 0 to get Hamiltonian-limit glueball masses.
- Setting: lattice, 3+1D, SU(3).
- Relevance: the reference for "Hamiltonian-limit" SU(3) glueball numbers to compare an ED/TN calculation with.

**Pavel, "Expansion of the Yang-Mills Hamiltonian in spatial derivatives and glueball spectrum", Phys. Lett. B 685, 353 (2010).** arXiv:0912.5465 [arXiv-API]
- Strong-coupling expansion of the SU(2) YM Hamiltonian in the number of spatial derivatives (symmetric gauge), on boxes of size a: the free part is a sum of independent "Yang-Mills quantum mechanics of constant fields" per box, inter-box terms treated perturbatively in λ = g^{-2/3}; spin-0 glueball to order λ².
- Setting: continuum Hamiltonian discretized into boxes, 3+1D, SU(2).
- Relevance: an explicit "atoms on a lattice of points" picture; tight-binding-like expansion of continuum YM.

**Raychowdhury, "Low energy spectrum of SU(2) lattice gauge theory", Eur. Phys. J. C 79, 235 (2019).** DOI 10.1140/epjc/s10052-019-6753-0 [Crossref; content not read]
- Low-energy spectrum of SU(2) KS Hamiltonian in the prepotential / loop basis.
- Setting: lattice Hamiltonian, SU(2), small lattices.
- Relevance: direct prior ED-style spectrum of SU(2) KS.

## Gauge-invariant bases, truncations, quantum simulation

**Byrnes, Yamamoto, "Simulating lattice gauge theories on a quantum computer", Phys. Rev. A 73, 022328 (2006).** arXiv:quant-ph/0510027 [arXiv-API]
- Qubit encodings of U(1), SU(2), SU(3) KS Hamiltonians with truncated link Hilbert spaces.
- Setting: lattice Hamiltonian, general d.
- Relevance: first explicit truncation scheme; the truncation error is the "controlled error" the program needs for a Result-A Hamiltonian truncation.

**Byrnes, Sriganesh, Bursill, Hamer, "DMRG approach to the massive Schwinger model", Phys. Rev. D 66, 013002 (2002).** arXiv:hep-lat/0202014 [arXiv-API]
- DMRG for a 1+1D lattice gauge theory (Hamiltonian, staggered fermions).
- Setting: lattice, 1+1D, U(1).
- Relevance: earliest DMRG gap extraction in an LGT; template for extrapolating the gap.

**Tagliacozzo, Celi, Lewenstein, "Tensor Networks for Lattice Gauge Theories with continuous groups", Phys. Rev. X 4, 041024 (2014).** arXiv:1405.4811 [arXiv-API]
- Gauge-invariant tensor network ansatz (including SU(2)) with truncated representations.
- Setting: lattice Hamiltonian, 1+1D and 2+1D.

**Zohar, Burrello, "A Formulation of Lattice Gauge Theories for Quantum Simulations", Phys. Rev. D 91, 054506 (2015).** arXiv:1409.3085 [arXiv-API]
- Representation (electric) basis formulation of non-Abelian link Hilbert spaces, local Gauss law.

**Dalmonte, Montangero, "Lattice gauge theories simulations in the quantum information era", Contemp. Phys. 57, 388 (2016).** arXiv:1602.03776 [arXiv-API]
**Bañuls et al., "Simulating Lattice Gauge Theories within Quantum Technologies", Eur. Phys. J. D (2020).** arXiv:1911.00003 [arXiv-API; journal unverified]
- Reviews of TN/ED/quantum simulation of LGTs.

**Klco, Stryker, Savage, "SU(2) non-Abelian gauge field theory in one dimension on digital quantum computers", Phys. Rev. D 101, 074512 (2020).** arXiv:1908.06935 [arXiv-API]
- SU(2) plaquette chain, truncated at j_max=1/2, on hardware; spectrum and dynamics.
- Setting: lattice, quasi-1D plaquette chain, SU(2).

**Raychowdhury, Stryker, "Loop, String, and Hadron Dynamics in SU(2) Hamiltonian Lattice Gauge Theories", Phys. Rev. D 101, 114502 (2020).** arXiv:1912.06133 [arXiv-API]
- Loop-string-hadron (LSH) formulation: Gauss law solved locally, gauge-invariant basis with only abelian-like constraints.
- Relevance: the natural basis for ED in this program.

**Davoudi, Raychowdhury, Shaw, "Search for efficient formulations for Hamiltonian simulation of non-Abelian lattice gauge theories", Phys. Rev. D 104, 074505 (2021).** arXiv:2009.11802 [arXiv-API]
- Compares angular-momentum, purely gauge-invariant (loop), and LSH bases for SU(2) cost.

**Ciavarella, Klco, Savage, "A Trailhead for Quantum Simulation of SU(3) Yang-Mills Lattice Gauge Theory in the Local Multiplet Basis", Phys. Rev. D 103, 094501 (2021).** arXiv:2101.10227 [arXiv-API]
- SU(3) plaquette chain in local irreps, truncated; spectrum.

**Atas et al., "SU(2) hadrons on a quantum computer", Nat. Commun. 12, 6499 (2021).** arXiv:2102.08920 [arXiv-API]
- 1+1D SU(2) with matter, gauge fields integrated out, VQE spectrum.

**Zache, González-Cuadra, Zoller, "Quantum and classical spin network algorithms for q-deformed Kogut-Susskind gauge theories", Phys. Rev. Lett. 131, 171902 (2023).** arXiv:2304.02527 [arXiv-API]
- q-deformation of KS (quantum group SU(2)_k) gives a finite local Hilbert space with exact gauge invariance; spin network basis; classical and quantum algorithms.
- Relevance: a regulator with a controlled limit k→∞; also connects KS to Levin-Wen string-nets (frustration-free models, see file 03).

**Hayata, Hidaka, "String-net formulation of Hamiltonian lattice Yang-Mills theories and quantum many-body scars in a nonabelian gauge theory".** arXiv:2305.05950 [arXiv-API; journal unverified]
- Maps truncated SU(2) KS onto a string-net (Levin-Wen-type) model; ED on small lattices.

**Yao, "SU(2) Gauge Theory in 2+1 Dimensions on a Plaquette Chain Obeys the Eigenstate Thermalization Hypothesis", Phys. Rev. D 108, L031504 (2023).** arXiv:2303.14264 [arXiv-API]
**Ebner, Müller, Schäfer, Seidl, "Eigenstate Thermalization in 2+1 dimensional SU(2) Lattice Gauge Theory", Phys. Rev. D 109, 014504 (2024).** arXiv:2308.16202 [arXiv-API]
- Full ED spectra of truncated SU(2) plaquette chains and small 2D lattices.
- Relevance: existing ED codes and spectra at j_max = 1/2, 1; gap numbers available for cross-checks.

**Ciavarella, Burbano, Bauer, "Efficient Truncations of SU(N_c) Lattice Gauge Theory for Quantum Simulation".** arXiv:2503.11888 [arXiv-API]
- 1/N_c-expanded Hamiltonian, Krylov-subspace truncation built from plaquette operators; consistent with Euclidean lattice at small couplings.
- Relevance: a truncation with large-N control; relevant to a large-N Result-A attempt.

**Balaji, Conefrey-Shinozaki, Draper, Elhaderi et al., "Quantum Circuits for SU(3) Lattice Gauge Theory".** arXiv:2503.08866 [arXiv-API]
- Pure SU(3) in 2 and 3 spatial dimensions in a streamlined electric basis.

**Siew, Chandrasekharan, Bhattacharya, "Asymptotic-freedom and massive glueballs in a qubit-regularized SU(2) gauge theory", Phys. Rev. D 113, 114522 (2026).** arXiv:2512.11068 [arXiv-API]
- Qubit-regularized SU(2) plaquette chain maps to the transverse-field Ising model in a field; tuned to a continuum limit with Ising CFT UV and E8 massive "glueballs" in the IR; √σ/m₁ = 0.249(1).
- Setting: lattice Hamiltonian, quasi-1D plaquette chain, SU(2) qubit regularization, continuum limit of the chain model (not of 3+1D YM).
- Relevance: **closest existing "spin chain gap = glueball gap" construction with a continuum limit.** Any plaquette-chain gap study here must cite and compare with it.

## Quantum-information-inspired ansatz

**Milsted, Osborne, "Quantum Yang-Mills theory: an overview of a programme", Phys. Rev. D 98, 014505 (2018).** arXiv:1604.01979 [arXiv-API]
- Temporal-gauge lattice Hamiltonian; quantum parallel transporter; MERA-type ground-state ansatz interpolating strong and weak coupling.
- Relevance: tensor-network ground-state program for continuum YM, a direct predecessor of the "DMRG/TN for gapped ground states" idea.
