# 05. θ-vacua, Bloch structure, small-volume YM

Line of attack 4. Tags as in `README.md`.

## θ-vacua as Bloch states

**Jackiw, Rebbi, "Vacuum Periodicity in a Yang-Mills Quantum Theory", Phys. Rev. Lett. 37, 172 (1976).** DOI 10.1103/PhysRevLett.37.172 [Crossref]
**Callan, Dashen, Gross, "The structure of the gauge theory vacuum", Phys. Lett. B 63, 334 (1976).** DOI 10.1016/0370-2693(76)90277-X [Crossref]
- Large gauge transformations (winding n ∈ π₃(SU(N)) = Z) act as a translation group on the classical vacua |n⟩; physical states are eigenstates of this translation with eigenvalue e^{iθ}: |θ⟩ = Σ e^{inθ}|n⟩. Instantons give the tunneling between neighbours.
- Setting: continuum, 3+1D, SU(N), infinite volume, semiclassical.
- Relevance: **the "θ = Bloch quasi-momentum" analogy is the original framing, not new.** The "lattice" is π₃ = Z of vacua, the "Brillouin zone" is θ ∈ [0, 2π). The band index structure gives E₀(θ) (vacuum energy density, topological susceptibility), not the glueball gap, which is an excitation within a fixed-θ sector.

**Witten, "Theta Dependence In The Large N Limit Of Four-Dimensional Gauge Theories", Phys. Rev. Lett. 81, 2862 (1998).** arXiv:hep-th/9807109 [arXiv-API]
- At large N, E(θ) is multi-branched: E(θ) = min_k f(θ + 2πk), with level crossing at θ = π. Not a single cosine band.
- Relevance: the "band" in θ is not tight-binding at large N (instantons suppressed as e^{-N}); a Bloch-band picture is correct only at small N/semiclassical. Useful for a Result-C statement.

**Gaiotto, Kapustin, Komargodski, Seiberg, "Theta, Time Reversal, and Temperature", JHEP (2017).** arXiv:1703.00501 [arXiv-API; journal unverified]
- Mixed 't Hooft anomaly between Z_N center (1-form) symmetry and time reversal at θ = π: SU(N) YM at θ=π cannot have a trivially gapped unique vacuum (it is gapped with broken T, or deconfined, or gapless).
- Relevance: **the YM analog of Lieb-Schultz-Mattis** (see file 03, Cho-Hsieh-Ryu). Gives an obstruction to a trivial gap at θ=π; says nothing against a gap at θ=0.

## Small and intermediate volume

**Lüscher, "Some analytic results concerning the mass spectrum of Yang-Mills gauge theories on a torus", Nucl. Phys. B 219, 233 (1983).** DOI 10.1016/0550-3213(83)90436-4 [Crossref]
- For small spatial torus L, reduces YM to an effective quantum mechanics of the constant (toron) abelian zero modes; glueball energies ∝ g^{2/3}/L at leading order.
- Setting: continuum Hamiltonian, 3+1D, SU(N), finite volume T³, weak coupling (small L).
- Relevance: **a rigorous-flavoured bound-state calculation of the gap in finite volume**; the anharmonic-oscillator spectrum ("quartic potential well") is a textbook electron-levels problem.

**Lüscher, Münster, "Weak-coupling expansion of the low-lying energy values in the SU(2) gauge theory on a torus", Nucl. Phys. B 232, 445 (1984).** DOI 10.1016/0550-3213(84)90038-5 [Crossref]
- Higher-order effective Hamiltonian for SU(2) on T³; energy levels from Rayleigh-Ritz.

**Koller, van Baal, "A non-perturbative analysis in finite volume gauge theory", Nucl. Phys. B 302, 1 (1988).** DOI 10.1016/0550-3213(88)90665-7 [Crossref]
- Includes the tunneling between gauge-equivalent minima (sphaleron barriers) into the T³ effective Hamiltonian; θ-type phases enter as boundary conditions across the fundamental domain; intermediate-volume glueball spectrum matched to lattice. (Title and DOI verified; content summarized from the later van Baal review arXiv:hep-ph/0008206 and memory, not reread.)
- Relevance: **this is already a Bloch-band calculation of the YM gap in finite volume.** A team "band-structure for YM in small volume" task would rederive this.

**van den Heuvel, van Baal, "Dynamics on the SU(2) fundamental domain", Nucl. Phys. Proc. Suppl. 42, 823 (1995).** arXiv:hep-lat/9411046 [arXiv-API]
**van den Heuvel, "Glueball Spectroscopy on S³", Phys. Lett. B 368, 124 (1996).** arXiv:hep-lat/9509019 [arXiv-API]
**van den Heuvel, "Glueballs on the three-sphere", Nucl. Phys. B 488, 282 (1997).** arXiv:hep-lat/9608101 [arXiv-API]
**van Baal, "Intermediate Volumes and the Role of Instantons".** arXiv:hep-th/9705112 [arXiv-API]
**van Baal, "QCD in a Finite Volume" (review).** arXiv:hep-ph/0008206 [arXiv-API]
- Effective Hamiltonian on the fundamental domain (Gribov region modulo residual copies); θ imposed as a phase relation between wavefunctions at the two sphalerons; Rayleigh-Ritz 0++ and 2++ masses vs coupling and θ.
- Setting: continuum Hamiltonian, 3+1D, SU(2), finite volume S³ or T³.
- Relevance: the most complete existing "particle in a periodic domain with Bloch boundary conditions" computation for YM; lines 3 and 4 meet here. The open piece is the large-volume limit, where the effective-QM reduction fails.
