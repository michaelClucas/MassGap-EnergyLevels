# 06. Interaction-driven gaps: dual superconductor, Savvidy, condensates, gap equations

Line of attack 5. Tags as in `README.md`.

## Dual superconductor

**'t Hooft, "Topology of the gauge condition and new confinement phases in non-abelian gauge theories", Nucl. Phys. B 190, 455 (1981).** DOI 10.1016/0550-3213(81)90442-9 [Crossref]
- Abelian projection: partial gauge fixing to the Cartan subgroup exposes magnetic monopoles; confinement = monopole condensation (dual Meissner effect).
- Setting: continuum, 3+1D, SU(N), heuristic.

**Mandelstam, "Vortices and quark confinement in non-Abelian gauge theories", Phys. Rep. 23, 245 (1976).** DOI 10.1016/0370-1573(76)90043-0 [Crossref]
- Dual superconductor proposal.

**Polyakov, "Quark confinement and topology of gauge theories", Nucl. Phys. B 120, 429 (1977).** DOI 10.1016/0550-3213(77)90086-4 [Crossref]
- Compact U(1) in 3D (and Georgi-Glashow): monopole plasma gives a photon mass and area law. Made rigorous by Göpfert-Mack (file 02).
- Relevance: the one case where "condensate ⇒ gap" is proved. Debye screening is an electron-gas mechanism.

**Seiberg, Witten, "Monopole Condensation, And Confinement In N=2 Supersymmetric Yang-Mills Theory", Nucl. Phys. B 426, 19 (1994).** arXiv:hep-th/9407087 [arXiv-API]
- Exact low-energy solution of N=2 SU(2); soft breaking to N=1 condenses monopoles, giving confinement with a mass gap.
- Setting: continuum, 3+1D, SU(2), supersymmetric (not pure YM).
- Relevance: dual superconductor demonstrated analytically, but only with SUSY. Shows the mechanism exists; not a pure-YM result.

## Savvidy vacuum and Nielsen-Olesen instability

**Savvidy, "Infrared instability of the vacuum state of gauge theories and asymptotic freedom", Phys. Lett. B 71, 133 (1977).** DOI 10.1016/0370-2693(77)90759-6 [Crossref]
- One-loop effective energy of a constant chromomagnetic field is below the perturbative vacuum: the trivial vacuum is unstable to forming a field condensate.
**Nielsen, Olesen, "An unstable Yang-Mills field mode", Nucl. Phys. B 144, 376 (1978).** DOI 10.1016/0550-3213(78)90377-2 [Crossref]
- The constant-field vacuum has a tachyonic mode (lowest Landau level of the charged gluon with spin aligned): the Savvidy vacuum itself is unstable ("spaghetti vacuum").
- Setting: continuum, 3+1D, SU(2), one loop.
- Relevance: formally similar to Landau levels of electrons in a magnetic field with Zeeman term g=2 driving an instability. **The electron-level analogy (Landau levels + Zeeman) is already the content of these papers.**

## Condensates and dynamical gluon mass

**Shifman, Vainshtein, Zakharov, "QCD and resonance physics. Theoretical foundations", Nucl. Phys. B 147, 385 (1979).** DOI 10.1016/0550-3213(79)90022-1 [Crossref]
- Gluon condensate ⟨(α_s/π)G²⟩ as a power correction in the OPE; numerical estimate from sum rules.

**Cornwall, "Dynamical mass generation in continuum quantum chromodynamics", Phys. Rev. D 26, 1453 (1982).** DOI 10.1103/PhysRevD.26.1453 [Crossref]
- Gauge-invariant (pinch technique) Schwinger-Dyson gap equation for the gluon propagator with a dynamical mass m(q²) decreasing at large q²; m(0) ≈ 500 ± 200 MeV.
- Setting: continuum, 3+1D, SU(3), truncated Schwinger-Dyson.
- Relevance: **"BCS-like gap equation for the gluon mass" already exists.** No error bound on the truncation.

**Aguilar, Binosi, Papavassiliou, "The gluon mass generation mechanism: a concise primer", Front. Phys. (2016).** arXiv:1511.08361 [arXiv-API; journal unverified]
- Schwinger mechanism (massless bound-state poles in the vertex) generates the gluon mass in the PT-BFM Schwinger-Dyson framework; matches lattice Landau-gauge propagators that saturate at q→0.

**Gubarev, Stodolsky, Zakharov, "On the Significance of the Quantity A Squared", Phys. Rev. Lett. 86, 2220 (2001).** arXiv:hep-ph/0010057 [arXiv-API]
- Minimum over gauge orbit of ∫A² is gauge invariant; proposed as a dimension-2 condensate.
**Verschelde, Knecht, Van Acoleyen, Vanderkelen, "The non-perturbative groundstate of QCD and the local composite operator A²", Phys. Lett. B 516, 307 (2001).** arXiv:hep-th/0105018 [arXiv-API]
- Local composite operator (LCO) formalism: renormalizable effective potential for ⟨A²⟩, a nonzero minimum, and a dynamical gluon mass. Structurally the same as the Gross-Neveu / BCS gap equation.
- Relevance: **prior art for a mean-field "Hubbard-Stratonovich + gap equation" treatment of a gluon mass**, in Landau gauge.

## Many-body (BCS) methods in Coulomb gauge

**Szczepaniak, Swanson, Ji, Cotanch, "Glueball Spectroscopy in a Relativistic Many-Body Approach to Hadron Structure", Phys. Rev. Lett. 76, 2011 (1996).** arXiv:hep-ph/9511422 [arXiv-API]
- Coulomb-gauge QCD Hamiltonian; **BCS variational vacuum** gives constituent gluon masses; glueballs in Tamm-Dancoff and RPA; quenched masses agree with lattice.
- Relevance: **direct prior art for "BCS gap for gluons → glueball spectrum".** See also file 08.

**Szczepaniak, Swanson, "Coulomb Gauge QCD, Confinement, and the Constituent Representation", Phys. Rev. D 65, 025012 (2001).** arXiv:hep-ph/0107078 [arXiv-API]
- Variational (Gaussian / BCS-type) vacuum in Coulomb gauge with the FP determinant (Gribov-type) included; self-consistent gap equation for the gluon dispersion; confining Coulomb kernel.
