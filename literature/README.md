# Literature

Prior work for each line of attack in `OVERSEER_PROMPT.md`. Built 2026-11-20 by the literature scout (task `work/20261120-001-literature-survey/`).

| File | Topic | Line of attack |
|---|---|---|
| `01-hamiltonian-lattice.md` | Kogut-Susskind, strong-coupling series, coupled cluster, variational, ED/TN/quantum simulation, gauge-invariant bases | 1 |
| `02-rigorous-lattice.md` | Osterwalder-Seiler, Balaban, constructive, Chatterjee, Shen-Zhu-Zhu Bakry-Émery gap, Cao-Nissim-Sheffield, unvetted claims | 1, 2, 3, 6 |
| `03-spectral-gap-methods.md` | AKLT, Knabe, martingale, Gosset-Mozgunov, Lemm-Mozgunov, stability, LSM, applications to gauge Hamiltonians | 2 |
| `04-configuration-space.md` | Feynman, Karabali-Kim-Nair, Gribov, Singer, Babelon-Viallet, Zwanziger, Moncrief/Mondal Bakry-Émery | 3 |
| `05-theta-bloch.md` | θ-vacua as Bloch states, GKKS anomaly at θ=π, Lüscher and van Baal finite volume | 4 |
| `06-interaction-driven.md` | Dual superconductor, Savvidy/Nielsen-Olesen, condensates, Cornwall, LCO A², Coulomb-gauge BCS | 5 |
| `07-stepping-stones.md` | O(N) sigma model, Haldane, 2+1D and large-N lattice spectra | 6 |
| `08-glueball-models.md` | Flux tube, bag, constituent gluon, Coulomb-gauge, holographic | Result B comparison |

## Verification tags

- **[arXiv-API]**: arXiv ID, title and authors checked against the arXiv API (export.arxiv.org) on 2026-11-20. Journal refs shown with this tag come from the arXiv metadata.
- **[Crossref]**: DOI, title, journal, volume and page checked against the Crossref API on 2026-11-20.
- **[unverified]** or "journal unverified": not checked; do not rely on it without checking.

Summaries are from abstracts (fetched for most recent papers) and from standard knowledge for classic papers. Where a summary was not checked against the paper it says so. Read the paper before using a specific number or theorem statement in a proof.

## Prior work that already does what a line of attack proposes

- Coupled-cluster (quantum chemistry) on the SU(2) KS Hamiltonian for glueballs: Schuette-Wichmann et al., hep-lat/9805011, hep-lat/0112015.
- Variational / Rayleigh-Ritz KS glueballs for SU(2)-SU(5): Carlsson-McKellar, hep-lat/0303016.
- ED of truncated SU(2) plaquette chains and small 2D lattices: Klco et al. 1908.06935, Yao 2303.14264, Ebner et al. 2308.16202, Hayata-Hidaka 2305.05950.
- Plaquette chain as a spin chain with a continuum limit and glueball-like masses: Siew-Chandrasekharan-Bhattacharya 2512.11068.
- Martingale gap proof for a (deformed, frustration-free) lattice gauge Hamiltonian: Gorantla-Huang 2409.10612 (Z2, 3+1D).
- Bakry-Émery curvature gives the lattice YM gap: Shen-Zhu-Zhu 2204.12737 (rigorous, strong coupling, Euclidean). Bakry-Émery on the continuum orbit space: Moncrief-Marini-Maitra 1809.06318, Mondal 2301.06996 (formal / conditional).
- YM as QM on the reduced configuration space with a gap from the measure: Feynman 1981, Karabali-Kim-Nair hep-th/9705087.
- θ as Bloch momentum: Jackiw-Rebbi 1976, Callan-Dashen-Gross 1976. Finite-volume Bloch-type spectrum: Koller-van Baal 1988, van den Heuvel 1996-97.
- LSM-type obstruction for YM: Gaiotto-Kapustin-Komargodski-Seiberg 1703.00501 (θ=π only).
- BCS vacuum for gluons and glueballs from it: Szczepaniak-Swanson-Ji-Cotanch hep-ph/9511422. Gap equation for a gluon mass: Cornwall 1982, Verschelde et al. hep-th/0105018.
- Glueballs as 1D Schrödinger bound states: holographic models (Csáki et al. hep-th/9806021, Brower-Mathur-Tan hep-th/0003115).
- Landau-level + Zeeman instability of gluons in a background field: Nielsen-Olesen 1978.
- [09-resultC-additions.md](09-resultC-additions.md): RDM/SDP/bootstrap lower bounds, AKLT certified bound, partial and claimed proofs (from the Result-C literature skeptic).
