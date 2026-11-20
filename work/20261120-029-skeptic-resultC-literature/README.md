# 20261120-029: Skeptic B (literature and novelty) on the candidate Result C

## Question

Does any published result contradict a "blocked" or "obstruction" row of `work/20261120-025-resultC-draft/claim.md`? Is its literature-check section accurate about what is already known? What is the literature status of the unvetted claimed proofs? Do the citations resolve?

## Method

Read `CLAUDE.md`, `OVERSEER_PROMPT.md` (Result C, literature-check stop condition), `claim.md`, `literature/*.md` and the cited work-folder READMEs. Searched the arXiv API (abstract queries, newest first), INSPIRE-HEP API (records, citing papers) and Crossref (bibliographic lookup), plus general web search for comments on the claimed proofs. Every ID in the tables below was resolved through the arXiv API, INSPIRE or Crossref on 2026-11-20. No compute beyond the HTTP lookups. Full texts were read only where stated; everything else is from abstracts.

## Verdict

**Survives with required changes.** No row's obstruction is contradicted by a published result: no rigorous mass gap for continuous-group lattice YM outside strong coupling or the 't Hooft regime, no Bakry-Emery or LSI for 3D/4D lattice YM at weak coupling, no rigorous 1/N gap for 2D CP^{N-1}, and no gap proof for a non-Abelian KS Hamiltonian beyond strong coupling turned up in 2020-2026. But one row states something the literature contradicts (row 6), one whole lower-bound family is missing from the "examined / not examined" lists, the AKLT numbers are presented without the stronger published bound, the novelty claim for certified KS intervals has unacknowledged close preprints, and the unvetted-proofs list is both misclassified and incomplete.

## Issues, ranked

### 1. (High) Row 6 is contradicted: quantum chemistry does have lower-bound machinery (SDP / RDM / bootstrap)

Row 6 says quantum-chemistry methods "add no lower-bound machinery beyond certified ED + Schur (row 1)". That is false as written. The variational 2-RDM method relaxes N-representability to semidefinite constraints and gives rigorous **lower** bounds on E_0 (Nakata et al., J. Chem. Phys. 114, 8282 (2001), doi:10.1063/1.1360199; Mazziotti, PRL 93, 213001 (2004), doi:10.1103/PhysRevLett.93.213001). The same positivity idea now gives:
- rigorous two-sided bounds on Wilson loops in infinite-volume lattice YM at all couplings, d = 2, 3, 4 (Kazakov-Zheng arXiv:2203.11360, PRD 107, L051501 (2023); SU(2) arXiv:2404.16925, JHEP 03 (2025) 099; SU(3) Guo-Li-Yang-Zhu arXiv:2502.14421, JHEP 12 (2025) 033). These bound observables, not the gap;
- rigorous thermodynamic-limit gap lower bounds for a non-frustration-free Hamiltonian, the 1+1D TFIM (Nancarrow-Xin arXiv:2211.03819, JHEP 2023, 052), and a convergent SDP gap-certificate hierarchy that needs neither frustration-freeness nor geometric locality (Rao arXiv:2510.08427, finite-dimensional qubits).

None of these bounds a YM gap, so no obstruction is refuted. But the claim's structure is "every examined family is blocked" plus an open list. This family is neither examined nor listed as open. It also weakens row 2(a) as a blanket statement: non-frustration-freeness blocks Knabe/GM/martingale, not every spin-system gap-certification method.

**Fix:** Row 6, "what blocks it": replace "These methods add no lower-bound machinery beyond certified ED + Schur (row 1)" with "Wavefunction methods (HF, CI, CC, DMRG) add no lower-bound machinery beyond certified ED + Schur (row 1). RDM/SDP relaxations (v2RDM; Nakata et al. 2001, Mazziotti 2004) give lower bounds on E_0, and SDP/bootstrap hierarchies give gap lower bounds without frustration-freeness (Nancarrow-Xin 2211.03819; Rao 2510.08427). These were not examined here." Add one bullet to "What the program did NOT examine": "SDP/bootstrap lower-bound hierarchies (v2RDM, lattice-YM positivity bootstrap 2203.11360 / 2404.16925 / 2502.14421, gap bootstrap 2211.03819 / 2510.08427 / 2610.18854). They are the one published route to rigorous lower bounds that needs neither frustration-freeness nor a truncation floor; no application to a gauge-theory gap was found." Add the same line to "What remains open".

### 2. (Medium-high) Novelty of the certified KS intervals: close preprints are not acknowledged

"What is ours" says the certified two-sided KS gap intervals have no prior ("no prior certified intervals found, `001`"). Two unrefereed preprints that the program already lists as "unvetted full proofs" overlap with rows 1 and 2:
- Mihajlovic arXiv:2611.11534 (8 Nov 2026). The abstract claims "explicit bounds ... on the entire physical vacuum complement", "finite-matrix criteria include all omitted representations", "quantitative fidelity obstructions for hard electric cutoffs, and necessary growth of retained representation spaces", "Schur estimates" isolating what a volume-uniform gap needs, and a certification "through Wilson coupling 0.354". The table of contents has a section titled "Certified infinite-representation one-plaquette benchmark" with "gap intervals". That is row 1's method (complement floor + Schur), row 1(a)'s obstruction (retained representation space must grow), and a certified plaquette interval. It is in transfer-matrix (Euclidean-time) form, not the same H_KS normalization. The body was not readable through the fetch tool, so the overlap in numbers is unconfirmed.
- Gauvin arXiv:2503.15539 (v3) claims a volume-independent gap >= 8a/3 for the SU(3) Hamiltonian at fixed spatial cutoff when b/a <= 1/648. That is a strong-coupling Hamiltonian gap certificate of the C6/Yarotsky type.

Neither is refereed and the claim must not rely on them. But a literature check that calls our intervals unprecedented has to name them.

**Fix:** In "What is ours", replace "(no prior certified intervals found, `001`)" with "(no refereed prior certified intervals found; the unrefereed preprints arXiv:2611.11534 (complement and Schur bounds, a certified one-plaquette benchmark, and a hard-cutoff obstruction requiring growth of retained representations) and arXiv:2503.15539 (fixed-cutoff SU(3) gap certificate at strong coupling) overlap with rows 1-2 and were not checked against our numbers)." Row 1(a): add "cf. the hard-electric-cutoff obstruction claimed in arXiv:2611.11534".

### 3. (Medium) The unvetted-proofs list is misclassified and incomplete

- arXiv:2611.11534 is not a claimed full proof. Its abstract ends: "A nontrivial four-dimensional spatial-cutoff construction and a positive spectral bound along it are not established; no exact continuum mass is asserted." Gauvin's 2503.15539 also says the continuum gap "remain[s] unproved". `literature/02` describes them correctly. STATE.md and claim.md call them claimed full proofs.
- Missing: **Faizal-Shabir arXiv:2606.19362, Fortschr. Phys. 74 (2026) e70097** (INSPIRE-listed, peer-reviewed journal, 0 citations). It claims a continuum 4D SU(N) YM with a mass gap and confinement via reflection positivity and multiscale analysis. It is part of a four-paper series that a web search attributes to IJGMMP 2026; that attribution is unverified. Because it is in a journal, the literature check must name it and state that it is not relied on.
- Status of Jacobsen arXiv:2506.00284 (physics.gen-ph): INSPIRE citation count 0. No comment, erratum or refutation was found (arXiv, INSPIRE, web). A Zenodo mirror exists (record 15786684), and web search attributes a NASA affiliation to it, unverified. Not accepted, not refuted.
- "The Zenodo record": `literature/02` names no specific Zenodo record, only "several Zenodo/ResearchGate/Academia records". Web search shows many such records (e.g. Zenodo 18906249 claims an SU(2) gap for beta >= 9587; 20547755 is a Lean U(1) strong-coupling gap). None has independent review. A Cambridge Open Engage SU(N) mass-gap preprint by another author is marked RETRACTED.

**Fix:** claim.md "NOT examined", last bullet: "Unvetted claimed full proofs: arXiv:2506.00284 (gen-ph, 0 INSPIRE citations, no published comment found) and arXiv:2606.19362 (Fortschr. Phys. 74 (2026) e70097, 0 citations, no published comment found); various Zenodo records (none independently reviewed). Explicitly partial and unrefereed: arXiv:2611.11534, 2503.15539 (see row 1). None is relied on or refuted." The overseer should make the same change to STATE.md's open list.

### 4. (Medium) AKLT numbers: the published bound is stronger, and the GM curve for AKLT is published

Rows 2 and 7 quote GM >= 0.3011 (n = 10, dense) and 0.3026 (n = 12, Lanczos) for the AKLT ring with no comparison. Rai-Kull-Emonts-Tura-Schuch (arXiv:2411.03680, Quantum 10, 2065 (2026), already in `literature/03`) prove Delta_AKLT >= 0.34976 in the thermodynamic limit (Sec. 4.1; the SDP output is turned into a certified solution, App. D). They also quote Knabe 0.24806 and plot the GM bound for AKLT (Fig. 2). Rozman-Molnar-Schuch arXiv:2607.19078 (2026) give improved martingale bounds. So the GM evaluation is a recomputation of a published curve, and the best published certificate is about 16% above ours. This does not touch the obstruction, since AKLT is calibration only. But the reader should not take 0.30 as the state of the art. It also removes the need to make the Lanczos-dependent n = 12 value carry weight.

**Fix:** Row 2, after "[N23]", add: "(recomputation of the Gosset-Mozgunov curve, which is plotted for AKLT in arXiv:2411.03680, Fig. 2. Best published certified bound: Delta_AKLT >= 0.34976, Rai et al., Quantum 10, 2065 (2026); Knabe gives 0.24806 there)". Under Known weaknesses, add that the Lanczos-dependent eps_12 is superseded by a published certified bound.

### 5. (Low-medium) RK attribution in the literature check

The claim credits the RK construction to a "Fendley-Fradkin remark". The general construction (stochastic-matrix-form RK Hamiltonian, ground state sqrt of the Boltzmann weight, spectrum equal to the classical relaxation spectrum) is Castelnovo-Chamon-Mudry-Pujol, Ann. Phys. 318, 316 (2005), doi:10.1016/j.aop.2005.01.006, and Henley cond-mat/0311345. Work 007 cites both, but claim.md does not. The identity gap(H_FF) = (1/4) x Poincare constant of the Wilson measure is the standard ground-state transform of the Langevin generator, so C6 = SZZ in Hamiltonian form (already conceded). 2D U(1) and SU(2) RK dynamics being diffusive is consistent with Chevyrev-Qu-Shen arXiv:2608.27828: the weak-coupling 3D U(1) Langevin dynamics converges to a stochastic heat equation (z = 2). That supports row 2(d).

**Fix:** Literature check: "the RK construction (Henley cond-mat/0311345; Castelnovo-Chamon-Mudry-Pujol, Ann. Phys. 318, 316 (2005); Fendley-Fradkin remark cond-mat/0502071)". Row 2(d): add "consistent with the z = 2 heat-equation limit of 3D U(1) Langevin dynamics, arXiv:2608.27828".

### 6. (Low) Smaller citation and wording fixes

- Row 7 / literature check: the O(3) rotor strong-coupling series has a published longer version, Hamer-Kogut-Susskind PRD 19, 3091 (1979), doi:10.1103/PhysRevD.19.3091 (015 marks it "[not checked]"). Cite it as prior art for [N49] and say whether the five coefficients agree.
- Row 5: "Froehlich-Spencer prove m = 0 there at large beta". FS's threshold is not explicit, so "there" (beta_XY = 2) is not literally covered by the theorem. Write "FS prove m = 0 for all sufficiently large beta (non-explicit threshold); beta_XY = 2 lies above the numerically known KT point".
- Row 5 / `literature/02`: Goepfert-Mack is no longer the only general all-coupling 3D confinement statement. Chatterjee arXiv:2602.00436 (2026) gives a self-contained proof of log-confinement for 3D lattice gauge theories whose group contains the central U(1) (Froehlich comparison + Glimm-Jaffe). It is not a gap result, and it supports row 5's "control only for Abelian / central-U(1) mechanisms".
- Known weaknesses: of the "from memory" citations, the following now resolve with matching bibliographic data (Crossref): Bazley-Fox PR 124, 483 (1961); Hainzl-Seiringer PRB 77, 184517 (2008); FHNS J. Geom. Anal. 17, 559 (2007); Savvidy PLB 71, 133 (1977); Nielsen-Olesen NPB 144, 376 (1978); Hasenfratz-Maggiore-Niedermayer PLB 245, 522 (1990); Goepfert-Mack CMP 82, 545 (1982); Aizenman-Holley IMA Vol. Math. Appl. 8, 1 (1987); Wu Ann. Probab. 34, 1960 (2006). Maehly 1952 (Helv. Phys. Acta) still does not resolve in Crossref. Theorem *content* for these was not re-read here; the bibliographic part of the weakness can be narrowed.
- `literature/02` journal refs marked unverified now resolve: Chatterjee CMP 366, 203 (2019), doi:10.1007/s00220-019-03353-3; Shen-Zhu-Zhu CMP 400, 805 (2023), doi:10.1007/s00220-022-04609-1; Gosset-Mozgunov J. Math. Phys. 57 (2016), doi:10.1063/1.4962337.
- Row 8: the FBSE ID 2503.03821 is published as Eur. Phys. J. C 85 (2025), doi:10.1140/epjc/s10052-025-14590-3. Add the journal ref.

## Searches with no contradicting hit (supports the rows)

| target | result |
|---|---|
| rigorous lattice YM gap beyond strong coupling, 2020-2026 | Only SZZ 2204.12737 (strong), Nissim 2510.22788 (U(N), 't Hooft regime), Cao-Nissim-Sheffield 2505.16585 / 2509.04688 (area law, not gap), Adhikari-Cao 2202.10375 (finite groups). Faria da Veiga-O'Carroll 2509.03513 (charge conjugation, Elitzur; no gap). Supports rows 1, 2. |
| Langevin / stochastic quantization gap, 3D/4D lattice YM | None. Chevyrev-Qu-Shen 2608.27828 is a 3D U(1) scaling limit (finite volume, no gap). Supports row 2(d). |
| Bakry-Emery / LSI for 3D/4D lattice YM at weak coupling | None. Only Piazza-Shen 2610.35286 (2D, torus, all beta; confirmed periodic). Supports row 2(c). |
| rigorous 1/N gap, 2D CP^{N-1} | None found (only formal 1/N work). Open item stands. |
| rigorous theta/Bloch statements for YM on T^3 | None new. Row 4 stands. |
| Gribov-region spectral bounds, Payne-Weinberger | Tedesco 2607.13228 (Gribov-horizon spectral flow of the FP operator, hedgehog threshold states; not an orbit-space Laplacian gap). Cucchieri-Mendes 1308.1283 bound FP eigenvalues. No PW-on-Gribov paper found. Row 3 stands; cite 2607.13228 as adjacent. |
| frustration-free / spectral-gap proofs for non-Abelian gauge Hamiltonians, 2024-2026 | None beyond Gorantla-Huang 2409.10612 (Z2) and the unrefereed preprints in issue 2. |
| certified / interval glueball or KS bounds | No refereed paper. Unrefereed 2611.11534 / 2503.15539 (issue 2). |
| functional glueball spectra after 2025 | Huber-Fischer-Sanchis-Alepuz EPJC 85 (2025) [=2503.03821], Carmo Terin 2601.19727, Chevalier-Mathieu review 2609.16790: all already in 022. Row 8 stands. |

## New references found (all IDs verified)

| reference | ID | bearing |
|---|---|---|
| Nakata et al., J. Chem. Phys. 114, 8282 (2001) | doi:10.1063/1.1360199 | contradicts row 6 wording (RDM SDP lower bounds) |
| Mazziotti, PRL 93, 213001 (2004) | doi:10.1103/PhysRevLett.93.213001 | same |
| Kazakov-Zheng, Bootstrap for lattice YM, PRD 107, L051501 (2023) | arXiv:2203.11360 | missing technique family; bounds observables, not gap |
| Kazakov-Zheng, finite-N SU(2), JHEP 03 (2025) 099 | arXiv:2404.16925 | same |
| Guo-Li-Yang-Zhu, SU(3) bootstrap, JHEP 12 (2025) 033 | arXiv:2502.14421 | same |
| Nancarrow-Xin, gap bootstrap (TFIM), JHEP 2023, 052 | arXiv:2211.03819 | gap lower bound without FF; weakens row 2(a) as a blanket statement |
| Rao, SDP gap certificates for qubit Hamiltonians (2025) | arXiv:2510.08427 | same, finite local dimension |
| Scheer-Khalaf, exact lower bounds on Laughlin gap (2026) | arXiv:2610.18854 | electron-level (FQH) SDP gap certificate in the thermodynamic limit (FF model) |
| Rai-Kull-Emonts-Tura-Schuch, Quantum 10, 2065 (2026) | arXiv:2411.03680 | AKLT >= 0.34976 certified; supersedes rows 2/7 numbers |
| Rozman-Molnar-Schuch (2026) | arXiv:2607.19078 | improved martingale bounds; row 2 |
| Mihajlovic (2026, unrefereed) | arXiv:2611.11534 | overlaps row 1 method, plaquette benchmark, and row 1(a) obstruction; partial, not a full-proof claim |
| Gauvin v3 (unrefereed) | arXiv:2503.15539 | fixed-cutoff SU(3) gap certificate at strong coupling; row 2; partial |
| Faizal-Shabir, Fortschr. Phys. 74 (2026) e70097 | arXiv:2606.19362 | claimed full proof in a journal; must be listed as unvetted |
| Castelnovo-Chamon-Mudry-Pujol, Ann. Phys. 318, 316 (2005) | doi:10.1016/j.aop.2005.01.006 | prior art for RK construction (row 2) |
| Chevyrev-Qu-Shen (2026) | arXiv:2608.27828 | supports row 2(d) (z = 2 limit) |
| Chatterjee, central-U(1) confinement in 3D (2026) | arXiv:2602.00436 | supports row 5; updates `literature/02` |
| Tedesco, Gribov-horizon spectral geometry (2026) | arXiv:2607.13228 | adjacent to row 3; no contradiction |
| Hamer-Kogut-Susskind, PRD 19, 3091 (1979) | doi:10.1103/PhysRevD.19.3091 | prior O(3) strong-coupling series (row 7) |

## Citation sample (claim.md and work folders)

Resolved and matching (arXiv API): 2610.35286, 2204.12737, 2401.13299, 2510.22788, 2505.16585, 2509.04688, 1512.00088, 2409.10612, 2411.03680, 2004.00415, 2110.09180, 2503.03821, 1809.06318, 2301.06996, 2308.09304, 0705.2024, 1301.1162, cond-mat/0311345, cond-mat/0502071, 2506.00284, 2611.11534, 2503.15539. Crossref: the 9 "from memory" items in issue 6 plus Kupiainen CMP 73, 273; KMR CMP 169, 121; Froehlich-Spencer CMP 81, 527; Szczepaniak et al. PRL 76, 2011; Durhuus-Froehlich CMP 75, 103; Hamer-Irving Z. Phys. C 27, 307; Dell'Antonio-Zwanziger CMP 138, 291; Koller-van Baal NPB 302, 1. Unresolved: Maehly 1952. Simon 1983 part I (Ann. IHP 38, 295) is not in Crossref; only part II resolves. Content claims were checked against abstracts only, except Piazza-Shen (torus, periodic: confirmed) and Rai et al. (AKLT numbers: confirmed from the HTML full text).

## Open issues

- The 2611.11534 body (sections 41, 45) could not be read through the fetch tool. Someone should compare its one-plaquette intervals and hard-cutoff theorem with [N6]/[N7] and row 1(a).
- Whether the C8c U(1) counterexample and the sharpness of the -8(d-1) Hessian constant appear in SZZ or Piazza-Shen remarks was not checked (full texts not read).
