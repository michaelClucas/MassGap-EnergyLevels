# 20261120-018: Skeptic audit of C10 (work 014) and C11 (work 015)

## Question

Do C10 (interaction-driven gaps, line 5) and C11 (toolkit calibration on rotor/AKLT/Heisenberg chains, line 6) survive an attempt to break them? Each needs real citations that say what is claimed, calculations that rerun and rederive, and a claim that does not overstate what is blocked.

## Method

- Reran `work/20261120-014-line5-interaction-gaps/interaction_gaps.py` (pinned, 1 thread). Output `rerun_014.log` is identical to the original `run.log`.
- Wrote `skeptic_checks.py`, which imports nothing from 014 or 015. It has its own rotor matrix elements, AKLT chain, PT code and DVR. Output is in `run.log`. It ran in 449 s on 1 thread and stayed under 1 GB. Command:
  `python tools/pin.py .venv/Scripts/python.exe -u work/20261120-018-skeptic-c10-c11/skeptic_checks.py`
- Checked citations on the web (arXiv, publisher and index pages, 2026-11-20). For Gosset-Mozgunov and Garcia-Saez-Murg-Wei I read the full text. For the others only the abstract or secondary sources were available, as noted per item.

## Calculation checks

| item | claimed (014/015) | independent result | verdict |
|---|---|---|---|
| Peierls gain, half-filled spinless ring, t(1 +- d) | a = 0.31812 fit, ~1/pi | Closed form E/N = -(2/pi) E(1-d^2) (complete elliptic integral of the 2nd kind) gives gain = (1/pi) d^2 [ln(1/d) + ln 4 - 1/2] + o(d^2). It matches quadrature to 1e-12, and ratio -> 1.000000 at d = 1e-3. b = (ln4 - 1/2)/pi = 0.2821 (fit gave 0.283). | survives |
| one-loop 4-gluon tadpole | Pi ~ (d-1+xi)(d-1)/d | Rederived by hand. Contracting the 4-gluon vertex with delta^{cd} D_rs kills the f^{abe}f^{cde} term. The other two terms give 2 C_A delta^{ab}(g_mn D_rr - D_mn), and D_rr = (d-1+xi)/k^2. Symmetric integration then gives (d-1+xi)(d-1)/d. | algebra survives; interpretation weakened (below) |
| Gribov Schwinger function | C(t) = (gamma/2sqrt2) e^{-a}(cos a - sin a) | The prefactor is **1/(2 sqrt2 gamma)**: checked by quadrature at gamma = 1, 2, and on dimensional grounds. 014's script uses gamma = 1, so its table is right and the README formula has a typo. Sign change at a = pi/4 confirmed. | survives (typo) |
| anharmonic QM, Gaussian gap / exact gap | 1.01x to 7e5x | sinc-DVR, converged to 1e-12: exact gaps 1.210356, 1.934122, 1.725658, 0.2968899, 1.508778e-2, 1.039557e-4, **2.077103e-7** (014 had 2.0790e-7). Ratios Omega_sym/exact are 1.009 to **7.18e5**. | survives. Caveat: the 7e5 uses the symmetric Hartree solution, which for mu <= -2 is not the variational minimum. The minimizing (broken) Gaussian predicts 0, so it errs in the other direction. |
| 2D spherical model | m ~ sqrt(32) e^{-2 pi beta} | W_2(m) vs (1/4pi) ln(32/m^2): agree to 1e-5 at m = 1e-2 and 1e-8 at m = 1e-4. | survives |
| AKLT eps_n (h = P_2, a projector) | eps_12 = 0.36436, GM >= 0.3026, Knabe >= 0.3008 | P_2 = S.S/2 + (S.S)^2/6 + 1/3 checked to be a projector. Dense, no Lanczos: eps_6, 8, 10 = 0.39845123, 0.37934913, 0.36984447, giving **GM >= 0.30107 at n = 10 with no Lanczos assumption**. Lanczos: eps_12 = 0.36435607, giving GM 0.30262 and Knabe 0.30079. Ring n = 12 gap 0.350120. | survives |
| GM / Knabe constants | as in 015 | Read from arXiv:1512.00088 v1, p. 5. Thm 1 (Knabe): n > 2, m > n, eps_m(ring) >= ((n-1)/(n-2))(eps_n - 1/(n-1)). Thm 3: n > 2, m > 2n, eps_m(ring) >= (5/6)((n^2+n)/(n^2-4))(eps_n - 6/(n(n+1))). The h_{i,i+1} are projectors, which matches 015. | survives |
| rotor strong-coupling series | gap/g = 1 - (2/3)x + (2/27)x^2 + (7/405)x^3 + 0.0026194x^4 | Own RSPT code (ring L = 5, M = 0, k = 0 magnon). l_max = 2 gives 1, -0.66666667, 0.07407407, 0.01728395. l_max = 3 gives x^4 = 0.00261937. E0: -(L/6)x^2. Analytic checks: the magnon hopping amplitude <0|n^a|1b> squared is 1/3, so the x^1 term is -2x/3. The vacuum second order is -(1/2)\|n.n\|0,0>\|^2 = -1/6 per bond. | survives |
| rotor certified intervals | e.g. L=2 open g=0.25 [0.1318,0.1318]; L=3 ring g=0.35 [0.12446,0.12824] | Own ED at larger l_max (12 for L = 2, 7 for L = 3): 0.131797, 0.288016 (L=2, g=0.5, interval [0.28799, 0.28802]), 0.131623 (L=3 open, [0.13109, 0.13204]), 0.126603 (L=3 ring). Each value lies inside its certified interval. | survives |
| Bakry-Emery Hessian | Hess(-n1.n2)(u,u) = (n1.n2)(\|u1\|^2+\|u2\|^2) - 2u1.u2 >= -2(...), per site -4, sharp | Formula rederived from the Riemannian Hessian on S^2. On an N = 6 ring the Neel configuration gives min eigenvalue exactly -4 per \|u\|^2 (sharp). Random configurations give >= -3.46. The ground-state-transform identity H_FF = sum L^2 - 2beta sum n.n + (beta^2/4) sum \|grad_i S\|^2 was rederived. The criterion form is CD(rho, inf): rho = Ric + Hess(beta S) >= 1 - 4beta. It is not sharp for the gap itself: at beta = 0 it gives 1, while the S^2 gap is 2. | survives |

## Citation checks

| citation | what is used | finding |
|---|---|---|
| Froehlich-Spencer, CMP 81, 527 (1981) | m = 0 for 2D XY at large beta | Abstract (IHES preprint, journal record): they prove a KT transition for the rotator, Villain, SOS and Z_n models, with power-law bounds at large beta. Used correctly. |
| Kopper-Magnen-Rivasseau, CMP 169, 121 (1995) | rigorous mass generation in large-N GN_2 | Index record plus secondary descriptions: they prove exponential decay of the 2-point function for large N. Full text not read. Used correctly. |
| **Goepfert-Mack, CMP 82, 545 (1981/82)** | C10: "3D compact U(1) **with Higgs**" | **Misattributed.** The paper treats pure 3D U(1) lattice gauge theory with the Villain action, with no matter field. It proves confinement for all couplings (title; Hamburg record; Chatterjee arXiv:2602.00436 abstract). The Higgs belongs to Polyakov's Georgi-Glashow argument, which is not rigorous. The obstruction itself still holds: the continuum limit of 3D compact U(1) is a free photon, and the control is a dilute monopole gas with fugacity ~ e^{-c beta}, which non-Abelian YM lacks. |
| Kupiainen, CMP 73, 273 (1980) | rigorous 1/N, O(N) | Springer record plus later reviews: a mass gap for large N above the spherical critical point, which in d <= 2 means all beta. As read from summaries, N >= N0(beta), not uniform in beta. Consistent with C10. |
| Kennedy-Lieb, PRL 59, 1309 (1987) | Peierls, classical static phonons | Abstract: they prove the minimizer is period 2 for the static SSH-type (Hueckel) model. Used correctly. |
| Hasenfratz-Maggiore-Niedermayer | m = (8/e) Lambda_MSbar, O(3) | The general-N abstract, PLB 245, 529 (1990), gives m = (8/e)^{1/(N-2)} / Gamma(1+1/(N-2)) Lambda_MSbar, which is 8/e at N = 3. The O(3)/O(4) companion is PLB 245, 522; I did not open it. Shin (hep-lat/9611006) states the derivation rests on the TBA, without a rigorous proof. C11's "integrability, not rigorous" is right. |
| Todo-Kato, PRL 87, 047203 (2001) | Haldane gap 0.41048(6) | Confirmed (arXiv cond-mat/9911047 abstract). |
| Garcia-Saez-Murg-Wei, arXiv:1308.3631 | AKLT gap ~0.350 | Confirmed: "Delta ~ 0.350", with H a sum of projectors onto pair spin 2, the same normalization as 015. |
| Nielsen-Olesen NPB 144, 376; Savvidy PLB 71, 133; Hainzl-Seiringer PRB 77, 184517; FHNS J. Geom. Anal. 17, 559 | as in 014 | Not re-checked by web here; 014 checked the bibliographic data with Crossref. The physics used (N-O tachyon E^2 = k_z^2 - B; HS = BCS functional gap and Tc asymptotics, not many-body) is standard and matches section G2, which was recomputed. These stay "cited, content from memory". |

## Scope of "no mechanism gives a gap lower bound for continuum YM"

- I found no rigorous large-N or 1/N result for a gauge theory that contradicts it:
  - 2D CP^{N-1}: searching turned up only formal 1/N work (Witten 1979; Campostrini-Rossi), no rigorous mass-gap proof. This matches 014's "conjectured open" tag, though the search was not exhaustive.
  - 't Hooft model: the bound-state equation is derived in planar perturbation theory, not rigorously from QCD_2.
  - Rigorous large-N lattice gauge results (Chatterjee 2019; Basu-Ganguly) are strong-coupling only.
  - Kupiainen and KMR are vector models.
- **Overstatement 1:** C10's setting lists d = 2. Continuum 2D YM has no propagating gluons and is exactly solvable (Migdal, Driver, Levy), so "mass gap" there is either trivial or ill-posed. The claim should be restricted to d = 3, 4.
- **Overstatement 2:** "no mechanism" should say "none of the mechanisms examined". What was examined is: BCS, Mott, Peierls, excitonic, dual superconductor, Gaussian/Hartree/variational, gap equations and vector large N.
- **Straw-man risk, section E.** The (d-1+xi) dependence holds only for the cutoff Hartree tadpole. No current program claims that object is a physical mass. Pinch-technique and Coulomb-gauge-variational masses are xi-independent or use a physical gauge, and E says nothing about them. They are covered by D (a variational calculation bounds E0 only) and F (positivity violation). E should be presented as an illustration that the Hartree mass is unphysical, not as a counterexample.
- **Counterexample A** applies the N = inf equation at N = 2. Nobody claims that equation is valid at N = 2. It supports only the narrow point "a massive gap-equation solution does not imply a gap", which is how the claim words it. It is acceptable as stated.

## Lanczos completeness in C11 (and C1, C5, C9)

This is a real gap in rigor, but a narrow and closable one.
- The Schur-complement lower bound uses lambda_k of a computed matrix A - X/(d - U_k), and also U_k < d. If eigsh skips an eigenvalue and returns lambda_{k+1} in place of lambda_k, the stated lower bound is too high. Residual bounds do not help: they show that *some* eigenvalue lies near each Ritz value, not that it is the k-th one.
- Upper bounds (min-max on the P space) are unaffected.
- Practical risk is low. The relevant states are nondegenerate in the M = 0, k = 0 sector, and the larger-l_max ED above agrees at 4 points.
- The fix is cheap. Run a Sylvester inertia count (sparse LDL^T of A - sigma I gives the exact number of eigenvalues below sigma), or use dense eigh where dim <= ~1e4. That covers every 015 row except L = 6 (73789) and L = 10 (184756).
- For AKLT the n = 10 dense computation already gives a Lanczos-free GM bound >= 0.30107. The n = 12 bound 0.3026 relies on Lanczos.

## Verdicts

- **C10: weakened (rewording; nothing broken).** Every recomputed number reproduces. Required corrections:
  - the Goepfert-Mack misdescription;
  - d = 2 dropped from the setting;
  - "no mechanism" scoped to the mechanisms examined;
  - E reframed as an illustration (Hartree tadpole) rather than evidence about gauge-invariant gap equations;
  - D's 7e5 attributed to the symmetric (non-minimizing) Gaussian;
  - Gribov prefactor typo in the 014 README.
- **C11: survives with minor rewording.**
  - Constants, ED gaps, PT coefficients, Hessian and the cited values are all confirmed.
  - The GM/Knabe bound is a bound on ring gaps for m > 2n (m > n), in units h = P_2.
  - The "l_max ~ sqrt(N/g)" floor scaling is a heuristic, not a proof.
  - Convergence of the infinite-chain strong-coupling series is not proved in 015. For finite chains it holds by analytic perturbation theory, since V is bounded. For the infinite chain it needs a cluster-expansion theorem, which is cited/standard.

## Replacement CLAIMS.md wording

**C10** (setting: continuum and lattice YM, d = 3, 4; toy models as stated in work 014):
> None of the interaction-gap mechanisms examined (BCS, Mott, Peierls, excitonic, dual superconductor, Gaussian/Hartree/BCS-type variational vacua, gap equations, vector large N) gives a gap lower bound for continuum YM. In condensed matter their control comes from a Fermi-surface log, a free small parameter, or a stable static background. YM has none of these at the gap scale: bosons have no log (the bosonic Hartree mass ~ sqrt(lambda) x cutoff is a counterterm-removable tadpole), the coupling at the gap scale is fixed, and the Savvidy background has the Nielsen-Olesen tachyon. Illustrations [proved/numerical]:
> - the leading-order large-N gap equation applied at N = 2 predicts m > 0 at all beta for 2D XY, but Froehlich-Spencer prove power-law decay at large beta;
> - the cutoff Hartree (4-gluon tadpole) gluon mass^2 is proportional to (d-1+xi), so it is not a pole mass;
> - the Gribov propagator's Schwinger function changes sign at gamma t = pi/(2 sqrt2), so gamma is not a spectral gap;
> - in anharmonic QM the symmetric Gaussian (Hartree) gap is 1.01x to 7.2e5x the exact gap, and the minimizing broken Gaussian predicts 0;
> - variational vacua bound E0 only.
>
> Controlled cases found, none of them continuum non-Abelian YM: vector large N (Kupiainen 1980 O(N); Kopper-Magnen-Rivasseau 1995 GN_2), Goepfert-Mack 3D compact U(1) lattice gauge theory with Villain action (pure U(1), confinement at all couplings via the dilute monopole gas; its continuum limit is a free photon), and strong-coupling lattice expansions. No rigorous 1/N mass-gap proof was found for 2D CP^{N-1} or the 't Hooft model.
>
> Status: proved/numerical (toy calculations, independently rechecked in work 018); cited (FS, KMR, GM, Kupiainen, Kennedy-Lieb checked against abstracts or records in work 018; others bibliographic only).

**C11** (setting: 1+1D O(3) quantum rotor chain, spin-1 Heisenberg and AKLT chains, finite rings/open chains, and the AKLT thermodynamic limit):
> - Rotor chain: certified ED lower bounds (floating point, Lanczos completeness assumed) reach g = 0.25 at L = 2, but only g = 2 at L = 8-10. The Q-space floor needs l_max growing with N/g [heuristic], and the dimension is exponential. Certified numbers are the finite-ring rotor tower (gap ~ g/L), not m(g) ~ e^{-2pi/g}/g.
> - AKLT (h = P_2): Gosset-Mozgunov Thm 3 with open-chain eps_10 = 0.369844 (dense ED) gives ring gap >= 0.3010 for all rings m > 20. With eps_12 = 0.364356 (Lanczos) it gives >= 0.3026 (Knabe: >= 0.3008), hence a thermodynamic-limit gap. The literature value is ~0.350 (Garcia-Saez-Murg-Wei).
> - Heisenberg and rotor chains are not frustration-free (Heisenberg frustration 0.594 per bond, N = 12). The deformation to AKLT is O(1) per bond.
> - RK parent H_FF of the rotor chain: Bakry-Emery gives gap(H_FF) >= 1 - 4beta. The Hessian bound -4 per site is sharp (Neel). This works only for beta < 1/4 (g > 2 when matched), and H_FF >= (2/g)H, so nothing transfers to H.
> - Strong-coupling series (infinite chain, x = 1/g^2): gap/g = 1 - (2/3)x + (2/27)x^2 + (7/405)x^3 + 0.0026194 x^4 + O(x^5). Independently reproduced; infinite-chain convergence radius not computed.
> - No technique reaches the O(3) continuum gap. Its known value m = (8/e) Lambda_MSbar (Hasenfratz-Maggiore-Niedermayer) rests on the TBA and is not rigorous.
>
> Status: proved (bounds as stated, floating point), numerical, cited; skeptic-checked in work 018.

## Open issues

- Inertia-count (or interval) certification of the Lanczos-based lower bounds in 015 (and C1, C5, C9) has not been done. The tool exists in scipy (sparse LU/LDL^T pivots); it needs a short script.
- Goepfert-Mack: the full text was not read to confirm it also states a mass gap (photon mass), in addition to confinement. 014 claims "photon mass^2 ~ e^{-S_0}"; that should be tagged as Polyakov's (non-rigorous) result unless GM is checked.
- The CP^{N-1} rigorous 1/N search was a single web search. A targeted search of mathematical-physics literature is needed before calling it open.
- Hainzl-Seiringer, FHNS, Nielsen-Olesen and Savvidy theorem statements were not re-read here.
- HMN O(3) paper (PLB 245, 522) not opened; only the general-N abstract (p. 529) was checked.
