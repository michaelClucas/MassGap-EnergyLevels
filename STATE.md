# State

Iteration: 2 (2026-11-20)

## Status
Iteration 2 done: data verified, iteration-1 ED claims audited, certified two-sided gap bounds, line 2 (frustration-free) assessed, lines 3 and 4 written up as obstructions. No result meets a stop condition. Lean toolchain still not installed (needs user OK to install `elan`).

Compute: CPUs 9-11 via `tools/pin.py`. Repin this iteration fixed 4 stray bash processes from subagent shells (python jobs were already pinned). A skeptic job (torus j_max = 2, ~4 GB RAM) and a certified-bounds job (~3.6 GB) ran concurrently; watch memory when two large torus jobs overlap.

## Settled
- Data (`data/`): all 808 rows verified against the papers' PDFs by an independent agent (`work/20261120-004-data-verify/`); 0 value errors. Two metadata errors fixed this iteration (logged in `data/SOURCES.md` and the 002 README erratum): AT2016 table numbers relabelled to printed B1-B12; Chen2006 arXiv date Oct 2005, journal PRD 73 (2006) 014516. Inter-source discrepancies are real (also printed by AT2021 Tabs 39-40): SU(3) m(0++)/sqrt(sigma) MP1999 +2.1 sigma vs AT2020; LTW2004 vs AT2021 SU(inf) 4.3 sigma.
  - Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
- Prior art (`literature/`, plus `work/20261120-007-*/literature-notes.md`): see iteration 1 list. New: Yarotsky (math-ph/0412040) Theorem 1 allows infinite-dimensional local spaces and covers strong-coupling KS with Gauss law. RK parent Hamiltonians for any real action are known (Castelnovo et al., Henley, Fendley-Fradkin). No paper citing 2409.10612 or 2204.12737 builds a non-Abelian FF gauge model. 2607.13228 (Gribov-horizon spectral geometry) relevant to line 3.
- SU(2) KS ED (C1-C3) audited by skeptic 005: normalization correct; C2 series extended and independently derived; C3 slope derived; ladder reproduced by independent code. Torus basis is zero-flux only; torelons lie lower.
- Certified bounds (C5, `src/certify.py`, `work/20261120-006-*/`): min-max upper + Feshbach lower bounds on the untruncated finite-lattice gap. Plaquette width 1e-11, ladder 1e-8 for g^2 >= 0.5, torus only g^2 >= 3.

## Live lines of attack (ranked)
1. Line 1 + 2 combined, stepping-stone theorem: 2+1D, free boundary, gap(H_FF) >= c min(1, 1/beta) for all beta (ground state is a plaquette product there). Needs certified single-link gap gamma_* and two-link correlation c(beta). Honest scope: H_FF has z = 2 RK dynamics and does not approach continuum YM, so at best a lattice stepping-stone; check against Result A novelty bar.
2. Line 1: certified bounds deeper into weak coupling (torus momentum blocking, j_max = 2, g^2 ~ 1; include flux sectors so the true lowest gap is bounded). Purpose: data for finite-size criteria and honest picture of truncation limits.
3. Line 6: O(3) quantum rotor chain as test bed: does any method here reproduce its known gap behavior?
4. Line 5 (interaction-driven): low priority, mostly published.

## Dead / Result-C candidates (pending skeptic)
- Line 4 (theta/Bloch), C4: blocked. Bloch decomposition gives no in-sector gap bound; band width grows with volume; band data and mass gap independent (exact 2D examples). Only gives theta- and flux-dependence on small T^3 (Koller-van Baal), not the 0++ gap.
- Line 3 (Payne-Weinberger on Gribov region), C7: blocked at including the potential (needs open log-concavity hypothesis; magnetic potential non-convex inside Omega); kinetic bound -> 0 in both limits.
- Line 2 as a route to the KS/continuum gap beyond strong coupling: KS frustration-free only with quasi-local terms (proved identity; no finite-range decomposition is conjectured); transferring H_FF gap to KS works only at strong coupling (Yarotsky already). String-net/quantum-double FF points are in the deconfined phase.

## Open questions
- Is the H_FF bound gap >= 1/2 - 2(d-1) beta correct, and is the "KS = H_FF + correction" identity right? (skeptic needed)
- Interval-arithmetic version of C5 bounds (currently floating-point).
- Gribov diameter scaling in d = 3 rests on L = 3, 4 only; diameters are lower estimates.
- AT2021 quotes error 0.15 on its Chen2006 conversion; our quadrature gives 0.10.
- Unvetted claimed full proofs (2506.00284, 2611.11534, Zenodo) remain unrelied on.

## Next planned tasks (iteration 3)
- Skeptic: audit C4 and C7 (agent 008) together.
- Skeptic: audit C5 (certified bounds lemma and code) and C6 (H_FF construction and Bakry-Emery bound).
- Mathematical physicist: attempt the 2+1D H_FF theorem gap >= c min(1, 1/beta); first check novelty against RK / Glauber-dynamics literature (gap of Langevin/heat-bath dynamics for 2D lattice gauge = product measure in 2D, likely known).
- Computational physicist: torus momentum blocking + flux sectors, certified bounds at g^2 ~ 1-2.
- Lean: blocked on user approval to install `elan`.
