# State

Iteration: 5 (2026-11-20)

## Status
Iteration 5 done. C10-C12 passed skeptics with rewording. Quantum-chemistry methods assessed (C13, pending skeptic), which closes the last technique family on the prompt's list. First Lean lemmas proved with no sorry (C2 arithmetic, C4 abstract steps, C5 E0-version Schur lemma). The program is converging on Result C. One open item could change the picture: post-2008 bound-state models of glueballs (functional Bethe-Salpeter, e.g. arXiv:2004.00415) reportedly agree well with lattice. A Bethe-Salpeter bound-state equation is the positronium-type electron-level method, so it must be checked for Result B (would be "already published", not ours) before a Result C is declared.

Lean: Lean v4.36.0 + mathlib v4.36.0 in `lean/`; `lake build` passes (4702 jobs), `#print axioms` shows only propext, Classical.choice, Quot.sound. mathlib v4.36.0 lacks Courant-Fischer (min-max), so eigenvalue-counting statements are not yet formalized.

Compute: CPUs 9-11 via `tools/pin.py`; repin found nothing to fix.

## Settled
- Data (`data/`): verified (iteration 2). Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
- Skeptic-passed claims: C1-C12. Pending: C13.
- Competitor glueball-model numbers verified against primary sources (`work/20261120-019-skeptic-c12/competitors_verified.csv`).

## Lines of attack (all assessed)
1. Line 1 (Hamiltonian lattice ED + certified bounds): numerics only, finite volume (C1-C3, C5, C9). Weak-coupling failure: truncation floor needs exponentially many states (also C11). Parked.
2. Line 2 (AKLT/Knabe/martingale/RK + Bakry-Emery): no new result (C6, C8); calibrated in C11.
3. Line 3 (configuration space, Gribov, Payne-Weinberger): obstruction (C7).
4. Line 4 (theta/Bloch): obstruction (C4), abstract part formalized.
5. Line 5 (interaction-driven): obstruction (C10).
6. Line 6 (stepping stones): calibration (C11): no technique reaches the O(3) continuum gap.
7. Quantum-chemistry methods: no lower bounds (C13, pending skeptic).
8. Result B: our constituent-gluon model refuted (C12); none of 7 published 1976-2008 models fits within lattice errors. Post-2008 functional/Bethe-Salpeter and holographic models unchecked.

## Open questions
- Do post-2008 bound-state glueball calculations (functional Bethe-Salpeter: Sanchis-Alepuz, Fischer, Huber et al.; recent holographic fits) reproduce lattice ratios within errors with fewer parameters than predictions? If yes, Result B exists in the literature (cite it; the loop may stop with the known result). Check parameter counts honestly (truncation choices, input propagators from lattice).
- Lanczos-completeness assumption in certified lower bounds (C1, C5, C9, C11): close with a sparse LDL^T inertia count.
- Goepfert-Mack full text: does it prove a photon mass, or only confinement? (C10)
- Rigorous 1/N mass generation for 2D CP^{N-1}: one web search found nothing; targeted search needed before calling it open.
- Courant-Fischer in Lean (needed for the k-eigenvalue Schur count).
- Unvetted claimed full proofs (2506.00284, 2611.11534, Zenodo) remain unrelied on.

## Next planned tasks (iteration 6)
- Literature scout + phenomenologist: post-2008 bound-state glueball models (functional Bethe-Salpeter 2004.00415 and related; holographic updates); extract their predictions, parameter counts and inputs; compare with `data/` against the Result-B bar.
- Skeptic: audit C13.
- Computational physicist: inertia-count certification of the Lanczos-based lower bounds (C5 ladder/torus, C9, C11 rotor).
- Overseer/drafter: draft the consolidated Result-C claim file (`work/<id>-resultC-draft/`), one row per technique with failure point and evidence, and extend `src/reproduce.py` so every number in it regenerates.
- Formalizer: the missing half of min-max from `eigenvectorBasis`, then the k-eigenvalue Schur count for C5.
