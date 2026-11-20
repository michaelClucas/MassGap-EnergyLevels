# State

Iteration: 6 (2026-11-20)

## Status
Iteration 6 done. All technique families on the prompt's list are assessed and every claim C1-C13 has passed a skeptic. No published electron-level glueball model meets Result B (pre- and post-2008, 12 models checked). A candidate Result C draft exists: `work/20261120-025-resultC-draft/claim.md` with `number_map.csv`. Not yet RESULT.md: the stop conditions still need (1) a fresh skeptic that attacks the consolidated claim file for one full iteration and fails, and (2) the full `python -m src.reproduce` run end to end.

Reproduce: `python -m src.reproduce --quick` passes (38 s, 29/29 quick numbers matched, rechecked by the overseer). Full mode untested end to end; estimated 2-2.5 h at 1 thread (019 `fit.py` alone ~75 min). 7 numbers are "recorded" (checked against committed files, not regenerated: torus j_max = 2 builds ~5 GB, rotor certification > 30 min).

Lean: `lake build` passes (5500 jobs), no sorry, axioms propext/Classical.choice/Quot.sound only. Proved: C2 strong-coupling arithmetic; C4 abstract Bloch steps; C5 Schur lemma, E0 version and k-eigenvalue count (via the easy half of Courant-Fischer, `lean/MassGap/MinMax.lean`). Not formalized: physics inputs (matrix elements, floors), full Courant-Fischer.

Compute: CPUs 9-11 via `tools/pin.py`; repin found nothing to fix.

## Settled
- Data verified (iteration 2). Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
- Claims C1-C13 skeptic-passed (see CLAIMS.md). Certified lower bounds in C5, C9 rechecked by inertia counts (work 024); C11 partly.
- Correction logged: C5 ladder gap at g^2 = 1 is 3.0620133 (006 README had a transcription error 3.0620125; results.csv is correct).
- Result B: none of the 12 models checked (7 pre-2008, 5 post-2008 incl. functional Bethe-Salpeter) fits within lattice errors. Best post-2008: functional BSE (Huber-Fischer-Sanchis-Alepuz), 0 fitted parameters, 1/9 within 1 lattice sigma, chi^2/n 104 (work 022).

## Lines of attack: all assessed (see claim draft table)
1. Bound-state ED + certified bounds (C1-C3, C5, C9): finite-volume numerics; fails at weak coupling (truncation floor).
2. Spectral-gap methods AKLT/Knabe/martingale/RK + Bakry-Emery (C6, C8, C11): strong coupling only or not new.
3. Configuration space / Gribov / Payne-Weinberger (C7): obstruction.
4. theta/Bloch (C4): obstruction.
5. Interaction-driven (C10): obstruction.
6. Stepping stones (C11): toolkit fails on the O(3) chain for the same reasons.
7. Quantum chemistry (C13): no lower-bound machinery beyond certified ED + Schur.
8. Result B (C12 + work 022): not met.

## Open (listed in the draft as open, not blocking)
- Explicit beta-dependence of the 2D Langevin gap (C8 conjecture Theta(1)).
- Rigorous 1/N mass generation for 2D CP^{N-1}.
- FBSE N-dependence and C = - glueballs (not computed in the literature).
- A weak-coupling comparison operator for Temple/Schur (none known, none ruled out).
- Goepfert-Mack full text (photon mass vs confinement only).
- Unvetted claimed full proofs (2506.00284, 2611.11534, Zenodo) remain unrelied on.

## Next planned tasks (iteration 7)
- Result-C skeptic A (technical): fresh agent given only claim.md and the repo; tries to break every row (hidden assumptions, overstated scope, finite-size effects, bounds that depend on lattice spacing, misread data) for one full iteration.
- Result-C skeptic B (literature/novelty): checks the literature-check section and that every "blocked" statement is not contradicted by a published result.
- Full reproduce run in the background (`python -m src.reproduce`, pinned, logged), ~2.5 h.
- If both skeptics fail to break it and full reproduce passes: write RESULT.md, commit, push, notify, stop the loop.
