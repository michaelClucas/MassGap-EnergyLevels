# State

Iteration: 3 (2026-11-20)

## Status
Iteration 3 done: skeptics passed C4-C7 with corrections; line 2's 2D target found already known in qualitative form (Piazza-Shen); torus gap computed in all flux sectors with certified bounds down to g^2 = 1.5 (zero flux). No result meets a stop condition. Lean toolchain not installed (needs user OK to install `elan`).

Compute: CPUs 9-11 via `tools/pin.py`; repin found nothing to fix this iteration. RAM: torus zero-flux j_max = 2 build peaks ~5.4 GB; flux-sector j_max = 2 needs more; don't run two such jobs at once.

## Settled
- Data (`data/`): verified (iteration 2). Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
- Prior art: see `literature/` and work 001, 007, 011. Key: rigorous lattice gap only at strong coupling (Osterwalder-Seiler, Shen-Zhu-Zhu, Cao-Nissim-Sheffield); Yarotsky covers strong-coupling KS with Gauss law; Piazza-Shen 2610.35286 gives a volume-uniform LSI for 2D lattice SU(N) at all beta.
- SU(2) KS small-lattice spectra and certified bounds (C1, C2, C3, C5, C9). Finite-volume, fixed-spacing only.
- C6: the RK-type H_FF and its Bakry-Emery gap are correct but not new.

## Lines of attack
1. Line 1 (Hamiltonian lattice, ED + certified bounds): alive as numerics only. Small tori are dominated by torelons and geometry; no continuum content. Certified bounds fail at g^2 <= 1 on 2x2 because the floor d drops below E0. Further pushes cost a lot of RAM for little information. Status: parked unless a later task needs specific bounds.
2. Line 2 (AKLT/Knabe/martingale/RK): effectively dead for new results. 2D qualitative gap known (Piazza-Shen); KS is frustration-free only with quasi-local terms; H_FF -> KS transfer only at strong coupling (already Yarotsky); H_FF dynamics (z = 2) has no YM continuum limit. Possible small result: explicit beta-dependence of the 2D gap (C8 conjecture Theta(1)), but it is a statement about Langevin dynamics, not Yang-Mills.
3. Line 3 (configuration space / Gribov / PW): Result-C obstruction, skeptic passed (C7).
4. Line 4 (theta/Bloch): Result-C obstruction, skeptic passed (C4).
5. Line 5 (interaction-driven: Mott/BCS/dual superconductor/gap equations): untested by the team; mostly published as models. Needs an obstruction write-up or a Result-B attempt.
6. Line 6 (stepping stones: O(3) rotor chain, Haldane, 2D sigma model): untested. Use it to check which electron-level proof techniques reproduce known gaps and where they fail.
7. Result B (quantitative model): not yet attempted. Data is ready.

## Open questions
- C8 (011) and C9 (012) need a skeptic.
- Is a Result C covering every technique the realistic end state? It needs lines 5 and 6 covered and the quantum-chemistry methods (HF, CI, CC, DMRG) assessed.
- Interval-arithmetic version of certified bounds (currently floating point beyond the plaquette).
- Unvetted claimed full proofs (2506.00284, 2611.11534, Zenodo) remain unrelied on.

## Next planned tasks (iteration 4)
- Phenomenologist: Result-B attempt. Electron-level bound-state model (e.g. two- or three-constituent-gluon Schrödinger/Salpeter in a linear + Coulomb potential, or a band/Mott picture) fit to glueball ratios with honest parameter count; compare with flux-tube, bag, constituent-gluon, holographic models on the same data (4D SU(3), SU(inf), 2+1D SU(3)).
- Mathematical physicist: line 5 assessment: mean-field (BCS/Mott/gap-equation) gap mechanisms; where error control fails (no small parameter, gauge dependence of the gluon mass); obstruction or open sub-question.
- Mathematical/computational physicist: line 6 stepping stone: O(3) quantum rotor chain (and spin-1 Heisenberg/AKLT) as a test bed; which of Temple/certified ED, Knabe, Bakry-Emery, variational bounds reproduce the known gap, and what blocks each from reaching the continuum O(3) sigma-model limit.
- Skeptic: audit C8 and C9.
- Lean: blocked on user approval to install `elan`.
