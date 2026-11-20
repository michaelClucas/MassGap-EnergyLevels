# State

Iteration: 4 (2026-11-20)

## Status
Iteration 4 done. Every line of attack and the first Result-B attempt have now come back negative. The program is converging on Result C (obstruction covering every technique examined). Remaining before a Result C can stop the loop: skeptic on C10-C12, quantum-chemistry methods (HF, CI, CC, DMRG/tensor networks) assessed, a consolidated claim file, `python -m src.reproduce` covering every number in it, and Lean lemmas for key steps.

Lean: `elan` 4.4.0 installed (user approved 2026-11-20). `lean/` pinned to Lean v4.36.0 and mathlib v4.36.0; smoke test builds with `lake build` (`work/20261120-017-lean-setup/`). No real lemmas formalized yet.

Compute: CPUs 9-11 via `tools/pin.py`; repin found nothing to fix.

## Settled
- Data (`data/`): verified (iteration 2). Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
- Prior art: `literature/` plus work 001, 007, 011, 013, 014. Rigorous lattice gap only at strong coupling (Osterwalder-Seiler, Shen-Zhu-Zhu, Cao-Nissim-Sheffield, Yarotsky for KS); Piazza-Shen 2610.35286 volume-uniform LSI for 2D SU(N).
- Skeptic-passed: C1-C9 (with rewording).

## Lines of attack (all assessed)
1. Line 1 (Hamiltonian lattice ED + certified bounds): numerics only, finite volume (C1, C2, C3, C5, C9). Certified bounds fail at weak coupling because the truncation floor needs exponentially many states (same failure on the rotor chain, C11). Parked.
2. Line 2 (AKLT/Knabe/martingale/RK + Bakry-Emery): no new result (C6, C8). Works only at strong coupling or for frustration-free deformations with z = 2 dynamics. Calibrated on AKLT (works) and rotor/Heisenberg chains (fails for the same reason) in C11.
3. Line 3 (configuration space, Gribov, Payne-Weinberger): obstruction (C7).
4. Line 4 (theta/Bloch): obstruction (C4).
5. Line 5 (interaction-driven): obstruction (C10, pending skeptic). Possible side question: rigorous 1/N mass generation in 2D CP^{N-1} (not checked against literature; belongs to line 6).
6. Line 6 (stepping stones): toolkit calibration (C11, pending skeptic). No technique reaches the O(3) continuum gap; its known value comes from integrability.
7. Result B: first attempt refuted (C12). Constituent-gluon model is a known class and fits worse than Kaidalov-Simonov. No electron-level model in the literature fits within lattice errors.

## Open questions
- C10, C11, C12 need a skeptic; C12's competitor numbers (Mathieu-Buisseret-Semay Table 3, Kaidalov-Simonov Table 6) need verification.
- Quantum-chemistry methods (HF, CI, coupled cluster, DMRG/tensor networks) not yet assessed as gap-proof techniques; required for Result C to cover every technique in the prompt.
- Lambda_MSbar r0 not in `data/` (needed for kappa = m/Lambda in C10).
- Interval-arithmetic certified bounds (currently floating point beyond the plaquette).
- Unvetted claimed full proofs (2506.00284, 2611.11534, Zenodo) remain unrelied on.

## Next planned tasks (iteration 5)
- Skeptic: audit C10 (line 5) and C11 (line 6).
- Skeptic/phenomenologist: audit C12, verify competitor-model numbers from primary sources.
- Mathematical physicist: quantum-chemistry methods (HF, CI, CC, DMRG/tensor networks) as gap techniques: what each gives (upper bounds, variational, no lower bound), where they fail for YM; cite prior uses on lattice gauge theory.
- Formalizer: first Lean lemmas: (a) abstract Bloch-decomposition step of C4 (direct sum of operators: spectrum is the union, no in-sector gap bound from band data) or min-max upper bound; (b) Feshbach/Schur-complement lower-bound lemma of C5 in finite dimensions. Report precisely what mathlib lacks.
- Overseer: after these, draft the consolidated Result-C claim file and extend `src/reproduce.py` to cover every number.
