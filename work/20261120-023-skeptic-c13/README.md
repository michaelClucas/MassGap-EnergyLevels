# 20261120-023 Skeptic audit of C13 (quantum-chemistry methods)

## Question

Does C13 (`work/20261120-020-quantum-chemistry-methods/`) hold up: HUM logic, the selected-CI underestimate, the HF statements, CC/EOM-CC, the MPS "non-monotone gap error", the area-law counterexample, Temple circularity, and the four citations written from memory?

Setting of every KS number: 2+1D Hamiltonian lattice, SU(2), a = 1, finite lattice (single plaquette, open 2-plaquette ladder), fixed g^2, H = (g^2/2) sum_l j_l(j_l+1) - (2/g^2) sum_p Tr U_p. Spin-chain numbers: finite open/periodic 1D chains.

## Method

- Reran `qc_methods.py` unchanged (copy in `rerun/` with only the repo-root path depth adjusted; 1 thread, pinned). Compared all CSVs with 020's.
- `skeptic_c13.py` (rerun: `python tools/pin.py .venv/Scripts/python.exe work/20261120-023-skeptic-c13/skeptic_c13.py`, 1 thread, about 10 min, < 1 GB). Output `summary.txt`, `run.log`, `hf_check.csv`, `ci_selected.csv`, `mps_check.csv`, `chain_nondegenerate.csv`:
  - A. HF: builds f(U_p1) f(U_p2) pointwise (f = plaquette ground state, tj <= 6), fits it in the ladder spin-network basis (tj <= 12) from 3x-oversampled random link configurations, evaluates <H_ladder>.
  - B. Selected CI without exact eigenvectors (two-root CIPSI with Epstein-Nesbet importance, from the vacuum), state-averaged selection (|v0|^2 + |v1|^2), separate state-specific spaces, and 020's E1-only selection.
  - C. MPS: degeneracy diagnostic of 020's Schmidt weights; 020's recipe with deterministic tie-breaking; a self-consistent DMRG-style run (environments from the Ritz vectors of the previous step, iterated to a fixed point, no exact states used).
  - D. Non-degenerate D = 1 gapless chain.
  - E. Sign patterns of 020's CC/EOM-CC table.
- Citations checked against Crossref and arXiv abstract pages on 2026-11-20.

## Verdict per sub-claim

**(1) HUM and "no bound in either direction": correct as stated, but incomplete in a way that matters.** Delta(V) = Delta + (E_1(V) - E_1) - (E_0(V) - E_0) with both brackets >= 0 (proved). But two Ritz values combined with one certified bound from elsewhere do give rigorous gap bounds (proved, one line each):
- L_1 <= E_1 certified  =>  Delta >= L_1 - E_0(V) (lower bound; the Ritz E_0(V) is the needed upper half).
- L_0 <= E_0 certified  =>  Delta <= E_1(V) - L_0 (upper bound).

C5 is exactly the first line: Ritz E_0 plus a Schur/Feshbach lower bound on E_1. The Feshbach/Schur step is itself a quantum-chemistry technique (Löwdin partitioning, Löwdin, J. Math. Phys. 3, 969 (1962) [Crossref]). So "none gives a gap lower bound" is too strong. The accurate statement: CI/HF/MPS Ritz values supply the upper half of a rigorous gap lower bound; the lower bound on E_1 must come from a comparison argument (partitioning/intermediate problems/Temple with certified beta), which in KS is C5's electric floor.

**(2) Selected-CI underestimate: reproduced, but a straw man.** Rerun reproduces 0.045263 vs 3.062013 at ladder g^2 = 1, n = 12 and all sign changes (`ci_truncation.csv` agrees to 2e-13). The construction keeps the n largest components of the exact E_1 eigenvector and reads E_0 from the same space; E_0(V) is then far above E_0 (-0.045 vs -3.348), which is the whole underestimate. No practitioner computes E_0 in a space selected for E_1. Honest selections, same lattice and couplings, sizes 4-64 (`ci_selected.csv`):

| g^2 | two-root CIPSI (no exact input) | state-averaged | state-specific (E_0 from E_0 space, E_1 from E_1 space) | 020 E1-only |
|---|---|---|---|---|
| 0.5 | +1.21 ... +0.141, all > 0 | +4.26 ... +0.062, all > 0 | +2.12 ... +0.065, all > 0 | +2.75 ... -2.96 ... -0.25 |
| 1.0 | +0.076 ... +3.8e-5, all > 0 | +0.78 ... +2.5e-5, all > 0 | +0.22 ... +2.5e-5, all > 0 | -3.02 at n = 12 |
| 2.0, 3.0 | all > 0 | all > 0 | all > 0 | all > 0 at these n |

So every reasonable CI selection tested overestimates the gap here (observation, finite lattice). Non-monotonicity survives in CIPSI (g^2 = 0.5: +1.21, +1.70, +0.79, +1.32, ...; g^2 = 1: +0.076, +0.19, +0.30, +0.21, ...). The structural point stands (HUM gives no sign), but the numerical "underestimates by up to 3.0" must be labelled as coming from an E_1-only selection.

**(3) HF claims: both correct (proved here independently), with one scope fix.**
- Link product state: if prod_l phi_l(U_l) is invariant under a gauge transformation g at a vertex x, then for each link at x the transformed factor equals c_l(g) phi_l (a product vector equals another product vector only factorwise up to scalars). c_l is a continuous unitary character of SU(2), hence trivial (SU(2) is perfect), so phi_l is invariant under left (or right) translation and is constant by transitivity. Every link has an endpoint, so the state is the all-j = 0 strong-coupling vacuum. Proved.
- E_HF = 2 E0(plaq): checked pointwise, not just asserted (020's script sets E_HF := 2 E0(plaq) without evaluating it). Fit residual <= 2.4e-15, norm 1 to 12 digits, <H_ladder> - 2 E0(plaq, tj <= 6) between -5.3e-15 and 0 at all four g^2 (`hf_check.csv`). Identity holds at every truncation, so it is algebraic, as the proof says. The cross term vanishes because sum_c (d_c/(d_a d_b)) C_c = C_a + C_b.
- Scope fix: "best gauge-invariant product state" depends on the factorization. It is the best product over the two plaquette holonomies in a maximal-tree gauge with base point on the shared link; another choice of loop variables (e.g. plaquette 1 and the outer loop) gives a different "HF". Symmetry-projected HF is not a product state and was not tested.

**(4) CC / EOM-CC: correct with wording fixes.** E_CC < E_0 at n = 1, 3, 5, 7, 9 at all g^2 where nonzero, and also at even n (g^2 = 0.5: n = 10; g^2 = 1: n = 8, 10). EOM-CC gap error does not strictly alternate: signs `+----+-+-+` (g^2 = 0.5), `---+-+-+--` (g^2 = 1), `-+-+-+-` (g^2 = 2, 3). Replace "alternates" with "takes both signs". CC numbers reproduce bit-for-bit.

**(5) MPS: 020 numbers not reproducible; the qualitative claim survives honest targeting.**
- Rerun with 1 thread changes `mps_ladder.csv` by up to O(1) (e.g. g^2 = 0.5, E1-target, D = 8: -2.20 in 020, +2.61 now; g^2 = 2 averaged D = 8: -5.0e-3 in 020, +9.2e-3 now). Cause: the E1-targeted density matrix has exactly degenerate Schmidt pairs (relative spacing 6e-16 at g^2 = 0.5 and 2), and truncating at D splits a degenerate pair, so the kept vector is chosen by rounding noise. A real DMRG keeps whole degenerate multiplets. The averaged-target g^2 = 0.5 row quoted in 020 does reproduce.
- Self-consistent DMRG-style run (no exact states): gap error still non-monotone in D and of both signs. g^2 = 0.5 averaged: -2.17, -2.17, +3.31, +3.31, +0.89, +1.19, +0.75, +0.17 (D = 2..24); g^2 = 1 E1-target: ..., +0.175, -0.56, -1.12, -0.32; g^2 = 2 averaged: ..., +0.143, +1.5e-3, -4.3e-6, +8.5e-8. So the claim survives, with new numbers.

**(6) Area-law converse: correct, cited correctly after one fix.**
- Ferromagnetic Heisenberg ring: product state exact, finite-L gap above the ground multiplet 1 - cos(2 pi/L) (magnon epsilon(k) = 1 - cos k ~ k^2/2). Correct, but the ground space is (L+1)-fold degenerate, so "gap" means gap above the multiplet. A non-degenerate example removes the objection: H = sum_i (1/4 - S_i.S_{i+1}) + (1/2 - S^z_1) on an open chain. All terms are PSD, |up...up> is the unique zero-energy state (kernel = symmetric subspace with site 1 up), D = 1, and the gap is 1 - cos(pi/(2L+1)) -> 0: ED gives gap (2L+1)^2 = 5.76, 5.55, 5.43, 5.35, 5.29 at L = 6-14 (-> pi^2/2 = 4.93) (`chain_nondegenerate.csv`). The critical TFIM is also non-degenerate with poly(L) bond dimension (Verstraete-Cirac), so it is the more relevant counterexample; it was already in 020.
- Hastings, J. Stat. Mech. (2007) P08024, arXiv:0705.2024 [arXiv]: area law for gapped local 1D systems and MPS approximability. Correct.
- 020 cites Arad-Kitaev-Landau-Vazirani arXiv:1301.1162 [arXiv] as giving "bond dimension poly(L) (exp(O(1/Delta)) prefactors)". AKLV actually prove S = O(log^3 d / Delta) and a sublinear bond dimension exp(O(log^{3/4} n / Delta^{1/4})). Fix the wording. The polynomial-time algorithm with poly bond-dimension MPS output is Landau-Vazirani-Vidick, arXiv:1307.5143 [arXiv] (published Nature Physics 11, 566 (2015) [memory, not checked]); add it if the poly(L) statement is kept.

**(7) Temple: 3 of 88 verified; "only non-circular route" not established.** Rerun gives the same 3/88 false with beta = CI E_2 and 0/88 with exact E_2 (88 = 2 lattices x 4 couplings x k = 2..12). The bound used, E_1 >= eta - sigma^2/(beta - eta) for E_1 < beta <= E_2, follows from <(H - E_1)(H - beta)> >= 0 (proved). The CLAIMS row says only that C5's floor closes the circle, which is fine. The README's sentence "no better comparison operator is known at weak coupling" is fair, but other standard closures exist (Weinstein-Aronszajn / Bazley-Fox intermediate problems, Lehmann-Goerisch with a base problem; Bazley-Fox, Phys. Rev. 124, 483 (1961) [memory, Crossref rate-limited]). Each needs a comparison operator H >= H_0 with a computable spectrum. The C5 floor (H_E - ||H_B|| on the complement) is one such H_0. It has not been shown that every H_0 fails at weak coupling. The claim must not say or imply "only".

**(8) Citations from memory.**
- Kato, "On the upper and lower bounds of eigenvalues", J. Phys. Soc. Japan 4, 334-339 (1949), doi:10.1143/JPSJ.4.334 [Crossref]. Verified (Kato-Temple inclusion bounds; content as used).
- Lehmann: ZAMM 29, 341-356 (1949), doi:10.1002/zamm.19502911005, and part II ZAMM 30, 1-16 (1950), doi:10.1002/zamm.19500300101 [Crossref]. Maehly, Helv. Phys. Acta 25, 547 (1952): no Crossref record. Secondary sources confirm Lehmann and Maehly independently gave the lower-bound method, but the primary paper was not verified. Mark Maehly [secondary].
- Fannes-Nachtergaele-Werner, "Finitely correlated states on quantum spin chains", CMP 144, 443-490 (1992), doi:10.1007/BF02099178 [Crossref]. Verified bibliographically; the parent-Hamiltonian gap for pure (injective) finitely correlated states is the content used.
- Bartlett-Musiał, "Coupled-cluster theory in quantum chemistry", Rev. Mod. Phys. 79, 291-352 (2007), doi:10.1103/RevModPhys.79.291 [Crossref]. Verified; non-variational CC is standard.

## Replacement CLAIMS.md wording (C13)

| C13 | quantum-chemistry methods (HF, CI, coupled cluster, EOM-CC, DMRG/MPS, Temple/Kato) | SU(2) KS single plaquette and open ladder, 2+1D, finite lattice, g^2 = 0.5-3; 1D spin chains; abstract statements as cited | None of these methods gives a gap lower bound by itself. CI/HF/MPS Ritz values bound each eigenvalue from above (Hylleraas-Undheim-MacDonald), so the Ritz gap alone has no bound in either direction; with a certified lower bound L_1 <= E_1 from elsewhere, Delta >= L_1 - E_0(V) is rigorous, which is how C5 works (Ritz E_0 + Schur/Löwdin-partitioning bound on E_1). On the ladder every CI selection tested that does not use the exact E_1 alone (vacuum-depth, j_max, two-root CIPSI, state-averaged, state-specific) overestimates the gap, non-monotonically for CIPSI (observation); an E_1-only selection underestimates by up to 3.0 (0.045 vs 3.062 at g^2 = 1) and changes sign. A gauge-invariant product over links is the strong-coupling vacuum (proved); the best product over the two plaquette holonomies (maximal-tree gauge, base point on the shared link) has E_HF = 2 E0(plaquette) (proved, checked pointwise to 5e-15). CC energy falls below E0 at every odd truncation n <= 9 where the error exceeds 1e-12, and at some even ones; EOM-CC gap error takes both signs. MPS gap error is non-monotone in bond dimension and takes both signs, also with self-consistent DMRG-style targeting (no exact-state input). Area law: gap => efficient MPS in 1D (Hastings; Arad-Kitaev-Landau-Vazirani), not the converse (critical TFIM; open ferromagnetic Heisenberg chain with boundary field: unique D = 1 ground state, gap ~ pi^2/(2(2L+1)^2)). Temple lower bounds need an a priori beta <= E_{k+1}: with CI values as beta, 3 of 88 'bounds' are false. Closing the regress needs a comparison operator with computable spectrum (C5's electric floor is one; it fails at weak coupling); no weak-coupling comparison operator is known, none is ruled out. | proved (structural statements, HF identities, Temple logic, D = 1 example); numerical (finite lattice); cited (HUM, Kato 1949, Lehmann 1949, Löwdin 1962, Fannes-Nachtergaele-Werner 1992, Bartlett-Musiał 2007, Hastings 2007, AKLV 2013: checked; Maehly 1952 via secondary sources); skeptic passed with rewording (023) |

## Open issues

- 020's `mps_ladder.csv` is not reproducible at the quoted digits (degenerate Schmidt pairs split by rounding). 020's README table should cite `mps_check.csv` here (deterministic, self-consistent) instead.
- 020's `qc_methods.py` asserts E_HF = 2 E0(plaq) without evaluating it; `skeptic_c13.py` part A is the evaluation.
- Symmetry-projected HF and other loop-variable factorizations not tested (do not affect the no-lower-bound conclusion).
- Bazley-Fox and LVV journal references not Crossref-checked (rate limit); Maehly 1952 primary not found online.
- Whether an intermediate-problem comparison operator better than H_E - ||H_B|| exists at weak coupling is open. It is the only route this family leaves, and it is the C5/C9/C11 question.
- `rerun/` holds the 1-thread rerun of 020's script (path depth changed only) and its outputs, used for the comparison above.
