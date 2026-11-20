# 20261120-028: technical skeptic of the Result-C draft (025)

## Question
Does the candidate Result C in `work/20261120-025-resultC-draft/claim.md` survive a technical attack: are the per-row blocks demonstrated, does it cover every technique examined and every concept in `OVERSEER_PROMPT.md`, do its numbers trace and recompute, does the Lean statement coverage match, and is the wording safe?

## Method
- Read `CLAIMS.md`, `STATE.md`, `OVERSEER_PROMPT.md`, `claim.md`, `number_map.csv`, and the READMEs of 001, 006, 008, 009, 011, 012, 014, 015, 016, 018 where a row depended on them.
- Number trace: 13 rows of `number_map.csv` drawn with `random.seed(28); random.sample(range(69), 13)`, each match string looked up in the committed output. 6 numbers recomputed with independent code (below).
- `--quick` reproduce: run once in the repo as briefed (see issue M5, it collided with the running full reproduce), then rerun on a copy of `src/ data/ work/ tools/` in a scratch folder outside the repo (no `artifacts/` sharing).
- Lean: read `lean/MassGap.lean` and every statement in `lean/MassGap/*.lean`; grepped for `sorry`/`admit`/`axiom`; compared `.olean` times to sources. Not rebuilt.

## Verdict
**Survives with required changes.** No calculation in the table was found wrong, and every sampled number traces and the recomputed ones agree. But the headline says every technique is "blocked ... shown by explicit calculation or counterexample", and several rows only show failure of our instantiation, or rest on a heuristic, an assertion, or an open hypothesis. One row misstates a result (U(1)), three printed intervals are rounded inward against the claim's own "holds at printed precision", and the Feynman / Karabali-Kim-Nair configuration-space gap named in the prompt is neither examined nor listed. All fixable by rewording; none needs new computation.

## Issues, by severity

### High
**H1. Headline overstates the strength of the blocks.** `claim.md` Answer: "each is blocked from the continuum by an identified structural property, shown by explicit calculation or counterexample." Row by row:
- Row 1: the explicit counterexample (N15, ladder at g^2 = 0.2) shows the *crude/split electric floor* fails at weak coupling. Row 6 itself says no other weak-coupling comparison operator is ruled out. The scaling "j_max ~ g^-2" has no derivation anywhere in the repo (see M2).
- Row 2: "KS is not frustration-free" is proved only for the natural decomposition, conjectured otherwise; "(d) RK dynamics ... not Lorentz-invariant YM" is asserted.
- Row 3: block (a) is an *open* hypothesis (log-concavity of J psi_0^2, Moncrief/Mondal), so the route is open, not blocked; (d) holds only for the diameter (kinetic-only) form, see M1.
- Row 5: argued, with toy counterexamples; "no free small parameter" is false as written for SU(N), see M4.
- Row 7: "finite radius" asserted, see L2.
Fix: add a column "block type" with values P (obstruction proved for the technique as applied to YM), S (shown for our instantiation only), A (argued / heuristic), O (open / conditional). My reading: 1 = S + A; 2 = P (no gap transfers from H_FF to KS; Bakry-Emery only at strong coupling) + conjectured (no FF form of KS); 3 = O (conditional on log-concavity) + S (kinetic-only bound -> 0, continuum small T^3); 4 = P; 5 = A with toy counterexamples; 6 = P (upper-bound structure) + O (comparison operator); 7 = S + A; 8 = numerical (Result B not met; not a technique block). Replace the Answer's "each is blocked ... shown by explicit calculation or counterexample" with: "each failure is shown by explicit calculation or counterexample for the route as the program instantiated it; rows marked A or O name what would have to change for the route to work."

**H2. Row 2 misstates the U(1) result.** `claim.md`: "U(1) analogue of the all-beta theorem is false: gap(H_FF) <= 1.2014e-3 at beta = 6". The "all-beta theorem" in that row is Piazza-Shen (volume-uniform gap > 0 at every beta). The U(1) bound is positive at every beta (~2 pi beta^2 e^{-2 beta}), so it does not refute volume-uniform positivity. What 011 proves is that the *rate* c min(1, 1/beta) (target T) fails for U(1) (011 item 3; 016 (c); CLAIMS C8c, which words it correctly). Settings also differ: Piazza-Shen is the torus, the U(1) bound is free-boundary boxes. Fix: "The U(1) analogue of the rate c min(1, 1/beta) is false: on free-boundary boxes inf_L gap(H_FF) <= 1.2014e-3 at beta = 6 (~2 pi beta^2 e^{-2 beta}) [N19]; volume-uniform positivity for U(1) is neither shown nor refuted."

**H3. Coverage gaps against the prompt.**
- Feynman 1981 and Karabali-Kim-Nair (2+1D YM as configuration-space QM, gap m = e^2 N / 2 pi from the measure) are named in the prompt's "Honest framing" and summarized in 001, but `claim.md` neither examines them nor lists them under "did NOT examine" (grep: 0 hits). Fix: add to "did NOT examine": "Feynman's bounded-configuration-space heuristic and the KKN 2+1D gauge-invariant-variable gap (leading order, not rigorous): summarized in 001, not analysed."
- The Moncrief / Mondal log-concavity hypothesis (conditional Bakry-Emery on the orbit space), on which row 3 depends, is absent from "What remains open". Fix: add it, with the note that a strongly log-concave J psi_0^2 would give a diameter-free bound, so row 3(d) does not block that variant.

### Medium
**M1. Row 3 drops the Omega-vs-Lambda correction from skeptic 009.** 009 ("Direction error: Lambda vs Omega") showed the lattice diameters are of Omega, Lambda is a subset of Omega, so pi^2/D(Omega)^2 is not an upper bound on the Lambda quantity the lemma uses; only the continuum T^3 chord gives a correctly directed bound. Row 3 still says "(d) diam grows with L, so the kinetic bound vanishes as L -> inf and a -> 0" and quotes N26/N30 without the caveat. Fix: "(d) [continuum small T^3 only] diam(Lambda) >= 2 pi sqrt(3L) gives Delta_PW <= g^2/(24L) [N29]; lattice numbers N26, N30 are for Omega (contains Lambda) and are not upper bounds on the Lambda constant."

**M2. "j_max ~ g^-2" (row 1) is unsupported and probably has the wrong form.** No README derives it. With the crude floor of 006 (d = (g^2/2) c_Q - 4 N_p/g^2) and E1 = -4 N_p/g^2 + O(N_p) (zero-point energy, lattice units), clearing E1 needs c_Q ~ N_p/g^2, i.e. j_max ~ sqrt(N_p)/g. This is the same form 015 uses for the rotor (l_max ~ sqrt(N/g) in its normalization). The volume factor is what makes the dimension exponential, and the claim leaves it out. Fix: "needs j_max >~ sqrt(N_p)/g with the norm-bound floor (heuristic), so dim P ~ j_max^{#links}".

**M3. Three printed intervals are rounded inward.** The row 1 status line says "every quoted interval holds at its printed precision". From `016/run_c9.log`: Delta_all in [2.0292728, 2.0342348], printed [2.02927, 2.03423] (upper end rounded down). Delta_00 in [4.1279338, 4.1280120], printed [4.12793, 4.12801] (upper end rounded down). gap_00 >= 0.769384, printed 0.7694 (lower bound rounded up; same in CLAIMS C9). Fix: [2.02927, 2.03424] (as CLAIMS C9 already has), [4.12793, 4.12802], >= 0.7693.

**M4. Row 5 omits 1/N.** "YM has none ... a free small parameter" is false as written: 1/N is free for SU(N). 014 handles it correctly ("YM is a matrix model, planar limit not a saddle point"; single-matrix QM is the only fermion reduction). Fix: add "1/N is free, but the N = inf limit is the planar sum governed by loop equations, not a gap equation of collective fields; rigorous large-N YM results are strong-coupling only (014 sec. 6, 018)."

**M5. reproduce.py is not safe to run concurrently, and my quick run damaged the running full run.** `prepare()` does `shutil.rmtree(artifacts/<folder>)` for every task folder, and quick and full share `artifacts/`. The briefed command `python -m src.reproduce --quick` in the repo deleted the full run's (027) sandbox outputs for 002, 004, 005, 007, 008, then crashed with `PermissionError` on `artifacts/20261120-010-skeptic-certify-hff/hff_checks.log` (held open by 027). The briefing asked for both "run --quick" and "do not touch artifacts/", which conflict for this script. I did not touch artifacts/ after that. Consequence: 027's final number check will report N1-N5, N9, N31 as "file missing". Fix: after 027 finishes, `python tools/pin.py .venv/Scripts/python.exe -m src.reproduce --only 20261120-002 20261120-004 20261120-005 20261120-007 20261120-008` then `... -m src.reproduce --check-only`. Code fix: a lock file in `artifacts/` that `main()` takes before `prepare()` and refuses to run if held. The quick run itself passes when isolated: on a scratch copy, 15 scripts exit 0 in 69 s, number check 29 ok, 0 failed (`quick.log` in this folder is the failed in-repo attempt).

### Low
- **L1.** Row 5: "Froehlich-Spencer prove m = 0 there" (there = beta_XY = 2). FS give a non-explicit beta_0; 014 says no explicit beta_0 is claimed. Fix: "the gap equation predicts m > 0 at every beta (1.0563e-2 at beta_XY = 2), while FS prove m = 0 for all sufficiently large beta."
- **L2.** Row 7: "strong-coupling series has finite radius" contradicts C11 ("convergence radius not computed"). Fix: drop it or mark conjectured.
- **L3.** Row 2 (c): "the Hessian constant -8(d-1) is attained [N17]": N17 is d = 2 only; d = 3 is listed as open. Say "attained in d = 2".
- **L4.** Row 2: compares the volume-uniform bound at beta = 1000 with the single plaquette at beta = 200. Quote both at one beta.
- **L5.** Row 2 (d) "z = 2 ... not Lorentz-invariant YM": mark as argued.
- **L6.** `CLAIMS.md` C5 still says the k-eigenvalue Schur count is not formalized; `lean/MassGap/MinMax.lean` has `schur_count` and `schur_count_sorted` (k = 1 gives the E1 bound). Overseer to update.
- **L7.** The prompt's "Kogut-Susskind as tight-binding on the group manifold" framing is never named. Row 1 is that model (single plaquette in the character basis is a tight-binding chain j -> j +- 1/2 with on-site 2 g^2 j(j+1)). One sentence fixes it. Literature summary for row 2 omits Gorantla-Huang 2025 (martingale gap for a FF deformation of 3+1D Z2 gauge theory, in 001), the closest prior art, and Lemm-Mozgunov.
- **L8.** The Gosset-Mozgunov prefactor 5(n^2+n)/(6(n^2-4)) in `015/spin1.py` was not checked against the paper by me; N22/N23 values depend on it, positivity of the bound does not.

## Wording
No sentence claims the Millennium problem solved or refuted, or a gap proved in the continuum. "Candidate C" and "draft" are used throughout. The one overreach is H1 ("blocked ... shown").

## Coverage table (prompt concepts and work folders)

| concept (OVERSEER_PROMPT) | in claim | status in claim | note |
|---|---|---|---|
| Schrodinger bound states / variational | rows 1, 6, 8 | examined | |
| min-max | row 1, Lean | examined | easy half of Courant-Fischer formalized |
| Temple | row 6 | examined | |
| Weinstein | "not examined" | listed | ok |
| Hardy inequalities | "not examined" | listed | ok |
| Bakry-Emery | rows 2, 7 | examined | |
| Payne-Weinberger / Gribov | row 3 | examined | M1, H3 |
| Bloch / periodic potential / Bragg | row 4 | examined | Bragg named only |
| tight-binding | row 4 (theta wells) | partial | KS-as-tight-binding not named (L7) |
| Mott, BCS, Peierls, excitonic | row 5 | examined | M4, L1 |
| Haldane gap / spin-1 chain | row 7 | examined | |
| AKLT | rows 2, 7 | examined | L8 |
| O(3) sigma model | row 7 | examined (cited TBA) | |
| Lieb-Schultz-Mattis | "not examined" | listed | ok |
| martingale, Knabe | row 2 | examined (007, 015) | |
| HF, CI, CC, DMRG/TN | row 6 | examined | PEPS listed not examined |
| Feynman / Karabali-Kim-Nair 2+1D | absent | **missing** | H3 |
| Gribov-Zwanziger gluon mass | row 5 (Gribov propagator) | examined | |
| Osterwalder-Seiler strong coupling | literature section | cited | |
| Moncrief/Mondal log-concavity | row 3 (a) | missing from "open" | H3 |
| Luscher / van Baal small volume | "not examined" | listed | ok |
| SU(3) KS, 3+1D | "not examined" | listed | ok |

Work folders: every research folder 003, 005-016, 018-020, 022-024, 026 maps to a row; 001 (literature), 002/004 (data), 017/021 (Lean setup/lemmas), 025 (this draft), 027 (full reproduce, running) are support. 011 is used (C8a/b/c) but never cited by folder number in `claim.md`; add `011` to the row 2 references.

## Number trace (seed 28, 13 rows)

| label | value in claim | script | committed output checked | found | recomputed here |
|---|---|---|---|---|---|
| N10 | 3.8715812 | 003 run_ed.py | 003/results.csv | yes | yes: 3.8715812223 (own character-basis ED, j <= 40) |
| N11 | 3.4005352 | 003 run_ed.py | 003/results.csv | yes | |
| N14 | 0.7694 | 016 check_c9.py | 016/run_c9.log (0.769384) | yes | printed value rounds a lower bound up (M3) |
| N15 | -7.861 | 012 analyze.py | 012/results.csv | yes | |
| N17 | -8.0000 | 010 hff_checks.py | 010/hff_checks.log | yes | |
| N23 | 0.36435607; 0.3026 | 015 spin1.py | 015/spin1_aklt_knabe.csv | yes | |
| N26a | 72.816 | 009 skeptic_checks.py | 009/skeptic_summary.txt | yes | |
| N26b | 150.411 | 009 skeptic_checks.py | 009/skeptic_summary.txt | yes | |
| N27 | -3.014e-2 | 009 skeptic_checks.py | 009/skeptic_summary.txt | yes | |
| N41b | +0.076 ... +3.8e-5 | 023 skeptic_c13.py | 023/summary.txt | yes | |
| N51 | 3.1323, 0.7187, -3.3044 | 019 fit.py | 019/run.log (map says fit.log, the sandbox name) | yes | |
| N52a | 94.0, 0/11 | 019 fit.py | 019/run.log | yes | |
| N56a | 1/9, chi2/n 104 (5.44) | 022 score.py | 022/score.log | yes | |

Extra recomputations (independent code, 1 thread):
- N6: plaquette gap g^2 = 1 = 3.2690991855 (claim 3.2690992; exact interval 3.26909918552).
- N9: series at g^2 = 10 gives 15.0037327365 vs ED 15.0037327366.
- N22: AKLT open chain n = 10 by sparse Lanczos (P_2 projectors): eps_10 = 0.3698444724; Knabe 0.291075, GM 0.301067 (claim 0.36984447, 0.3011).
- N40: Gribov Schwinger function at gamma t = 1.5: closed form e^{-a}(cos a - sin a)/(2 sqrt 2), a = t/sqrt 2, gives -4.70521e-2, quadrature agrees (claim -4.705e-2).
- N48: e^{-2 pi/0.35}/0.35 = 4.566e-8 (claim 4.57e-8).
All agree.

## Lean
`lean/MassGap.lean` imports Basic, BlochObstruction, SchurBound, StrongCoupling, MinMax. No `sorry`, `admit` or `axiom` in any file. All `.olean` files are newer than their sources (not rebuilt). Statement-level match with `claim.md`: C2 (first-order 3x/2 and second-order 56/(15 x^3), x = g^2, matrix elements as hypotheses-by-definition: matches "arithmetic only"); C4 (spectrum of a block sum = union, `band_data_do_not_bound_in_sector_gap` with 2x2 blocks, `overlapping_bands_example`: matches "finite-dimensional, two fibers"); C5 (`schur_lower_bound`, `feshbach_form_bound`, `schur_count`, `schur_count_sorted`: matches "k-eigenvalue count"); MinMax (easy half of Courant-Fischer: matches). No overstatement found; only CLAIMS.md is stale (L6).

## Open issues
- 027's number check needs the repair in M5 before it can pass.
- L8 (GM prefactor) not checked against the paper.

## Files
- `quick.log`: the in-repo `--quick` attempt that hit the concurrent full run (PermissionError).
