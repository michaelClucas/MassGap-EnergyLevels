# 20261120-025: Result-C draft and reproduce pipeline

## Question
Prepare what a Result C needs so the overseer can decide once the open checks are in: a claim file covering every technique family examined, and a `python -m src.reproduce` that regenerates every number in it.

## Method
- Read `CLAIMS.md`, `STATE.md` and every work README (001-021), then folded in 022 (post-2008 bound-state models), 023 (C13 skeptic), 024 (inertia counts) and 026 (Lean min-max) as relayed by the overseer. Each number in `claim.md` was taken from a work-folder output file, not from a README.
- `src/reproduce.py` rewritten: tasks run on copies of the work-folder scripts under `artifacts/<folder>/` (gitignored, same depth as `work/`, so `parents[2]` is still the repo root); committed outputs are never overwritten. After the runs it checks every row of `number_map.csv` (literal match string in the regenerated output, or in the committed file for rows marked `recorded`). Existing tasks (003 validate/run_ed, 006 certify_gaps) are kept.
- `number_map.csv`: label (as in `claim.md`), value as quoted, script, output file, match string, mode (`quick`, `full`, `recorded`), note.

## Result
- `claim.md`: candidate Result C, technique table (8 families), open items, not-examined list, literature summary, weaknesses, reproduce commands.
- `--quick`: 72 s, 15 scripts, number check 29 ok / 0 failed (`run.log`).
- Full mode: not run (estimate above 45 min). Every full-mode match string was checked against the committed outputs of the original runs (all 40 found), so the map is consistent; regeneration in full mode is untested end to end.
- Lean is not built by reproduce.py: `cd lean && lake build`.

## Runtime
- Quick: 72 s (1 thread, pinned).
- Full, estimated at 1 thread while other programs share CPUs 9-11: about 2-2.5 h. Breakdown: 019 `fit.py` measured 4 of 7 fits in 42 min (013 README quoted ~10 min on a quieter machine), so ~75 min here; 006 `certify_gaps.py` full ~10+ min; 023 `skeptic_c13.py` ~10 min; 018 ~7.5 min; 009 ~7 min; 003 ~4 min; 020 a few min; 016 check_c8 ~2 min; 010 plaq_exact, 015 spin1/sc_series/bakry_emery, 019 limits not timed (expected minutes). Above the 45-min threshold, so the full run was not started; the partial timing run was stopped.

## Verdict
Draft ready for the overseer. Every number in `claim.md` maps to a script and output; rows marked `recorded` (2x2 torus j_max 2 builds, rotor certification run) are checked against committed outputs because the producing runs need ~5 GB or > 30 min.

## Open issues
- Transcription error found: `CLAIMS.md` C5 / 006 README give the ladder gap at g^2 = 1 as 3.0620125; `006/results.csv` gives 3.0620133 (Delta in [3.0620132814, 3.0620132822]).
- `recorded` rows are not regenerated (torus builds: RAM; `015/run_rotor.py`: > 30 min and never completed in 015).
- 008 `gribov_scope.py` was dropped from reproduce: its numbers are superseded by 009 and not used in the claim.
