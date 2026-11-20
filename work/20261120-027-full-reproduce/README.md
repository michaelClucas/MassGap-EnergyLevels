# 027: full reproduce run

**Question.** Does `python -m src.reproduce` (full mode) regenerate every number in the Result-C claim (`work/20261120-025-resultC-draft/number_map.csv`)?

**Method.** Full mode, 1 thread, pinned via `tools/pin.py` (`run.log`, 2100 s of script runs). During the run a skeptic's concurrent `--quick` run deleted the artifact folders of 002, 004, 005, 007, 008 (quick and full shared `artifacts/`; a lock now prevents this). Those five were rerun with `--only` (`repair.log`), then all numbers were rechecked with `--check-only` (`check.log`).

**Result.** All scripts exited 0. First check: 61/69 ok; 7 failures were the deleted files; 1 (N8, ladder Delta_lo at j_max 6) was a 16-digit match string against a regenerated value differing by 2e-15 (BLAS rounding: 3.0620132814171055 vs recorded 3.062013281417123); the match was truncated to 14 digits. After the repair and recheck: **69/69 ok, 0 failed** (`check.log`).

**Verdict.** Reproducibility stop condition met.

**Open issues.** Seven numbers are mode `recorded` (checked against committed outputs, not regenerated): 2x2 torus j_max = 2 builds (~5 GB each) and the long rotor certification run.
