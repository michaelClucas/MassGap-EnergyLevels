"""Regenerate the program's numbers and check them against the Result-C draft.

Run from the repo root, pinned to this program's CPUs, one thread:

    OMP_NUM_THREADS=1 python tools/pin.py .venv/Scripts/python.exe -m src.reproduce --quick
    OMP_NUM_THREADS=1 python tools/pin.py .venv/Scripts/python.exe -m src.reproduce          (full)

Options: --quick (fast subset, ~2-3 min), --only <substring> (run matching tasks only),
--check-only (skip running, only check numbers), --list (print the task table).

How it works
- Every task is a script in a work folder. The script is copied (all *.py files of the folder,
  plus listed input files) to artifacts/<folder>/ and run from there, with cwd = repo root.
  artifacts/ is gitignored and sits at the same depth as work/, so each script's
  Path(__file__).parents[2] is still the repo root and data/ and src/ resolve as usual.
  Committed outputs in work/ are never overwritten. stdout+stderr go to
  artifacts/<folder>/<script>[_args].log.
- After the runs, every row of work/20261120-025-resultC-draft/number_map.csv whose mode
  matches is checked: its `match` string must occur in the regenerated output file under
  artifacts/<folder>/ (modes quick/both/full), or in the committed file (mode "recorded": the
  producing run is too large to repeat here, see the row's note). The run ends with a
  pass/fail count and exits nonzero on any failure.
- Only one run at a time: main() creates artifacts/.reproduce.lock exclusively (PID inside) and
  removes it at exit; a second run exits with code 2. A stale lock from a killed run must be
  deleted by hand after checking that the PID is gone.
- Full mode runs everything in TASKS; quick mode runs the tasks marked "quick". Full-mode
  runtime is stated in work/20261120-025-resultC-draft/README.md.

Lean proofs are not built here. Build them separately with:

    cd lean && lake build

(Lean v4.36.0 + mathlib v4.36.0; see work/20261120-021-lean-lemmas/README.md.)
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SANDBOX = ROOT / "artifacts"
NUMBER_MAP = ROOT / "work" / "20261120-025-resultC-draft" / "number_map.csv"

# (folder, script, args, mode, extra input files copied next to the script)
# mode "quick": runs in quick and full mode; "full": full mode only.
# An input named "quick:<file>" is copied in quick mode only.
# Args in a 2-tuple (quick_args, full_args) differ between modes.
TASKS = [
    ("20261120-002-glueball-data", "ratios.py", [], "quick", []),
    ("20261120-004-data-verify", "verify_ratios.py", [], "quick", []),
    ("20261120-005-skeptic-ks-ed", "plaq_independent.py", [], "quick", []),
    ("20261120-007-frustration-free-gauge", "checks.py", [], "quick", []),
    ("20261120-008-theta-obstruction-gribov-scope", "theta_bands.py", [], "quick", []),
    ("20261120-010-skeptic-certify-hff", "hff_checks.py", [], "quick", []),
    ("20261120-014-line5-interaction-gaps", "interaction_gaps.py", [], "quick", []),
    ("20261120-016-skeptic-c8-c9", "check_c9.py", [], "quick", []),
    # 016 check_c9 reads ../20261120-012-torus-sectors/{results,refine}.csv (recorded block bounds)
    ("20261120-012-torus-sectors", None, [], "quick", ["results.csv", "refine.csv"]),
    ("20261120-019-skeptic-c12", "fit.py", [], "full", ["competitors_verified.csv"]),
    # quick mode compares the committed fit parameters (fits.csv); full mode refits first
    ("20261120-019-skeptic-c12", "compare.py", [], "quick", ["competitors_verified.csv", "quick:fits.csv"]),
    ("20261120-019-skeptic-c12", "rank.py", [], "quick", []),
    # 022 imports model.py from ../20261120-019-skeptic-c12 (copied above)
    ("20261120-022-post2008-boundstate-glueballs", "score.py", [], "quick", ["models.csv"]),
    # pre-existing reproduce tasks (iteration 1-2); --quick passed through in quick mode
    ("20261120-003-su2-ks-ed", "validate.py", (["--quick"], []), "quick", []),
    ("20261120-003-su2-ks-ed", "run_ed.py", (["--quick"], []), "quick", []),
    ("20261120-006-certified-gap-bounds", "certify_gaps.py", (["--quick"], []), "quick", []),
    # full mode only
    ("20261120-010-skeptic-certify-hff", "plaq_exact.py", [], "full", []),
    ("20261120-016-skeptic-c8-c9", "check_c8.py", ["ABC"], "full", []),
    ("20261120-019-skeptic-c12", "limits.py", [], "full", []),
    ("20261120-015-line6-rotor-chain", "sc_series.py", [], "full", []),
    ("20261120-015-line6-rotor-chain", "bakry_emery.py", [], "full", []),
    ("20261120-015-line6-rotor-chain", "spin1.py", [], "full", []),
    ("20261120-020-quantum-chemistry-methods", "qc_methods.py", [], "full", []),
    ("20261120-023-skeptic-c13", "skeptic_c13.py", [], "full", []),
    ("20261120-009-skeptic-theta-gribov", "skeptic_checks.py", [], "full", []),
    ("20261120-018-skeptic-c10-c11", "skeptic_checks.py", [], "full", []),
]


def task_args(args, quick):
    if isinstance(args, tuple):
        return args[0] if quick else args[1]
    return args


def selected(quick, only):
    out = []
    for t in TASKS:
        folder, script, args, mode, inputs = t
        if quick and mode != "quick":
            continue
        if only and not any(o in f"{folder}/{script}" for o in only):
            continue
        out.append(t)
    return out


def prepare(tasks, quick):
    """Fresh sandbox copy of every folder used: all *.py (recursive) plus listed inputs."""
    inputs = {}
    for folder, _, _, _, files in tasks:
        keep = {f.split(":", 1)[1] if f.startswith("quick:") else f
                for f in files if quick or not f.startswith("quick:")}
        inputs.setdefault(folder, set()).update(keep)
    for folder, files in inputs.items():
        src, dst = ROOT / "work" / folder, SANDBOX / folder
        assert src.is_dir(), src
        if dst.exists():
            shutil.rmtree(dst)
        for p in src.rglob("*.py"):
            q = dst / p.relative_to(src)
            q.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, q)
        for f in files:
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src / f, dst / f)


def run(tasks, quick):
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env.setdefault(k, "1")
    env["PYTHONIOENCODING"] = "utf-8"
    failed, t_all = [], time.time()
    for folder, script, args, _, _ in tasks:
        if script is None:
            continue
        a = task_args(args, quick)
        log = SANDBOX / folder / (Path(script).stem + ("_" + "_".join(x.strip("-") for x in a) if a else "") + ".log")
        print(f"=== {folder}/{script} {' '.join(a)}", flush=True)
        t0 = time.time()
        with open(log, "w", encoding="utf-8") as fh:
            rc = subprocess.call([sys.executable, "-u", str(SANDBOX / folder / script), *a],
                                 cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        print(f"    exit {rc}, {dt:.0f} s, log {log.relative_to(ROOT).as_posix()}", flush=True)
        if rc != 0:
            failed.append(f"{folder}/{script}")
    print(f"=== runs done in {time.time() - t_all:.0f} s; failed: {failed or 'none'}", flush=True)
    return failed


def check(quick, only):
    modes = {"quick", "both", "recorded"} if quick else {"quick", "both", "full", "recorded"}
    if not NUMBER_MAP.exists():
        print(f"=== no {NUMBER_MAP.relative_to(ROOT).as_posix()}; nothing to check")
        return []
    rows = list(csv.DictReader(open(NUMBER_MAP, newline="", encoding="utf-8")))
    n_ok, bad, skipped = 0, [], 0
    for r in rows:
        if r["mode"] not in modes or (only and not any(o in r["output"] for o in only)):
            skipped += 1
            continue
        base = ROOT if r["mode"] == "recorded" else SANDBOX
        path = base / r["output"]
        text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else None
        if text is not None and r["match"] in text:
            n_ok += 1
        else:
            bad.append(f"{r['label']}: '{r['match']}' not in {path.relative_to(ROOT).as_posix()}"
                       + ("" if text is not None else " (file missing)"))
    print(f"=== number check ({'quick' if quick else 'full'}): {n_ok} ok, {len(bad)} failed, "
          f"{skipped} not in this mode", flush=True)
    for b in bad:
        print("    FAIL " + b, flush=True)
    return bad


LOCK = SANDBOX / ".reproduce.lock"


def acquire_lock():
    """Create artifacts/.reproduce.lock exclusively; return False if another run holds it.

    Quick and full mode share artifacts/ and prepare() deletes task folders there, so two
    overlapping runs corrupt each other. A lock left by a killed run is not removed
    automatically: the message names the PID so the user can check and delete it.
    """
    SANDBOX.mkdir(exist_ok=True)
    try:
        fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            info = LOCK.read_text(encoding="utf-8").strip()
        except OSError:
            info = "(unreadable)"
        print(f"=== another reproduce run holds {LOCK.relative_to(ROOT).as_posix()}: {info}\n"
              f"    Wait for it to finish. If no process with that PID is running (stale lock\n"
              f"    from a killed run), delete {LOCK.relative_to(ROOT).as_posix()} and rerun.",
              file=sys.stderr, flush=True)
        return False
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(f"pid {os.getpid()} started {time.strftime('%Y-%m-%d %H:%M:%S')} args {' '.join(sys.argv[1:])}\n")
    return True


def release_lock():
    try:
        LOCK.unlink()
    except FileNotFoundError:
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--list", action="store_true")
    o = ap.parse_args()
    tasks = selected(o.quick, o.only)
    if o.list:
        for folder, script, args, mode, inputs in TASKS:
            print(f"{mode:5s} {folder}/{script} {args} inputs={inputs}")
        return 0
    # --check-only also locks: it would read half-written outputs of a running run
    if not acquire_lock():
        return 2
    try:
        failed = []
        if not o.check_only:
            prepare(tasks, o.quick)
            failed = run(tasks, o.quick)
        bad = check(o.quick, o.only)
    finally:
        release_lock()
    return 1 if (failed or bad) else 0


if __name__ == "__main__":
    sys.exit(main())
