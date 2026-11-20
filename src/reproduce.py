"""Regenerate the program's numbers. Run from the repo root: python -m src.reproduce [--quick]

Current contents:
  work/20261120-003-su2-ks-ed: SU(2) Kogut-Susskind exact diagonalization (validation + gap scan).
  work/20261120-006-certified-gap-bounds: certified two-sided bounds on E0, E1, gap of the
    untruncated KS Hamiltonian (full mode includes the 2x2 torus, ~10+ min).
"""
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TASKS = [
    "work/20261120-003-su2-ks-ed/validate.py",
    "work/20261120-003-su2-ks-ed/run_ed.py",
    "work/20261120-006-certified-gap-bounds/certify_gaps.py",
]


def main():
    for script in TASKS:
        print(f"=== {script}", flush=True)
        argv = sys.argv
        sys.argv = [script] + [a for a in argv[1:] if a == "--quick"]
        runpy.run_path(str(ROOT / script), run_name="__main__")
        sys.argv = argv


if __name__ == "__main__":
    main()
