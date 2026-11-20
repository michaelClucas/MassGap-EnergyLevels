"""Regenerate the program's numbers. Run from the repo root: python -m src.reproduce [--quick]

Current contents:
  work/20261120-003-su2-ks-ed: SU(2) Kogut-Susskind exact diagonalization (validation + gap scan).
"""
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TASKS = [
    "work/20261120-003-su2-ks-ed/validate.py",
    "work/20261120-003-su2-ks-ed/run_ed.py",
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
