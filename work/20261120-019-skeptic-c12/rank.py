"""Skeptic 019: rank models from comparison_verified.csv by several metrics, incl. on the common low-lying set.

Run after compare.py (same folder). Prints a table; no lattice numbers here (they come via comparison_verified.csv,
which compare.py builds from data/).
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "comparison_verified.csv", newline="")))
COMMON = {("2++", "0"), ("0-+", "0"), ("0++", "1")}
models = {}
for r in rows:
    models.setdefault(r["model"], []).append(r)
print(f"{'model':40s} {'n':>3s} {'in1sig':>6s} {'chi2/n lat':>10s} {'chi2/n 5%':>9s} | common(2++,0-+,0++*): "
      f"{'n':>2s} {'in1':>3s} {'chi2/n lat':>10s} {'5%':>6s}")
for name, rs in models.items():
    pl = [float(r["pull_lattice"]) for r in rs]
    pb = [float(r["pull_budget5pct"]) for r in rs]
    cs = [r for r in rs if (r["JPC"], r["level"]) in COMMON]
    cpl = [float(r["pull_lattice"]) for r in cs]
    cpb = [float(r["pull_budget5pct"]) for r in cs]
    n, nc = len(pl), len(cpl)
    print(f"{name[:40]:40s} {n:3d} {sum(abs(x) <= 1 for x in pl):6d} {sum(x * x for x in pl) / n:10.1f} "
          f"{sum(x * x for x in pb) / n:9.2f} | {nc:2d} {sum(abs(x) <= 1 for x in cpl):3d} "
          f"{(sum(x * x for x in cpl) / nc if nc else float('nan')):10.1f} "
          f"{(sum(x * x for x in cpb) / nc if nc else float('nan')):6.2f}")
