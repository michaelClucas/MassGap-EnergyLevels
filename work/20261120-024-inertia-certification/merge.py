"""Merge inertia_c5.csv, inertia_c9.csv, inertia_c11.csv, inertia_c11_aklt.csv into inertia.csv and
print verdict counts per claim. Usage: python work/20261120-024-inertia-certification/merge.py"""
import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = []
for f in ("inertia_c5.csv", "inertia_c9.csv", "inertia_c11.csv", "inertia_c11_aklt.csv"):
    if (HERE / f).exists():
        rows += list(csv.DictReader(open(HERE / f)))
with open(HERE / "inertia.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
cnt = collections.Counter((r["claim"], r["system"].split(" L=")[0].split(", flux")[0].split(" n=")[0],
                           "weakened" if r["verdict"].startswith("verified at") else r["verdict"]) for r in rows)
for k, v in sorted(cnt.items()):
    print(*k, v, sep=" | ")
wk = [float(r["verdict"].split("- ")[1]) for r in rows if r["verdict"].startswith("verified at")]
if wk:
    print("max weakening", max(wk))
