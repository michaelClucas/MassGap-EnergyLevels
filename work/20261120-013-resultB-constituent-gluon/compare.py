"""Compare models on the same 4D SU(3) data (AT2020) as ratios m/m(0++), which removes each model's scale.

Lattice ratio errors: m/sqrt(sigma) errors in quadrature (correlations ignored, so slightly too large).
Reads competitors.csv (published model numbers, cited per row) and fits.csv (this work). Writes comparison.csv.
"""
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as M  # noqa: E402

HERE = Path(__file__).resolve().parent
LAT = M.load("glueball_su3_4d_athenodorou_teper_2020.csv", 3, 4, "AT2020")
MODEL_ERR = 0.05
STATES = [("2++", 0), ("0-+", 0), ("0++", 1), ("2-+", 0), ("1+-", 0), ("2++", 1), ("3++", 0), ("1--", 0),
          ("2--", 0), ("0-+", 1), ("2-+", 1), ("3+-", 0), ("1-+", 0), ("4++", 0)]

m0, e0 = LAT[("0++", 0)]
lat_ratio = {}
for s in STATES:
    v, e = LAT[s]
    lat_ratio[s] = (v / m0, (v / m0) * math.hypot(e / v, e0 / m0))

models = {}
with open(HERE / "competitors.csv", newline="") as f:
    for r in csv.DictReader(f):
        models.setdefault(r["model"], {})[(r["JPC"], int(r["level"]))] = float(r["value"])
with open(HERE / "fits.csv", newline="") as f:
    for r in csv.DictReader(f):
        if r["setting"] == "4D SU(3) AT2020" and r["model"] != "":
            models.setdefault("this work: " + r["variant"], {})[(r["JPC"], int(r["level"]))] = float(r["model"])

out = []
print(f"{'model':42s} " + " ".join(f"{j}{'*' * l:<5s}" for j, l in STATES))
print(f"{'lattice AT2020':42s} " + " ".join(f"{lat_ratio[s][0]:6.3f}" for s in STATES))
for name, d in models.items():
    base = d.get(("0++", 0))
    if base is None:
        continue
    cells, n, nin, chi_l, chi_b = [], 0, 0, 0.0, 0.0
    for s in STATES:
        if s in d:
            r = d[s] / base
            lv, le = lat_ratio[s]
            pl = (r - lv) / le
            pb = (r - lv) / math.hypot(le, MODEL_ERR * lv)
            n += 1
            nin += abs(pl) <= 1
            chi_l += pl ** 2
            chi_b += pb ** 2
            cells.append(f"{r:6.3f}")
            out.append(dict(model=name, JPC=s[0], level=s[1], model_ratio=round(r, 4), lattice_ratio=round(lv, 4),
                            lattice_err=round(le, 4), pull_lattice=round(pl, 2), pull_budget5pct=round(pb, 2)))
        else:
            cells.append("   -  ")
    print(f"{name[:42]:42s} " + " ".join(cells) + f"   n={n} in1sig={nin} chi2/n lat={chi_l / max(n, 1):.0f} "
          f"5%={chi_b / max(n, 1):.2f}")
with open(HERE / "comparison.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
