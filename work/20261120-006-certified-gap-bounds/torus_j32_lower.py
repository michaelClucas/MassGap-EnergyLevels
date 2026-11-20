"""Certified lower bounds for the periodic 2x2 zero-flux sector at j_max = 3/2 (dim 9418) via the
sparse LDL^T inertia route of src/certify.py (crude shell X, crude floor d). Slow: ~2-5 min per
factorization, two per coupling. Not in the quick reproduce path.

Usage (repo root): python work/20261120-006-certified-gap-bounds/torus_j32_lower.py [g2 ...]
Appends rows to results_periodic2x2_j32.csv; certify_gaps.py merges it into results.csv.
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402
from src import certify as C  # noqa: E402

OUT = Path(__file__).resolve().parent
FIELDS = ["lattice", "sector", "g2", "j_max", "dim", "E0_lo", "E0_hi", "E1_lo", "E1_hi",
          "Delta_lo", "Delta_hi", "d_floor", "d_source", "c_Q", "method"]


def main():
    g2s = [float(a) for a in sys.argv[1:]] or [3.0]
    lat = ks.periodic(2, 2)
    basis, cas, W = C.build_sector(lat, 3)
    maxsp = np.array([max(st[0]) for st in basis])
    iP = np.arange(len(basis))
    c_Q, _ = C.casimir_floor(lat, 3)
    nmax = C.plaq_count_at_max(lat, basis, maxsp <= 3, 3)
    f = OUT / "results_periodic2x2_j32.csv"
    new = not f.exists()
    with open(f, "a", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS)
        if new:
            w.writeheader()
        for g2 in g2s:
            t = time.time()
            r = C.certify(cas, W, g2, len(lat.plaqs), iP, c_Q, at_max_P=nmax)
            lo, hi = r["lo"], r["hi"]
            row = dict(lattice="periodic2x2", sector="gauge-invariant zero Z2-flux; all momenta",
                       g2=g2, j_max=1.5, dim=r["dim"], E0_lo=lo[0], E0_hi=hi[0], E1_lo=lo[1],
                       E1_hi=hi[1], Delta_lo=lo[1] - hi[0], Delta_hi=hi[1] - lo[0],
                       d_floor=r["d"], d_source=r["dsrc"], c_Q=c_Q,
                       method=r["method"] + "+sparse-ldl")
            w.writerow(row)
            fh.flush()
            print(f"g2={g2}: Delta in [{row['Delta_lo']:.6g}, {row['Delta_hi']:.6g}]  "
                  f"E0 [{lo[0]:.8g},{hi[0]:.8g}] E1 [{lo[1]:.8g},{hi[1]:.8g}] d={r['d']:.4g} "
                  f"{time.time()-t:.0f}s", flush=True)


if __name__ == "__main__":
    main()
