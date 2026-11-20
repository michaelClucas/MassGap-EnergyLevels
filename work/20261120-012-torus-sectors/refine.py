"""Refine certified lower bounds on the zero-flux k = 0 block (2x2 torus, SU(2) KS, 2+1D, finite
lattice): when the standard choice L = f(U) is poor (U close to the floor d), bisect for the largest
E < min(d, U) with lambda_k(A - X/(d-E)) >= E. Same Lemma as analyze.py; only the choice of E differs.
P = spins <= 3/2, exact shell X from spin-2 states, split floor with c_Q = 11.
Usage (repo root): python work/20261120-012-torus-sectors/refine.py g2 [g2 ...]
Appends to refine.csv.
"""
import csv
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import certify as C  # noqa: E402
from src import torus_sym as T  # noqa: E402

ART = ROOT / "artifacts" / "012-torus-sectors"
OUT = Path(__file__).resolve().parent


def main():
    with open(ART / "blocks_f00_tj4.pkl", "rb") as fh:
        reps, cas, W, S = pickle.load(fh)["blocks"][(0, 0)]
    maxsp = np.array([max(r[0]) for r in reps])
    iP, ish = np.where(maxsp <= 3)[0], np.where(maxsp == 4)[0]
    new = not (OUT / "refine.csv").exists()
    with open(OUT / "refine.csv", "a", newline="") as fh:
        wr = csv.writer(fh)
        if new:
            wr.writerow(["sector", "momentum", "g2", "j_max", "level", "lo", "hi_P", "d_floor", "method"])
        for g2 in [float(a) for a in sys.argv[1:]]:
            A, X = C.block_ops(cas, W, g2, iP, 4, shell=ish)
            d, t = T.split_floor(g2, 11.0, ts=np.linspace(0, 1, 201))
            Ad, Xd = A.toarray(), X.toarray()
            hiP = C.upper_bounds(A, 2) + C._eps_margin(A)
            for k in (0, 1):
                def f(E):
                    return C.lower_bound_k(Ad - Xd / (d - E), k)
                b = min(d, hiP[k]) - 1e-9
                a = f(b) if hiP[k] < d else min(-1e3, d - 1e3)
                for _ in range(16):
                    m = 0.5 * (a + b)
                    if f(m) >= m:
                        a = m
                    else:
                        b = m
                wr.writerow(["(0,0)", "(0,0)", g2, 1.5, k, f"{a:.10g}", f"{hiP[k]:.10g}", f"{d:.10g}",
                             f"schur-exact-shell bisect; split floor t={t:.3f}"])
                fh.flush()
                print(f"g2={g2} (0,0) k=0 E{k} >= {a:.6f} (U_P={hiP[k]:.6f}, d={d:.4f})", flush=True)


if __name__ == "__main__":
    main()
