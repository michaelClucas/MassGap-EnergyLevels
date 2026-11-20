"""Scan of the k=2 (C0 = 0) relativistic model: for each m, the kappa that fits m(0++), and the resulting m(2++).
Shows whether a C0 = 0 solution with m >~ sqrt(sigma) exists (the least-squares fit ran to the m = 0.3 bound)."""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as M  # noqa: E402

LAT = M.load("glueball_su3_4d_athenodorou_teper_2020.csv", 3, 4, "AT2020")
t0, t2 = LAT[("0++", 0)][0], LAT[("2++", 0)][0]


def lev(m, k, jpc):
    return M.by_channel(M.spectrum4d((m, k, 0.0), Lmax=0, three=False))[jpc][0][0]


print("lattice 0++, 2++:", t0, t2)
print("   m   kappa(0++ fit)   m(2++)   m(2++)-m(0++)")
for m in [0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5]:
    ks = np.linspace(0.0, 4.0, 21)
    f = [lev(m, k, "0++") - t0 for k in ks]
    sol = None
    for i in range(len(ks) - 1):
        if f[i] * f[i + 1] < 0:
            sol = brentq(lambda k: lev(m, k, "0++") - t0, ks[i], ks[i + 1], xtol=1e-4)
            break
    if sol is None:
        print(f"{m:5.2f}  no kappa in [0,4] (0++ range {min(f) + t0:.2f}..{max(f) + t0:.2f})")
    else:
        e2 = lev(m, sol, "2++")
        print(f"{m:5.2f}  {sol:8.4f}        {e2:7.3f}   {e2 - t0:7.3f}")
