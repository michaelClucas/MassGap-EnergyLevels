"""Grid convergence check: low states at given params on the default grid vs a larger, finer grid.
Usage: conv.py m kappa C0 m3 C03
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as M  # noqa: E402

m, kappa, C0, m3, C03 = map(float, sys.argv[1:6])
g2N = 1.0 / M.load_scalar("glueball_sun_3d_athenodorou_teper_2016.csv", 3, "sqrt_sigma/g2N")[0]
keys = [("0++", 0), ("0++", 1), ("2++", 0), ("0-+", 0), ("2-+", 0), ("3++", 0), ("1--", 0)]
a = M.by_channel(M.spectrum4d((m, kappa, C0)))
keys3 = [("0++", 0), ("0++", 1), ("2++", 0), ("0--", 0)]
b3 = M.by_channel(M.spectrum3d((m3, C03), 3, g2N))
M.GRID4 = M.Grid(R=20.0, n=1100)
M.GRID4_3 = M.Grid(R=20.0, n=1100)
M.GRID3 = M.Grid(R=22.0, n=1400)
M.GRID3_3 = M.Grid(R=22.0, n=1100)
b = M.by_channel(M.spectrum4d((m, kappa, C0)))
c3 = M.by_channel(M.spectrum3d((m3, C03), 3, g2N))
for k in keys:
    print("4D", k, round(a[k[0]][k[1]][0], 4), round(b[k[0]][k[1]][0], 4))
for k in keys3:
    print("3D", k, round(b3[k[0]][k[1]][0], 4), round(c3[k[0]][k[1]][0], 4))
