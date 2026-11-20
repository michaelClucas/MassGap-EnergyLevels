"""Certified two-sided gap bounds for the O(3) rotor chain, H = (g/2) sum L^2 - (1/g) sum n.n.
Run from repo root: python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/run_rotor.py
Writes rotor_certified.csv (one row per L, bc, lmax, g)."""
import csv, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rotor

HERE = os.path.dirname(os.path.abspath(__file__))
GS = [5.0, 3.0, 2.0, 1.5, 1.0, 0.7, 0.5, 0.35, 0.25]
CONF = [(2, "open", 4), (2, "open", 8),
        (3, "open", 5), (3, "per", 3), (3, "per", 5),
        (4, "open", 1), (4, "open", 2), (4, "open", 3), (4, "per", 1), (4, "per", 2), (4, "per", 3),
        (5, "open", 2), (5, "per", 2), (6, "open", 2), (6, "per", 2),
        (6, "per", 1), (8, "open", 1), (8, "per", 1), (10, "open", 1), (10, "per", 1)]
out = os.path.join(HERE, "rotor_certified.csv")
if len(sys.argv) > 1 and sys.argv[1] == "large":  # second process, coarse g grid, large L
    CONF = [(6, "per", 2), (8, "per", 1), (10, "per", 1), (10, "open", 1)]
    GS = [2.0, 1.0, 0.7, 0.5]
    out = os.path.join(HERE, "rotor_certified_large.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["L", "bc", "lmax", "g", "dim", "floor", "E0_lo", "E0_hi", "E1_lo", "E1_hi", "gap_lo", "gap_hi", "secs"])
    for L, bc, lmax in CONF:
        for g in GS:
            t = time.time()
            r = rotor.certify(L, bc, g / 2, 1 / g, lmax, depth=1 if L * lmax <= 12 or lmax == 1 and L <= 6 else 0)
            glo = r["L1"] - r["U0"] if np.isfinite(r["L1"]) else -np.inf
            row = [L, bc, lmax, g, r["dim"], f"{r['floor']:.6g}", f"{r['L0']:.10g}", f"{r['U0']:.10g}",
                   f"{r['L1']:.10g}", f"{r['U1']:.10g}", f"{glo:.10g}", f"{r['U1'] - r['L0'] if np.isfinite(r['L0']) else np.inf:.10g}",
                   f"{time.time() - t:.1f}"]
            w.writerow(row); f.flush()
            print(*row, flush=True)
