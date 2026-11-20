"""Check that flux (0,1) has the same spectrum as flux (1,0) (sigma_d maps one onto the other and
commutes with H, validate.py), at j_max = 3/2, per momentum block with k -> (ky, kx).
Usage (repo root): python work/20261120-012-torus-sectors/check_01.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
ART = ROOT / "artifacts" / "012-torus-sectors"


def main():
    b = {}
    for f in ("10", "01"):
        with open(ART / f"blocks_f{f}_tj3.pkl", "rb") as fh:
            b[f] = pickle.load(fh)["blocks"]
    worst = 0.0
    for k in b["10"]:
        for g2 in (1.0, 3.0):
            sp = []
            for f, kk in (("10", k), ("01", k[::-1])):
                reps, cas, W, S = b[f][kk]
                H = np.diag(0.5 * g2 * cas) - (2.0 / g2) * W.toarray()
                sp.append(np.linalg.eigvalsh(H)[:20])
            worst = max(worst, np.abs(sp[0] - sp[1]).max())
    print(f"(1,0) block k vs (0,1) block (ky,kx), lowest 20 levels, g2 = 1, 3: max diff {worst:.2e}")


if __name__ == "__main__":
    main()
