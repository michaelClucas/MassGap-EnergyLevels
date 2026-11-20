"""C11 (work 015/018): inertia check of the AKLT open-chain gap eps_n used in the Gosset-Mozgunov /
Knabe bounds. 1D spin-1 chain, open, n sites, H = sum_i P_2(i,i+1) (projectors), S_z = 0 sector
(every multiplet of integer total spin has an S_z = 0 member, so the S_z = 0 spectrum contains every
eigenvalue of H). H >= 0 and its kernel on the open chain is 4-dimensional (AKLT edge states,
cited), 2 of which lie in S_z = 0. The certificate needs: every eigenvalue of H other than the
kernel is >= eps_claimed, i.e. #eig(H|Sz=0 < eps_claimed) <= 2. Also reported: #eig < 1e-6
(should be 2: the kernel) as a consistency check.
Usage (repo root): python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c11_aklt.py <n> [<n> ...]
"""
import sys
import time

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W015 = ROOT / "work" / "20261120-015-line6-rotor-chain"
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE)); sys.path.insert(0, str(W015))
import spin1  # noqa: E402
from src import inertia as I  # noqa: E402
import common  # noqa: E402
import memguard  # noqa: E402

CLAIMED = {10: 0.36984447, 12: 0.36435607}  # 015 spin1_aklt_knabe.csv / 018 run.log (8 digits)
CLAIMS_MD = {12: 0.364356}  # value quoted in CLAIMS.md C11


def main():
    memguard.start(4600)
    rows = []
    for n in [int(a) for a in sys.argv[1:]]:
        t = time.time()
        H = spin1.chain(spin1.P2, n, False).tocsc()
        for s, what in ((CLAIMED[n], "eps_n (8 digits, 015)"),) + (((CLAIMS_MD[n], "eps_n (CLAIMS.md)"),) if n in CLAIMS_MD else ()) + ((1e-6, "kernel count"),):
            memguard.state["label"] = f"AKLT n={n} s={s}"
            res = I.at_most(H, s, 2)
            rows.append(dict(claim="C11", work="015/018", system=f"AKLT open chain n={n}",
                             setting="1D spin-1, open chain, S_z=0 sector, h = P_2", g="", j_or_l_max="",
                             level_k=2, quantity=what, claimed_lower=s, dim=H.shape[0],
                             lanczos_used="yes (eigsh)" if n > 10 else "no (dense)", matrix="H",
                             allowed=2, **common.res_fields(res),
                             note="lower bound on eps_n: all non-kernel eigenvalues >= s" if s > 1e-3 else
                             "consistency: count below 1e-6 should equal kernel dim 2"))
            print(f"n={n} dim {H.shape[0]} s={s}: count in [{res['n_lo']},{res['n_hi']}] -> {res['verdict']} "
                  f"(r={np.nanmax([res['r_lo'], res['r_hi']]):.1e}, {res['method']}) {time.time()-t:.0f}s", flush=True)
    common.write("inertia_c11_aklt.csv", rows)


if __name__ == "__main__":
    main()
