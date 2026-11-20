"""C5 (work 006): inertia check of every finite certified lower bound on the open 2-plaquette ladder
and the periodic 2x2 torus (zero flux, j_max <= 1), SU(2) KS, 2+1D, finite lattice, fixed g^2.
Every lower bound in src/certify.certify comes from lower_bound_k(M, k) (dense LAPACK for these
dims, n <= 622: no Lanczos was used). Here lower_bound_k is wrapped: each returned value v is
re-checked by a Bunch-Kaufman inertia count of M - vI ('#eig(M < v) <= k'). Calls are the nested
floor d (k = 0 on Q1) and E0, E1 (k = 0, 1 of A - X/(d - U_k)).
Usage (repo root): python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c5_ladder_torus.py
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from src import ks_su2 as ks  # noqa: E402
from src import certify as C  # noqa: E402
from src import inertia as I  # noqa: E402
import common  # noqa: E402

G2 = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
PLAN = {"ladder2": (ks.ladder2, 13, [2, 4, 8, 12]), "periodic2x2": (lambda: ks.periodic(2, 2), 3, [1, 2])}
calls = []
_orig = C.lower_bound_k


def wrapped(M, k):
    v = _orig(M, k)
    calls.append((M, k, v))
    return v


C.lower_bound_k = wrapped


def add(rows, name, tj, g2, k, lab, v, dim, matrix, res, note):
    rows.append(dict(claim="C5", work="006", system=name,
                     setting="SU(2) KS 2+1D finite lattice" + (" zero flux" if "periodic" in name else ""),
                     g=g2, j_or_l_max=tj / 2, level_k=k, quantity=lab, claimed_lower=v, dim=dim,
                     lanczos_used="no (dense LAPACK)", matrix=matrix, allowed=k, **common.res_fields(res), note=note))
    print(f"{name} j={tj/2} g2={g2} {lab}: lo={v:.10g} dim={dim} count in [{res['n_lo']},{res['n_hi']}] "
          f"allowed {k} -> {res['verdict']} (r={np.nanmax([res['r_lo'], res['r_hi']]):.1e})", flush=True)


def main():
    rows = []
    for name, (ctor, top, tjs) in PLAN.items():
        lat = ctor()
        basis, cas, W = C.build_sector(lat, top)
        maxsp = np.array([max(st[0]) for st in basis])
        c_top, _ = C.casimir_floor(lat, top)
        for tj in tjs:
            iP = np.where(maxsp <= tj)[0]
            c_Q, _ = C.casimir_floor(lat, tj)
            kw = dict(shell=np.where(maxsp == tj + 1)[0])
            iQ1 = np.where(maxsp > tj)[0]
            kw.update(iQ1=iQ1, c_Q_top=c_top,
                      at_top_Q1=C.plaq_count_at_max(lat, [basis[i] for i in iQ1], [True] * len(iQ1), top))
            for g2 in G2:
                calls.clear()
                t = time.time()
                r = C.certify(cas, W, g2, len(lat.plaqs), iP, c_Q, **kw)
                # 1. nested floor: re-check the k = 0 call on Q1; if only a weakened value verifies, use it
                d = r["d"]
                fl_calls = [c for c in calls if c[0].shape[0] != len(iP)]
                for M, k, v in fl_calls:
                    if r["dsrc"] != "nested" or not np.isfinite(v):
                        continue
                    res = I.at_most(M, v, 0)
                    if res["verdict"].startswith("verified at claimed - "):
                        d = max(v - float(res["verdict"].split("- ")[1]), 0.5 * g2 * c_Q - 4.0 * len(lat.plaqs) / g2)
                    elif res["verdict"] != "verified":
                        d = 0.5 * g2 * c_Q - 4.0 * len(lat.plaqs) / g2  # fall back to the crude floor
                    add(rows, name, tj, g2, 0, "floor d (nested, Q1)", v, M.shape[0], "A1 - X1/(d2 - U)", res,
                        f"d used below = {d:.12g}")
                # 2. E0, E1: Lemma hypothesis at E = recorded lower bound with the verified d
                A, X = C.block_ops(cas, W, g2, iP, len(lat.plaqs), shell=kw["shell"])
                for k in (0, 1):
                    E = r["lo"][k]
                    if not np.isfinite(E):
                        continue
                    res = I.at_most((A - X / (d - E)).toarray(), E, k) if E < d else dict(
                        n_lo=None, n_hi=None, delta=0.0, r_lo=0.0, r_hi=0.0, method="-", verdict="failed (E >= d)")
                    add(rows, name, tj, g2, k, f"E{k}", E, A.shape[0], "A - X/(d-E) at E = recorded lower", res,
                        f"d={d:.12g} ({r['dsrc']})")
                print(f"  {time.time()-t:.1f}s", flush=True)
    common.write("inertia_c5.csv", rows)


if __name__ == "__main__":
    main()
