"""C9 (work 012): inertia check of every finite certified lower bound on the SU(2) KS 2x2 torus
(2+1D, finite lattice, fixed g^2), per Z2 flux sector and momentum block, and of refine.csv.
The work-012 lower bounds used dense LAPACK on blocks of dim <= 2472 (no Lanczos); the Ritz upper
bounds from the j_max = 2 blocks used eigsh but are min-max upper bounds for any subspace.
Check, directly at the Lemma hypothesis (work 006): for the recorded lower bound E of level k,
E < d and #eig(A - X/(d - E) < E) <= k, with A, X rebuilt from the work-012 blocks and
d = max(crude, split floor) recomputed exactly as in 012/analyze.py.
Needs artifacts/024-inertia/torus_blocks.npz (extract_torus.py).
Usage (repo root): python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c9_torus.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
from src import certify as C  # noqa: E402
from src import torus_sym as T  # noqa: E402
from src import inertia as I  # noqa: E402
import common  # noqa: E402

W012 = ROOT / "work" / "20261120-012-torus-sectors"
KEY = {"(0,0)": "00", "(1,0)=(0,1)": "10", "(1,1)": "11"}
KM = {"(0,0)": "00", "(pi,0)": "10", "(0,pi)": "01", "(pi,pi)": "11"}


def load():
    z = np.load(ROOT / "artifacts" / "024-inertia" / "torus_blocks.npz")
    out = {}
    for name in {k[:-len("_cas")] for k in z.files if k.endswith("_cas")}:
        def m(nm):
            return sp.csr_matrix((z[f"{name}_{nm}_v"], (z[f"{name}_{nm}_r"], z[f"{name}_{nm}_c"])),
                                 shape=tuple(z[f"{name}_{nm}_shape"]))
        out[name] = (z[f"{name}_cas"], m("WP"), m("WS"))
    return out


def ops(blk, g2):
    cas, WP, WS = blk
    A = (sp.diags(0.5 * g2 * cas) - (2.0 / g2) * WP).toarray()
    X = ((4.0 / g2**2) * (WS.T @ WS)).toarray()
    return A, X


def floor(g2, c_Q):
    crude = 0.5 * g2 * c_Q - 4.0 * 4 / g2
    f, _ = T.split_floor(g2, c_Q, ts=np.linspace(0, 1, 201))
    return max(crude, f)


def main():
    blocks = load()
    lat = T.torus2x2()
    cQ = {"00": C.casimir_floor(lat, 3)[0], "10": C.casimir_floor(lat, 2)[0], "11": C.casimir_floor(lat, 2)[0]}
    print("c_Q", cQ, flush=True)
    todo = []
    for r in csv.DictReader(open(W012 / "results.csv")):
        if r["kind"] != "certified":
            continue
        for k, col in ((0, "E0_lo"), (1, "E1_lo")):
            if r[col] and np.isfinite(float(r[col])):
                todo.append(("results.csv", r["sector"], r["momentum"], float(r["g2"]), k, float(r[col]),
                             float(r["d_floor"]), r["method"].split(";")[0]))
    for r in csv.DictReader(open(W012 / "refine.csv")):
        todo.append(("refine.csv", "(0,0)", "(0,0)", float(r["g2"]), int(r["level"]), float(r["lo"]),
                     float(r["d_floor"]), "bisect (refine.py)"))
    rows = []
    cache = {}
    for src, sec, mom, g2, k, E, d_rec, meth in todo:
        fk = KEY[sec]
        name = f"f{fk}_k{KM[mom]}"
        if name not in blocks and fk == "10":  # (0,pi) of (1,0) stored as such; all four k present
            raise KeyError(name)
        if (name, g2) not in cache:
            cache = {(name, g2): ops(blocks[name], g2)}
        A, X = cache[(name, g2)]
        d = floor(g2, cQ[fk])
        assert abs(d - d_rec) < 1e-8 * max(1, abs(d)), (d, d_rec)
        if not E < d:
            verdict = dict(n_lo=None, n_hi=None, delta=0.0, r_lo=0.0, r_hi=0.0, method="-", verdict="failed (E >= d)")
        else:
            verdict = I.at_most(A - X / (d - E), E, k)
        rows.append(dict(claim="C9", work="012/" + src, system="periodic 2x2 torus, flux " + sec + ", k=" + mom,
                         setting="SU(2) KS 2+1D finite lattice", g=g2,
                         j_or_l_max=1.5 if fk == "00" else 1.0, level_k=k, quantity=f"E{k} of block",
                         claimed_lower=E, dim=A.shape[0], lanczos_used="no (dense LAPACK); upper bounds Ritz",
                         matrix="A - X/(d-E) at E = claimed lower", allowed=k, **common.res_fields(verdict),
                         note=f"{meth}; d={d:.10g}"))
        print(f"{src} {sec} k={mom} g2={g2} E{k} >= {E:.10g}: dim {A.shape[0]}, count in "
              f"[{verdict['n_lo']},{verdict['n_hi']}] allowed {k} -> {verdict['verdict']} "
              f"(r={np.nanmax([verdict['r_lo'], verdict['r_hi']]):.1e})", flush=True)
    common.write("inertia_c9.csv", rows)


if __name__ == "__main__":
    main()
