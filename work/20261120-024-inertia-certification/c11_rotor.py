"""C11 (work 015): inertia check of every finite certified lower bound of the O(3) rotor chain whose
matrices exceeded the dense limit (n > 1500, eigsh used), 1+1D Hamiltonian lattice, finite chain,
fixed g. For each row (L, bc, l_max, g) of rotor_certified*.csv:
  1. Floor d: recomputed as in rotor.certify. Each depth-1 sub-certificate L0(H_alpha) (eigsh) is
     checked at the Lemma hypothesis: #eig(A_a - X/(d_a - L0) < L0) = 0. Unverified ones are dropped
     (Anderson bound only) and d recomputed.
  2. Each recorded E_k lower bound (k = 0, 1): E < d and #eig(A - X/(d - E) < E) <= k.
Usage (repo root): python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/c11_rotor.py <mindim> <maxdim>
Appends to inertia_c11.csv (rows already present are skipped).
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W015 = ROOT / "work" / "20261120-015-line6-rotor-chain"
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE)); sys.path.insert(0, str(W015))
import rotor  # noqa: E402
from src import inertia as I  # noqa: E402
import common  # noqa: E402
import memguard  # noqa: E402

OUTF = HERE / "inertia_c11.csv"
ALPHAS = (0.3, 0.6, 0.85)


def emit(row):
    new = not OUTF.exists()
    with open(OUTF, "a", newline="") as fh:
        w = csv.DictWriter(fh, common.FIELDS)
        if new:
            w.writeheader()
        w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in row.items()})


def on_abort(label, mb):
    L, bc, lmax, g, k, dim = label
    emit(dict(claim="C11", work="015", system=f"O(3) rotor chain L={L} {bc}", setting="1+1D lattice, finite chain, M=0",
              g=g, j_or_l_max=lmax, level_k=k, quantity=f"E{k}" if k >= 0 else "floor", dim=dim, lanczos_used="yes",
              verdict="too large to factor", note=f"SuperLU fill exceeded RAM guard ({mb:.0f} MB > 4600 MB)"))


def floor_d(L, bc, gk, gv, lmax, depth, label):
    """Returns (d, notes) with every sub-certificate verified by inertia."""
    T, V, X = rotor.parts(L, bc, lmax)
    c = (lmax + 1) * (lmax + 2)
    best, notes = -np.inf, []
    for al in ALPHAS:
        lb = rotor.anderson_E0_lower(L, bc, gk * (1 - al), gv)
        if depth > 0:
            gk2 = gk * (1 - al)
            sub = rotor.certify(L, bc, gk2, gv, lmax, 0, ALPHAS, want_k=1)
            E = sub["L0"]
            if np.isfinite(E) and E > lb:
                d2 = sub["floor"]
                A2 = gk2 * T - gv * V
                memguard.state["label"] = label
                res = I.at_most(A2 - 2 * gv**2 * X / (d2 - E), E, 0, dense_max=1500) if E < d2 else dict(verdict="failed")
                notes.append(f"a={al}: subL0={E:.8g} {res['verdict']}")
                if res["verdict"] == "verified":
                    lb = E
        best = max(best, al * gk * c + lb)
    return best, notes


KEYROWS = "--key" in sys.argv
_key = None


def keyrows(_r, _fn):
    """Per (L, bc, l_max): the row with the smallest g whose gap lower bound is positive (the
    015 summary table entry) and the g = 2 row."""
    global _key
    if _key is None:
        best = {}
        _key = set()
        for fn in ("rotor_certified.csv", "rotor_certified_large.csv"):
            for r in csv.DictReader(open(W015 / fn)):
                c = (int(r["L"]), r["bc"], int(r["lmax"]))
                g = float(r["g"])
                if float(r["gap_lo"]) > 0 and (c not in best or g < best[c]):
                    best[c] = g
                if g == 2.0:
                    _key.add(c + (g,))
        _key |= {c + (g,) for c, g in best.items()}
    return _key


def main():
    mindim, maxdim = int(sys.argv[1]), int(sys.argv[2])
    done = set()
    if OUTF.exists():
        done = {(r["system"], r["g"], r["j_or_l_max"], r["level_k"]) for r in csv.DictReader(open(OUTF))}
    memguard.start(4600, on_abort)
    rows = []
    for fn, large in (("rotor_certified.csv", False), ("rotor_certified_large.csv", True)):
        for r in csv.DictReader(open(W015 / fn)):
            dim = int(r["dim"])
            if not (mindim <= dim <= maxdim) or dim <= 1500:
                continue
            ks = [k for k in (0, 1) if np.isfinite(float(r[f"E{k}_lo"]))]
            if ks and (not KEYROWS or (int(r["L"]), r["bc"], int(r["lmax"]), float(r["g"])) in keyrows(r, fn)):
                rows.append((int(r["L"]), r["bc"], int(r["lmax"]), float(r["g"]), dim, ks,
                             {k: float(r[f"E{k}_lo"]) for k in ks}, float(r["floor"])))
    for L, bc, lmax, g, dim, ks, Elo, frec in rows:
        system = f"O(3) rotor chain L={L} {bc}"
        if all((system, f"{g:.10g}", str(lmax), str(k)) in done for k in ks):
            continue
        t = time.time()
        gk, gv = g / 2, 1 / g
        depth = 1 if L * lmax <= 12 or lmax == 1 and L <= 6 else 0
        d, notes = floor_d(L, bc, gk, gv, lmax, depth, (L, bc, lmax, g, -1, dim))
        T, V, X = rotor.parts(L, bc, lmax)
        A = gk * T - gv * V
        Xs = 2 * gv**2 * X
        for k in ks:
            E = Elo[k]
            memguard.state["label"] = (L, bc, lmax, g, k, dim)
            if E < d:
                res = I.at_most(A - Xs / (d - E), E, k, dense_max=1500)
            else:
                res = dict(n_lo=None, n_hi=None, delta=0.0, r_lo=0.0, r_hi=0.0, method="-", verdict="failed (E >= d)")
            row = dict(claim="C11", work="015", system=system, setting="1+1D lattice, finite chain, M=0 sector",
                       g=g, j_or_l_max=lmax, level_k=k, quantity=f"E{k}", claimed_lower=E, dim=dim,
                       lanczos_used="yes (eigsh)", matrix="A - X/(d-E) at E = claimed lower", allowed=k,
                       **common.res_fields(res),
                       note=f"d={d:.8g} (csv {frec:g}, depth {depth}); " + "; ".join(notes))
            emit(row)
            print(f"L={L} {bc} lmax={lmax} g={g} E{k}>={E:.10g} dim {dim}: count in [{res['n_lo']},{res['n_hi']}] "
                  f"allowed {k} -> {res['verdict']} (r={np.nanmax([res['r_lo'], res['r_hi']]):.1e}, {res['method']}) "
                  f"d={d:.6g} csv {frec:g} | {'; '.join(notes)} | {time.time()-t:.1f}s", flush=True)


if __name__ == "__main__":
    main()
