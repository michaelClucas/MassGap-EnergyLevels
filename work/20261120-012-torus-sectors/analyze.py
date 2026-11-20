"""2x2 torus, SU(2) KS (2+1D, finite lattice, a = 1): spectra per Z2 electric-flux sector and
momentum, C4v irreps at k = (0,0), (pi,pi) in the zero-flux sector, and certified bounds.

Inputs: artifacts/012-torus-sectors/blocks_f??_tj{top}.pkl from build.py (top = largest 2*j_max).
(0,1) is obtained from (1,0) by the diagonal reflection sigma_d (a unitary symmetry of H that
maps the two sectors, see validate.py), so it is not rebuilt at top = 4.

Certified bounds per (sector, k) block (Lemma of work/20261120-006, applied to the block; P, Q and
the shell are translation invariant, so the Lemma holds block by block):
  P = all states with spins <= (top-1)/2, exact shell X from the spin-top/2 states (dense, dims
  <= 2500), floor d = max(crude, split) with c_Q = casimir_floor(top-1);
  upper bounds: Ritz values of the top-j_max block (min-max), min with the P-block ones.
If U_k >= d the standard choice fails; then bisection for the largest E < d with
lambda_k(A - X/(d-E)) >= E (still the Lemma; gives a weaker but valid bound).
Usage (repo root): python work/20261120-012-torus-sectors/analyze.py <top> [--quick]
"""
import csv
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import certify as C  # noqa: E402
from src import torus_sym as T  # noqa: E402

ART = ROOT / "artifacts" / "012-torus-sectors"
OUT = Path(__file__).resolve().parent
KN = {(0, 0): "(0,0)", (1, 0): "(pi,0)", (0, 1): "(0,pi)", (1, 1): "(pi,pi)"}
SECTORS = [(0, 0), (1, 0), (1, 1)]
SEC_NAME = {(0, 0): "(0,0)", (1, 0): "(1,0)=(0,1)", (1, 1): "(1,1)"}
FIELDS = ["kind", "sector", "momentum", "irrep", "g2", "j_max", "dim", "E0", "E1",
          "E0_lo", "E0_hi", "E1_lo", "E1_hi", "d_floor", "method"]
NPLAQ = 4


def lowest(H, k):
    if H.shape[0] <= C.DENSE_MAX:
        return sla.eigh(H.toarray(), subset_by_index=[0, k - 1])
    w, v = spla.eigsh(H, k=k, which="SA", tol=1e-11)
    o = np.argsort(w)
    return w[o], v[:, o]


def ham(cas, W, g2):
    return (sp.diags(0.5 * g2 * cas) - (2.0 / g2) * W).tocsr()


def irrep_of(sx, sd, dim):
    if dim == 1:
        return {(1, 1): "A1", (-1, -1): "A2", (1, -1): "B1", (-1, 1): "B2"}.get(
            (int(round(sx)), int(round(sd))), "?")
    if dim == 2 and abs(sx) < 1e-6 and abs(sd) < 1e-6:
        return "E"
    return f"deg{dim}(sx={sx:.2f},sd={sd:.2f})"


def point_ops(lat, reps, k):
    """sigma_x, sigma_d on a k-block (k invariant under both): R|r_k> = sum_s a_s chi_k(h_s)|rep_s,k>."""
    pos = {r: i for i, r in enumerate(reps)}
    mats = {}
    for name, op in (("sx", T.sigma_x()), ("sd", T.sigma_d())):
        rows, cols, vals = [], [], []
        for i, r in enumerate(reps):
            for s, a in T.apply_op(lat, op, r):
                rr, h, _ = T.canon(s)
                if rr not in pos:  # |rr_k> = 0 (character nontrivial on its stabilizer)
                    continue
                rows.append(pos[rr]); cols.append(i); vals.append(np.real(a) * T.chi(k, h))
        mats[name] = sp.csr_matrix((vals, (rows, cols)), shape=(len(reps),) * 2)
    return mats


def classify(mats, vecs, w, tol=1e-6):
    labels = [None] * len(w)
    i = 0
    while i < len(w):
        j = i + 1
        while j < len(w) and abs(w[j] - w[i]) < tol * max(1, abs(w[i])):
            j += 1
        V = vecs[:, i:j]
        sx, sd = (np.trace(V.T @ (mats[n] @ V)) for n in ("sx", "sd"))
        lab = irrep_of(sx, sd, j - i)
        if j == len(w) and j - i == 1 and lab == "?":
            lab = "?(cut)"
        for m in range(i, j):
            labels[m] = lab
        i = j
    return labels


def certify_block(cas, W, reps, top, g2, c_Q, kk, floor_cache, useful_min=-np.inf):
    maxsp = np.array([max(r[0]) for r in reps])
    iP = np.where(maxsp <= top - 1)[0]
    ish = np.where(maxsp == top)[0]
    A, X = C.block_ops(cas, W, g2, iP, NPLAQ, shell=ish)
    if (g2, c_Q) not in floor_cache:
        crude = 0.5 * g2 * c_Q - 4.0 * NPLAQ / g2
        f, t = T.split_floor(g2, c_Q, ts=np.linspace(0, 1, 201))
        floor_cache[(g2, c_Q)] = (max(crude, f), "split(t=%.3f)" % t if f > crude else "crude")
    d, dsrc = floor_cache[(g2, c_Q)]
    hiP = C.upper_bounds(A, kk) + C._eps_margin(A)
    Hs = ham(cas, W, g2)
    hiT = C.upper_bounds(Hs, kk) + C._eps_margin(Hs)
    hi = np.minimum(hiP, hiT)
    lo = np.full(kk, -np.inf)
    meth = []
    Ad = A.toarray()
    Xd = X.toarray()

    def f(E, k):
        return C.lower_bound_k(Ad - Xd / (d - E), k)

    for k in range(kk):
        # largest E < d with f(E) >= E; f(E) <= lambda_k(A) <= hiP[k], so E <= min(d, hiP[k])
        if d <= useful_min:  # any bound obtained would lie below the ground-state upper bound
            meth.append("none(d<=E0_hi)")
            continue
        b = min(d, hiP[k]) - 1e-9
        a = f(b, k) if hiP[k] < d else -np.inf  # standard choice (L = f(U)) when U < d
        if np.isfinite(a):  # standard choice; bisection would gain only O(|X|/(d-U)^2)
            lo[k] = min(a, b)
            meth.append("std")
            continue
        if not np.isfinite(a):
            a = min(-1e3, d - 1e3)
        for _ in range(14):
            m = 0.5 * (a + b)
            if f(m, k) >= m:
                a = m
            else:
                b = m
        lo[k] = a
        meth.append("bisect")
    return lo, hi, d, dsrc, len(iP), meth


def fmt(v):
    return f"{v:.10g}" if isinstance(v, (float, np.floating)) else v


def main():
    top = int(sys.argv[1])
    quick = "--quick" in sys.argv
    g2s = [3.0] if quick else [1.0, 1.5, 2.0, 3.0]
    g2s_ed = [3.0] if quick else [1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    lat = T.torus2x2()
    print(f"top tj = {top} (j_max = {top/2})", flush=True)
    rows = []
    data = {}
    for fl in SECTORS:
        tf = top if (ART / f"blocks_f{fl[0]}{fl[1]}_tj{top}.pkl").exists() else top - 1
        with open(ART / f"blocks_f{fl[0]}{fl[1]}_tj{tf}.pkl", "rb") as fh:
            data[fl] = pickle.load(fh)
        data[fl]["top"] = tf
        data[fl]["c_Q"] = C.casimir_floor(lat, tf - 1)[0]
        print(f"flux {fl}: dim {data[fl]['dim']}, blocks "
              f"{[len(b[0]) for b in data[fl]['blocks'].values()]}, j_max {data[fl]['top']/2}, "
              f"c_Q(P) = {data[fl]['c_Q']}", flush=True)
    pops = {k: point_ops(lat, data[(0, 0)]["blocks"][k][0], k) for k in ((0, 0), (1, 1))}
    # ---------------- ED spectra by block
    spec = {}
    for g2 in g2s_ed:
        t = time.time()
        for fl in SECTORS:
            for k, (reps, cas, W, S) in data[fl]["blocks"].items():
                cls = fl == (0, 0) and k in pops
                nev = 6 if cls else 2
                w, v = lowest(ham(cas, W, g2), nev)
                labels = classify(pops[k], v, w) if cls else [""] * nev
                spec[(g2, fl, k)] = (w, labels)
                rows.append(dict(kind="ed", sector=SEC_NAME[fl], momentum=KN[k],
                                 irrep="/".join(labels[:2]) if cls else "", g2=g2, j_max=data[fl]["top"] / 2,
                                 dim=len(reps), E0=w[0], E1=w[1],
                                 method="translation-blocked ED; E0,E1 = two lowest in block"))
                if cls:
                    for m in range(2, nev):
                        rows.append(dict(kind=f"ed-level{m}", sector=SEC_NAME[fl], momentum=KN[k],
                                         irrep=labels[m], g2=g2, j_max=top / 2, dim=len(reps),
                                         E0=w[m], method=f"level {m} in block"))
        print(f"g2={g2}: ED done {time.time()-t:.0f}s", flush=True)
    for g2 in g2s_ed:
        w00, l00 = spec[(g2, (0, 0), (0, 0))]
        E0 = w00[0]
        levels = sorted((e, fl, k, (lab[m] if lab[m] else ""))
                        for (gg, fl, k), (w, lab) in spec.items() if gg == g2
                        for m, e in enumerate(w[:2]))
        assert abs(levels[0][0] - E0) < 1e-9, "global ground state not in (0,0) k=0"
        e1 = levels[1]
        e1_00 = min(lv for lv in levels[1:] if lv[1] == (0, 0))
        print(f"g2={g2}: E0={E0:.6f}  Delta_all {e1[0]-E0:.6f} (flux {e1[1]}, k={KN[e1[2]]}); "
              f"Delta_00 {e1_00[0]-E0:.6f} (k={KN[e1_00[2]]} {e1_00[3]}); "
              f"k=0 levels {[round(x - E0, 5) for x in w00[1:]]} {l00[1:]}", flush=True)
        rows.append(dict(kind="gap-ed", sector="all", momentum=KN[e1[2]], irrep=f"flux {e1[1]}",
                         g2=g2, j_max=top / 2, E0=E0, E1=e1[0], method=f"Delta_all = {e1[0]-E0:.8g}"))
        rows.append(dict(kind="gap-ed", sector="(0,0)", momentum=KN[e1_00[2]], irrep=e1_00[3],
                         g2=g2, j_max=top / 2, E0=E0, E1=e1_00[0], method=f"Delta_00 = {e1_00[0]-E0:.8g}"))
        rows.append(dict(kind="gap-ed", sector="(0,0)", momentum="(0,0)", irrep=l00[1], g2=g2,
                         j_max=top / 2, E0=E0, E1=w00[1], method=f"Delta_00(k=0) = {w00[1]-E0:.8g}"))
    # ---------------- certified bounds
    floor_cache = {}
    for g2 in g2s:
        cert = {}
        e0hi = -np.inf
        for fl in SECTORS:
            for k, (reps, cas, W, S) in data[fl]["blocks"].items():
                kk = 2 if (fl, k) == ((0, 0), (0, 0)) else 1
                t = time.time()
                tf = data[fl]["top"]
                lo, hi, d, dsrc, dimP, meth = certify_block(cas, W, reps, tf, g2, data[fl]["c_Q"], kk,
                                                            floor_cache, e0hi)
                if (fl, k) == ((0, 0), (0, 0)):
                    e0hi = hi[0]
                cert[(fl, k)] = (lo, hi)
                r = dict(kind="certified", sector=SEC_NAME[fl], momentum=KN[k], g2=g2,
                         j_max=(tf - 1) / 2, dim=dimP, E0_lo=lo[0], E0_hi=hi[0], d_floor=d,
                         method=f"schur-exact-shell[{','.join(meth)}]; floor {dsrc}; "
                                f"upper from j_max {tf/2}")
                if kk == 2:
                    r.update(E1_lo=lo[1], E1_hi=hi[1])
                rows.append(r)
                print(f"g2={g2} flux {fl} k={KN[k]}: E0 in [{lo[0]:.6f}, {hi[0]:.6f}]"
                      + (f" E1 in [{lo[1]:.6f}, {hi[1]:.6f}]" if kk == 2 else "")
                      + f" d={d:.4f} ({dsrc}) {meth} {time.time()-t:.0f}s", flush=True)
        (lo00, hi00) = cert[((0, 0), (0, 0))]
        exc = [(lo00[1], hi00[1], (0, 0))]
        exc += [(lo[0], hi[0], fl) for (fl, k), (lo, hi) in cert.items() if (fl, k) != ((0, 0), (0, 0))]
        for name, sel in (("all", lambda fl: True), ("(0,0)", lambda fl: fl == (0, 0)),
                          ("(0,0) k=0", None)):
            ex = [(lo00[1], hi00[1], (0, 0))] if sel is None else [e for e in exc if sel(e[2])]
            l1, h1 = min(e[0] for e in ex), min(e[1] for e in ex)
            dlo, dhi = l1 - hi00[0], h1 - lo00[0]
            rows.append(dict(kind="gap-certified", sector=name, g2=g2, j_max=(top - 1) / 2,
                             E0_lo=lo00[0], E0_hi=hi00[0], E1_lo=l1, E1_hi=h1,
                             method=f"Delta in [{dlo:.8g}, {dhi:.8g}]"
                                    + ("" if dlo > 0 else " ; no certified positive lower bound")))
            print(f"g2={g2} certified Delta_{name} in [{dlo:.6g}, {dhi:.6g}]", flush=True)
    name = "results_quick.csv" if quick else "results.csv"
    with open(OUT / name, "w", newline="") as fh:
        wr = csv.DictWriter(fh, FIELDS)
        wr.writeheader()
        for r in rows:
            wr.writerow({k: fmt(v) for k, v in r.items()})
    print(f"wrote {name}")


if __name__ == "__main__":
    main()
