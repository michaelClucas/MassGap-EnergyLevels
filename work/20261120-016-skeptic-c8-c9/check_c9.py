"""Skeptic checks for C9 (work 012), independent of src/torus_sym.py and src/certify.py.

Setting: SU(2) KS, 2+1D, periodic 2x2 lattice (8 links, 4 plaquettes), finite volume, fixed g^2,
H = (g^2/2) sum_l j(j+1) - (2/g^2) sum_p Tr U_p.

1. c_Q: min sum_l j(j+1) over gauge-invariant spin configurations with max spin >= j_Q
   (own enumeration; vertex condition: sum of 2j even and 2 max <= sum, the condition for a
   nonzero SU(2) invariant in a 4-fold tensor product). Also per Z2 flux sector.
2. Plaquette count / link sharing on the 2x2 torus (each link in exactly 2 distinct plaquettes,
   each plaquette has 4 distinct links and 4 distinct corners).
3. Split floor d(g^2, c_Q) = max_t (1-t)(g^2/2)c_Q + 4 e(t), e(t) = lowest eigenvalue of the
   character-basis chain (diag t g^2 j(j+1), off-diag -2/g^2), dense truncation 2j <= 600
   (numerical cross-check only; truncation gives an upper estimate of e(t), the certified code
   uses a Schur tail bound).
4. Recompute every certified gap interval from the raw block bounds in results.csv / refine.csv.
"""
import csv
import itertools
from pathlib import Path

import numpy as np
import scipy.linalg as sla

HERE = Path(__file__).resolve().parent
W012 = HERE.parent / "20261120-012-torus-sectors"
L = 2


def vid(x, y):
    return (x % L) + L * (y % L)


def xl(x, y):
    return 2 * vid(x, y)


def yl(x, y):
    return 2 * vid(x, y) + 1


LEGS = [[xl(x, y), yl(x, y), xl(x - 1, y), yl(x, y - 1)] for y in range(L) for x in range(L)]
PLAQ = [[xl(x, y), yl(x + 1, y), xl(x, y + 1), yl(x, y)] for y in range(L) for x in range(L)]
CORNERS = [[vid(x, y), vid(x + 1, y), vid(x + 1, y + 1), vid(x, y + 1)] for y in range(L) for x in range(L)]
CUTX = [xl(0, y) for y in range(L)]
CUTY = [yl(x, 0) for x in range(L)]


def geometry():
    cnt = np.zeros(8, int)
    for p in PLAQ:
        assert len(set(p)) == 4
        for l in p:
            cnt[l] += 1
    assert all(len(set(c)) == 4 for c in CORNERS)
    # flux conservation: each plaquette meets each cut in 0 or 2 links
    mx = [len(set(p) & set(CUTX)) for p in PLAQ]
    my = [len(set(p) & set(CUTY)) for p in PLAQ]
    return cnt.tolist(), mx, my


def c_Q_enum(tjQ, tjcap=6):
    grids = np.indices([tjcap + 1] * 8, dtype=np.int16).reshape(8, -1)
    ok = grids.max(axis=0) >= tjQ
    for legs in LEGS:
        sub = grids[legs]
        s = sub.sum(axis=0)
        ok &= (s % 2 == 0) & (2 * sub.max(axis=0) <= s)
    G = grids[:, ok]
    cas = (G.astype(float) * (G + 2) / 4).sum(axis=0)
    fx = G[CUTX].sum(axis=0) % 2
    fy = G[CUTY].sum(axis=0) % 2
    out = {"any": cas.min()}
    for f in itertools.product((0, 1), repeat=2):
        m = (fx == f[0]) & (fy == f[1])
        out[f] = cas[m].min() if m.any() else np.inf
    return out


def e_plaq(kin, g2, J2=600):
    tj = np.arange(J2 + 1)
    return sla.eigh_tridiagonal(kin * tj * (tj + 2) / 4, -2.0 / g2 * np.ones(J2), select="i",
                                select_range=(0, 0), eigvals_only=True)[0]


def split_floor(g2, cQ, nt=2001):
    best = (-np.inf, None)
    for t in np.linspace(0, 1, nt):
        e = -4.0 / g2 if t == 0 else e_plaq(t * g2, g2)
        f = (1 - t) * 0.5 * g2 * cQ + 4 * e
        if f > best[0]:
            best = (f, t)
    return best


def f(x):
    return float(x) if x not in ("", None) else np.nan


def intervals():
    rows = [r for r in csv.DictReader(open(W012 / "results.csv")) if r["kind"] == "certified"]
    out = []
    for g2 in sorted({float(r["g2"]) for r in rows}):
        R = [r for r in rows if float(r["g2"]) == g2]
        gr = [r for r in R if r["sector"] == "(0,0)" and r["momentum"] == "(0,0)"][0]
        E0lo, E0hi, E1lo, E1hi = f(gr["E0_lo"]), f(gr["E0_hi"]), f(gr["E1_lo"]), f(gr["E1_hi"])
        others = [r for r in R if r is not gr]
        # global E0 must sit in the ground block: every other block's lower bound above E0_hi
        e0_in_ground = all(f(r["E0_lo"]) > E0hi for r in others)
        lo_all = min([E1lo] + [f(r["E0_lo"]) for r in others]) - E0hi
        hi_all = min([E1hi] + [f(r["E0_hi"]) for r in others]) - E0lo
        z = [r for r in others if r["sector"] == "(0,0)"]
        lo_00 = min([E1lo] + [f(r["E0_lo"]) for r in z]) - E0hi
        hi_00 = min([E1hi] + [f(r["E0_hi"]) for r in z]) - E0lo
        out.append((g2, e0_in_ground, lo_all, hi_all, lo_00, hi_00, E0lo, E0hi, E1lo, E1hi, gr["d_floor"]))
    ref = list(csv.DictReader(open(W012 / "refine.csv")))
    return out, ref


if __name__ == "__main__":
    cnt, mx, my = geometry()
    print(f"plaquettes per link {cnt}; plaquette-cut intersections x {mx} y {my}")
    for tjQ, lab in [(4, "max spin >= 2 (P = j<=3/2)"), (3, "max spin >= 3/2 (P = j<=1)")]:
        c = c_Q_enum(tjQ)
        print(f"c_Q {lab}: any sector {c['any']}; per flux " +
              ", ".join(f"{k}: {v}" for k, v in c.items() if k != "any"))
    print("# split floor (dense chain 2j<=600, t grid 2001) vs d_floor in results.csv")
    for g2, cQ, ref in [(1, 11, -7.861029116), (1.5, 11, -0.657391941), (2, 11, 4.319398566),
                        (3, 11, 12.0458277), (1, 7.5, -8.702509926), (1.5, 7.5, -2.676666269),
                        (2, 7.5, 1.273744133), (3, 7.5, 7.098860449)]:
        d, t = split_floor(g2, cQ)
        crude = 0.5 * g2 * cQ - 16.0 / g2
        print(f"g2={g2} c_Q={cQ}: split={d:.7f} (t={t:.4f}) csv={ref:.7f} diff={d-ref:+.1e} crude={crude:.4f}")
    print("# gap intervals recomputed from raw certified block bounds (results.csv)")
    out, ref = intervals()
    for (g2, ok, la, ha, l0, h0, a, b, c, d, dfl) in out:
        print(f"g2={g2}: E0 in ground block certified: {ok}; Delta_all in [{la:.7f}, {ha:.7f}]; "
              f"Delta_00 in [{l0:.7f}, {h0:.7f}]  (E0 [{a:.7f},{b:.7f}], E1 [{c:.7f},{d:.7f}])")
    print("# g2 = 1.5 with refine.csv lower bounds")
    r0, r1 = ref[0], ref[1]
    g15 = [o for o in out if o[0] == 1.5][0]
    E0hi_j2, E0hi_P = g15[7], float(r0["hi_P"])
    E1lo = float(r1["lo"])
    rows = [r for r in csv.DictReader(open(W012 / "results.csv"))
            if r["kind"] == "certified" and float(r["g2"]) == 1.5 and r["sector"] == "(0,0)" and r["momentum"] != "(0,0)"]
    kmin = min(float(r["E0_lo"]) for r in rows)
    print(f"E1 >= {E1lo}; min E0_lo over (0,0) k != 0 = {kmin:.6f}")
    print(f"gap_00 >= {min(E1lo, kmin) - E0hi_j2:.6f} using E0_hi = {E0hi_j2:.9f} (j_max 2 Ritz, results.csv)")
    print(f"gap_00 >= {min(E1lo, kmin) - E0hi_P:.6f} using E0_hi = {E0hi_P:.9f} (P Ritz, refine.csv hi_P)")
    print(f"upper: E1_hi - E0_lo(refine) = {g15[9] - float(r0['lo']):.6f}")
