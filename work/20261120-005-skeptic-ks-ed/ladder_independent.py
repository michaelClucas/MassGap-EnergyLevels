"""Independent open 2-plaquette ladder ED (no import from src/).

Reduction: the 4 corner vertices are 2-valent, so each outer 3-link path carries one spin and its
3 link Casimirs are equal. Gauge-invariant states are functions on the theta graph: outer paths
with holonomies P1, P2 (spins a, b, weight 3 each) and the shared middle link (spin c, weight 1).
Gauge-fixing the middle link to 1 leaves f(P1, P2) invariant under simultaneous conjugation, with
    f_abc = N sum_M sum C^{cM}_{a al, b be} C^{cM}_{a al', b be'} D^a_{al al'}(P1) D^b_{be be'}(P2),
    N^2 = (2a+1)(2b+1)/(2c+1).
E_mid^2 = Casimir of the simultaneous right action = c(c+1). Plaquette traces are Tr P1, Tr P2.
Reflection is (P1, P2) -> (P2, P1), i.e. (a, b, c) -> (b, a, c) with sign +1 (two CG swaps).
Matrix elements of Tr P1 come from the CG series and Haar orthogonality (Racah-formula CGs coded
here), then are checked pointwise at random group elements.
H = (g^2/2)[3a(a+1) + 3b(b+1) + c(c+1)] - (2/g^2)(Tr P1 + Tr P2).
Spins are stored as twice-values (ints).
"""
import csv
import functools
print = functools.partial(print, flush=True)
import itertools
from functools import lru_cache
from math import factorial, sqrt

import numpy as np
from scipy.linalg import eigh, expm


@lru_cache(None)
def cg(j1, m1, j2, m2, J, M):
    """<j1 m1 j2 m2|J M>, args are twice-values. Racah formula, Condon-Shortley phase."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2 or (j1 + j2 + J) % 2:
        return 0.0
    if (j1 + m1) % 2 or (j2 + m2) % 2 or (J + M) % 2:
        return 0.0
    f = lambda t2: factorial(t2 // 2)
    pre = sqrt((J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J) / f(j1 + j2 + J + 2))
    pre *= sqrt(f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2))
    s = 0.0
    for k in range(0, 400):
        a = [k, (j1 + j2 - J) // 2 - k, (j1 - m1) // 2 - k, (j2 + m2) // 2 - k,
             (J - j2 + m1) // 2 + k, (J - j1 - m2) // 2 + k]
        if a[1] < 0 or a[2] < 0 or a[3] < 0:
            break
        if min(a) < 0:
            continue
        den = 1
        for x in a:
            den *= factorial(x)
        s += (-1) ** k / den
    return pre * s


def ms(tj):
    return list(range(-tj, tj + 1, 2))


@lru_cache(None)
def T_tensor(ta, tb, tc):
    """T[al, be, al', be'] = sum_M C^{cM}_{a al b be} C^{cM}_{a al' b be'}."""
    da, db = ta + 1, tb + 1
    T = np.zeros((da, db, da, db))
    for i, al in enumerate(ms(ta)):
        for k, be in enumerate(ms(tb)):
            c1 = cg(ta, al, tb, be, tc, al + be)
            if c1 == 0:
                continue
            for i2, al2 in enumerate(ms(ta)):
                be2 = al + be - al2
                if abs(be2) > tb:
                    continue
                k2 = ms(tb).index(be2)
                T[i, k, i2, k2] += c1 * cg(ta, al2, tb, be2, tc, al + be)
    return T


def basis(tjmax):
    return [(ta, tb, tc) for ta in range(tjmax + 1) for tb in range(tjmax + 1)
            for tc in range(abs(ta - tb), min(ta + tb, tjmax) + 1, 2)]


def trace_on_first(ta, tb, tc, ta2, tc2):
    """<f_{a2 b c2}| Tr P1 |f_{abc}>, using D^{1/2}_{ss} D^a_{al al'} = sum_a2 C C D^{a2}."""
    T = T_tensor(ta, tb, tc)
    T2 = T_tensor(ta2, tb, tc2)
    da2 = ta2 + 1
    X = np.zeros((da2, tb + 1, da2, tb + 1))
    for s in (-1, 1):
        for i, al in enumerate(ms(ta)):
            c1 = cg(1, s, ta, al, ta2, s + al)
            if c1 == 0:
                continue
            mu = ms(ta2).index(s + al)
            for i2, al2 in enumerate(ms(ta)):
                c2 = cg(1, s, ta, al2, ta2, s + al2)
                if c2 == 0:
                    continue
                mu2 = ms(ta2).index(s + al2)
                X[mu, :, mu2, :] += c1 * c2 * T[i, :, i2, :]
    N = sqrt((ta + 1) * (tb + 1) / (tc + 1))
    N2 = sqrt((ta2 + 1) * (tb + 1) / (tc2 + 1))
    return N * N2 * np.sum(T2 * X) / ((ta2 + 1) * (tb + 1))


def build(tjmax):
    T_tensor.cache_clear()  # keep memory bounded across builds
    B = basis(tjmax)
    idx = {s: i for i, s in enumerate(B)}
    n = len(B)
    W = np.zeros((n, n))
    for (ta, tb, tc), i in idx.items():
        for da, dc in itertools.product((1, -1), repeat=2):
            ta2, tc2 = ta + da, tc + dc
            if (ta2, tb, tc2) in idx:
                v = trace_on_first(ta, tb, tc, ta2, tc2)
                W[idx[(ta2, tb, tc2)], i] += v                 # Tr P1
                W[idx[(tb, ta2, tc2)], idx[(tb, ta, tc)]] += v  # Tr P2 via f_abc(P1,P2) = f_bac(P2,P1)
    cas = np.array([3 * ta * (ta + 2) / 4 + 3 * tb * (tb + 2) / 4 + tc * (tc + 2) / 4 for ta, tb, tc in B])
    P = np.zeros((n, n))
    for (ta, tb, tc), i in idx.items():
        P[idx[(tb, ta, tc)], i] = 1.0
    return B, cas, W, P


# ---- pointwise check with explicit D matrices built here
def spin_mats(tj):
    """Standard spin matrices in the ascending-m basis ms(tj)."""
    m = np.array(ms(tj)) / 2
    j = tj / 2
    jp = np.zeros((tj + 1, tj + 1))
    for k in range(tj):  # J+ |m_k> = c |m_{k+1}>
        jp[k + 1, k] = sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    return (jp + jp.T) / 2, (jp - jp.T) / 2j, np.diag(m)


def Dmat(tj, v):
    jx, jy, jz = spin_mats(tj)
    return expm(-1j * (v[0] * jx + v[1] * jy + v[2] * jz))


def f_val(ta, tb, tc, P1, P2):
    N = sqrt((ta + 1) * (tb + 1) / (tc + 1))
    return N * np.einsum('abcd,ac,bd->', T_tensor(ta, tb, tc), Dmat(ta, P1), Dmat(tb, P2))


def pointwise_check(tjmax=4, npts=6, seed=1):
    rng = np.random.default_rng(seed)
    B, cas, W, P = build(tjmax + 1)
    idx = {s: i for i, s in enumerate(B)}
    err = 0.0
    for _ in range(npts):
        P1, P2 = rng.normal(size=3) * 2, rng.normal(size=3) * 2
        vals = np.array([f_val(*s, P1, P2) for s in B])
        trs = np.trace(Dmat(1, P1)) + np.trace(Dmat(1, P2))
        for s in B:
            if max(s) > tjmax:
                continue
            err = max(err, abs(trs * vals[idx[s]] - W[:, idx[s]] @ vals))
    return err


def norm_check(npts=4000, seed=2):
    """Monte Carlo Haar check of orthonormality on a few states (rough, ~1/sqrt(npts))."""
    rng = np.random.default_rng(seed)

    def haar():
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        th = np.arccos(np.clip(q[0], -1, 1)); n = q[1:] / (np.linalg.norm(q[1:]) + 1e-300)
        return 2 * th * n
    S = [(0, 0, 0), (1, 1, 0), (1, 1, 2), (2, 2, 2), (2, 1, 1)]
    G = np.zeros((len(S), len(S)), complex)
    for _ in range(npts):
        P1, P2 = haar(), haar()
        v = np.array([f_val(*s, P1, P2) for s in S])
        G += np.outer(v.conj(), v)
    return np.abs(G / npts - np.eye(len(S))).max()


def gaps(B, cas, W, P, g2):
    evp, Up = eigh(P)
    H = 0.5 * g2 * np.diag(cas) - (2 / g2) * W
    Qp, Qm = Up[:, np.isclose(evp, 1)], Up[:, np.isclose(evp, -1)]
    wp = np.linalg.eigvalsh(Qp.T @ H @ Qp)
    wm = np.linalg.eigvalsh(Qm.T @ H @ Qm)
    return (Qp.shape[1], Qm.shape[1]), wp, wm


if __name__ == "__main__":
    print(f"pointwise identity (TrP1+TrP2) f = sum W f, max error: {pointwise_check():.2e}")
    print(f"Monte Carlo orthonormality (4000 Haar pts, expect ~2e-2): {norm_check():.2e}")
    B, cas, W, P = build(8)
    print(f"|W-W^T| = {abs(W - W.T).max():.1e}, |[P,W]| = {abs(P @ W - W @ P).max():.1e}")

    ref = {}
    with open("work/20261120-003-su2-ks-ed/results.csv") as fh:
        for r in csv.DictReader(fh):
            if r["lattice"] == "ladder2":
                ref[(r["sector"], float(r["g2"]), float(r["j_max"]))] = (int(r["dim"]), float(r["E0"]), float(r["E1"]))
    worst, nmatch = 0.0, 0
    for tjmax in range(1, 17):
        B, cas, W, P = build(tjmax)
        for g2 in (0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
            dims, wp, wm = gaps(B, cas, W, P, g2)
            for name, w, d in (("P+", (wp[0], wp[1]), dims[0]), ("P-", (wp[0], wm[0]), dims[1])):  # csv: P- row stores (E0 global, lowest P-)
                key = (name, g2, tjmax / 2)
                if key in ref:
                    d0, e0, e1 = ref[key]
                    assert d0 == len(B), ("dim mismatch", key, len(B), d0)  # results.csv stores the total dim
                    worst = max(worst, abs(w[0] - e0), abs(w[1] - e1)); nmatch += 1
    print(f"vs results.csv: {nmatch} (sector, g2, j_max) points, dims identical, max |E - E_ref| = {worst:.2e}")

    print("\nconvergence beyond the original j_max=8 (gaps measured from E0 in P+):")
    for tjmax in (16, 18, 20):
        B, cas, W, P = build(tjmax)
        row = []
        for g2 in (0.2, 0.3, 0.5):
            dims, wp, wm = gaps(B, cas, W, P, g2)
            row.append(f"g2={g2}: P+ {wp[1]-wp[0]:.8f} P- {wm[0]-wp[0]:.8f}")
        print(f"  j_max={tjmax/2:4.1f} dim={len(B)}: " + " | ".join(row))

    print(f"\nweak coupling, j_max=10; harmonic 2sqrt3={2*sqrt(3):.6f}, sqrt3+sqrt5={sqrt(3)+sqrt(5):.6f}")
    for g2 in (0.03, 0.05, 0.1, 0.2):
        dims, wp, wm = gaps(B, cas, W, P, g2)
        gp, gm = wp[1] - wp[0], wm[0] - wp[0]
        print(f"  g2={g2}: P+ {gp:.6f} (2sqrt3-D)/g2={(2*sqrt(3)-gp)/g2:.4f};  "
              f"P- {gm:.6f} (sqrt3+sqrt5-D)/g2={(sqrt(3)+sqrt(5)-gm)/g2:.4f}")
