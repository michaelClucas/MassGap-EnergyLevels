"""Independent exact-rational certification of E0, E1 of the UNTRUNCATED single-plaquette SU(2) KS
Hamiltonian (2+1D, one plaquette, finite lattice, a = 1, fixed g^2), and comparison with the
certified intervals of work/20261120-006 (results_plaq1.csv).

Written from scratch (does not import src/). Gauge-invariant states of one plaquette are class
functions of U_p; orthonormal basis chi_j, j = 0, 1/2, 1, ... (Peter-Weyl). All four links carry spin
j, so sum_l E_l^2 = 4 j(j+1); Tr U chi_j = chi_{j+1/2} + chi_{j-1/2} (Clebsch-Gordan with spin 1/2).
    H = 2 g^2 j(j+1) delta - (2/g^2) [adjacent j]      (tridiagonal, H as in src/ks_su2.py)

Bounds, all in exact rational arithmetic (fractions.Fraction), no floating point:
 - upper: lambda_k(H) <= lambda_k(A_n) (min-max), A_n = truncation to j <= n/2. Certified by a Sylvester
   inertia count: #neg pivots of (A_n - E) >= k+1  =>  lambda_k(A_n) < E.
 - lower: Feshbach/Schur. Q = span{chi_j, j > n/2}; QHQ >= 2 g^2 j1(j1+1) - (2/g^2)*2 =: d,
   j1 = (n+1)/2 (|Tr U| <= 2). PHQ couples only chi_{n/2} to chi_{(n+1)/2} with weight -2/g^2,
   so PHQ QHP = (4/g^4) e e^T exactly. If E < d and #neg pivots of (A_n - (4/g^4)/(d-E) e e^T - E)
   <= k, then lambda_k(H) >= E.  (Lemma of 006 README; re-proved in this folder's README.)
Bisection on E over rationals to width TOL.
"""
import csv
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOL = F(1, 10**13)


def neg_pivots(diag, off):
    """Sylvester: number of negative eigenvalues of the symmetric tridiagonal matrix (diag, off).
    Exact LDL^T recurrence; returns None on an exactly zero pivot (caller aborts)."""
    cnt, prev = 0, None
    for i, a in enumerate(diag):
        p = a if i == 0 else a - off[i - 1] ** 2 / prev
        if p == 0:
            return None
        cnt += p < 0
        prev = p
    return cnt


def mats(g2, n):
    diag = [2 * g2 * F(t, 2) * (F(t, 2) + 1) for t in range(n + 1)]
    off = [-2 / g2] * n
    return diag, off


def count_upper(g2, n, E):
    diag, off = mats(g2, n)
    return neg_pivots([a - E for a in diag], off)


def count_lower(g2, n, E, d):
    diag, off = mats(g2, n)
    diag = [a - E for a in diag]
    diag[-1] -= (4 / g2**2) / (d - E)
    return neg_pivots(diag, off)


def bisect(pred, lo, hi):
    """pred monotone: False at lo, True at hi. Returns (lo, hi) bracketing the switch."""
    while hi - lo > TOL:
        mid = (lo + hi) / 2
        mid = F(round(mid * 2**50), 2**50)  # keep denominators bounded
        if mid <= lo or mid >= hi:
            mid = (lo + hi) / 2
        while True:  # exact zero pivot: move the test point (it stays strictly inside (lo, hi))
            try:
                r = pred(mid)
                break
            except AssertionError:
                mid += (hi - mid) / 1000
        if r:
            hi = mid
        else:
            lo = mid
    return lo, hi


def certify(g2, n, k):
    """Return (L, U): exact rationals with L <= lambda_k(H) <= U."""
    j1 = F(n + 1, 2)
    d = 2 * g2 * j1 * (j1 + 1) - 4 / g2
    big = 2 * g2 * F(n, 2) * (F(n, 2) + 1) + 10 / g2 + 10
    # upper: smallest E with count_upper >= k+1 (count_upper(E) >= k+1 => lambda_k(A) < E)
    def up(E):
        c = count_upper(g2, n, E)
        assert c is not None, "exact zero pivot"
        return c >= k + 1
    _, U = bisect(up, -big, big)
    if U >= d:
        return None, U
    # lower: largest E < d with count_lower(E) <= k (=> lambda_k(H) >= E)
    def low_fails(E):
        c = count_lower(g2, n, E, d)
        assert c is not None, "exact zero pivot"
        return c > k
    L, _ = bisect(low_fails, -big, min(U, d - TOL))
    assert not low_fails(L)
    return L, U


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40  # 2*j_max for the exact reference
    g2s = [F(1, 5), F(1, 2), F(1), F(3, 2), F(2), F(3), F(5), F(10)]
    ref = {}
    out = Path(__file__).resolve().parent / "plaq_exact.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["g2", "two_jmax", "E0_lo", "E0_hi", "E1_lo", "E1_hi", "Delta_lo", "Delta_hi"])
        for g2 in g2s:
            L0, U0 = certify(g2, nmax, 0)
            L1, U1 = certify(g2, nmax, 1)
            ref[float(g2)] = (L0, U0, L1, U1)
            w.writerow([float(g2), nmax] + [f"{float(x):.15g}" for x in
                                            (L0, U0, L1, U1, L1 - U0, U1 - L0)])
            print(f"g2={float(g2):5}: E0 in [{float(L0):.15g}, {float(U0):.15g}]  "
                  f"E1 in [{float(L1):.15g}, {float(U1):.15g}]  Delta in "
                  f"[{float(L1-U0):.13g}, {float(U1-L0):.13g}]", flush=True)
    # compare with 006: each of their intervals must contain our (exact) interval for the truth
    bad, n = 0, 0
    worst = 0.0
    for r in csv.DictReader(open(ROOT / "work/20261120-006-certified-gap-bounds/results_plaq1.csv")):
        g2 = float(r["g2"])
        L0, U0, L1, U1 = ref[g2]
        for (lo_s, hi_s, L, U) in (("E0_lo", "E0_hi", L0, U0), ("E1_lo", "E1_hi", L1, U1)):
            lo, hi = float(r[lo_s]), float(r[hi_s])
            n += 1
            # sound iff lo <= true <= hi; true in [L, U]; violation certain iff lo > U or hi < L
            if lo > float(U) or hi < float(L):
                bad += 1
                print(f"  VIOLATION g2={g2} j_max={r['j_max']} {lo_s}: [{lo}, {hi}] vs exact "
                      f"[{float(L)}, {float(U)}]")
            if lo > -1e300:
                worst = max(worst, lo - float(L))  # >0 would mean their lower bound above ours
        # gap
        lo, hi = float(r["Delta_lo"]), float(r["Delta_hi"])
        n += 1
        if lo > float(U1 - L0) or hi < float(L1 - U0):
            bad += 1
            print(f"  VIOLATION Delta g2={g2} j_max={r['j_max']}: [{lo}, {hi}]")
    print(f"006 plaq1 intervals checked: {n}; violations: {bad}; "
          f"max (their lower - exact lower) = {worst:.3e}")


if __name__ == "__main__":
    main()
