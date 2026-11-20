"""Skeptic checks for C8 (work 011), written independently of hff2d.py. numpy/scipy + stdlib only.

Setting: 2D lattice gauge theory, Wilson measure, free-boundary box, H_FF = ground-state transform of
the link-gradient Langevin generator, gap(H_FF) = (1/4) x Poincare constant of
E(f) = sum_l int |grad_l f|^2 dpi_beta (work 007 fact 6). Unit S^3 metric (SU(2)), unit circle (U(1)).

A. Single plaquette SU(2): gap = lambda_1 of -(1/w)(w f')' on (0, pi), w = sin^2 e^{beta cos}.
   Independent cell-centred finite-volume discretisation.
B. U(1) Kramers test function, L -> infinity value (1/4) 4 Var(g')/Var(g), with the closed form
   Var(g') = 2 pi lam (1 - 2 pi lam) and Gauss-Legendre quadrature for Var(g).
C. SU(2) nested Wilson loops: moments a = I2/I1, c1 = I3/I1 from Bessel power series in 80-digit
   Decimal arithmetic; eigenvector from numpy; the Rayleigh quotient of that explicit vector is then
   re-evaluated in Decimal, so the printed number is the exact (to ~1e-50) quotient of an explicit
   test function, i.e. a genuine upper bound on gap(H_FF), independent of eigensolver conditioning.
"""
import sys
from decimal import Decimal as Dm, getcontext

import numpy as np
import scipy.linalg as sla
from scipy.special import ive

getcontext().prec = 80


# ---------------- A
def plaquette_gap(beta, n=20000):
    th_max = np.pi if beta < 20 else min(np.pi, 30 / np.sqrt(beta))
    h = th_max / n
    c = (np.arange(n) + 0.5) * h
    e = np.arange(1, n) * h
    lw = lambda t: 2 * np.log(np.sin(t)) + beta * (np.cos(t) - 1)
    wc = np.exp(lw(c)); we = np.exp(lw(e))
    kd = np.zeros(n); kd[:-1] += we / h; kd[1:] += we / h
    ko = -we / h
    s = 1 / np.sqrt(wc * h)
    vals = sla.eigh_tridiagonal(kd * s * s, ko * s[:-1] * s[1:], select="i", select_range=(0, 1),
                                eigvals_only=True)
    return vals[1]


# ---------------- B
def u1_bound_Linf(beta, npan=400, nq=40):
    x, wq = np.polynomial.legendre.leggauss(nq)
    edges = np.linspace(-np.pi, np.pi, npan + 1)
    F = np.concatenate([0.5 * (b - a) * x + 0.5 * (a + b) for a, b in zip(edges[:-1], edges[1:])])
    W = np.concatenate([0.5 * (b - a) * wq for a, b in zip(edges[:-1], edges[1:])])
    rho = np.exp(beta * (np.cos(F) - 1)) / (2 * np.pi * ive(0, beta))
    inv = np.exp(-beta * (np.cos(F) - 1)) * (2 * np.pi * ive(0, beta))
    lam = 2 * np.pi / (W @ inv)
    var_gp = 2 * np.pi * lam * (1 - 2 * np.pi * lam)
    # g(F) = (F + pi) - lam * int_{-pi}^F 1/rho : cumulative Gauss sums per panel + partial panel
    G = np.empty_like(F)
    acc = 0.0
    for p, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        sl = slice(p * nq, (p + 1) * nq)
        for i, Fi in enumerate(F[sl]):
            xs = 0.5 * (Fi - a) * x + 0.5 * (Fi + a)
            part = 0.5 * (Fi - a) * (wq @ (np.exp(-beta * (np.cos(xs) - 1)) * 2 * np.pi * ive(0, beta)))
            G[sl][i] = (Fi + np.pi) - lam * (acc + part)
        acc += W[sl] @ inv[sl]
    Eg = W @ (rho * G)
    var_g = W @ (rho * (G - Eg) ** 2)
    i0 = ive(0, beta) * np.exp(beta)
    return var_gp / var_g, beta / i0 ** 2, 2 * np.pi * beta ** 2 * np.exp(-2 * beta), 1 - 2 * np.pi * lam


# ---------------- C
def besseli_dec(nu, x):
    x = Dm(x); t = (x / 2) ** nu
    f = Dm(1)
    for k in range(1, nu + 1):
        f *= k
    t = t / f
    s = t; k = 0
    while True:
        k += 1
        t = t * (x / 2) ** 2 / (k * (k + nu))
        s += t
        if k > 10 and t < s * Dm(10) ** (-85):
            return s


def su2_ratios(beta):
    if beta == 0:
        return Dm(0), Dm(0)
    I1 = besseli_dec(1, beta)
    return besseli_dec(2, beta) / I1, besseli_dec(3, beta) / I1


def su2_ratios_quad(beta, n=200000):
    t = (np.arange(n) + 0.5) * np.pi / n
    w = np.sin(t) ** 2 * np.exp(beta * (np.cos(t) - 1))
    return (w @ np.cos(t)) / w.sum(), (w @ (4 * np.cos(t) ** 2 - 1)) / w.sum() / 3


def su2_nested(beta, smax):
    a, c1 = su2_ratios(beta)
    sides = list(range(1, smax + 1, 2)); n = [s * s for s in sides]; N = len(sides)
    E2 = [(1 + 3 * c1 ** k) / 4 for k in n]
    D = [4 * s * (1 - e2) for s, e2 in zip(sides, E2)]
    Cm = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            k, m = min(i, j), max(i, j)
            Cm[i][j] = (a ** (n[m] - n[k]) if m != k else Dm(1)) * E2[k] - a ** (n[k] + n[m])
    Cf = np.array([[float(Cm[i][j]) for j in range(N)] for i in range(N)])
    Dd = np.array([float(x) for x in D])
    sd = 1 / np.sqrt(Dd)
    _, V = np.linalg.eigh(sd[:, None] * Cf * sd[None, :])
    v = [Dm(float(z)) for z in sd * V[:, -1]]
    num = sum(v[i] ** 2 * D[i] for i in range(N))
    den = sum(v[i] * v[j] * Cm[i][j] for i in range(N) for j in range(N))
    single = min(D[i] / Cm[i][i] for i in range(N))
    return float(num / den / 4), float(single / 4), float(den)


if __name__ == "__main__":
    secs = sys.argv[1] if len(sys.argv) > 1 else "ABC"
    if "A" in secs:
        print("# A single plaquette SU(2): gap(H_FF) = lambda_1; claim 3 at beta=0, ~2beta large", flush=True)
        for beta in [0, 1, 5, 10, 50, 200]:
            print(f"beta={beta:6g}  gap={plaquette_gap(beta):.5f}  2beta={2*beta}", flush=True)
    if "B" in secs:
        print("# B U(1) Kramers bound, L->inf: Var(g')/Var(g); compare beta/I0^2, 2 pi beta^2 e^-2beta", flush=True)
        for beta in [2, 4, 6, 10, 15, 20]:
            r, b1, b2, egp = u1_bound_Linf(beta)
            print(f"beta={beta:3d}  bound={r:.4e}  beta/I0^2={b1:.4e}  2pi b^2 e^-2b={b2:.4e}  E g'={egp:.4f}",
                  flush=True)
    if "C" in secs:
        print("# C SU(2) nested loops; moments Bessel series (Decimal) vs midpoint quadrature", flush=True)
        for beta in [2, 10, 100]:
            a, c1 = su2_ratios(beta); aq, cq = su2_ratios_quad(beta)
            print(f"beta={beta}: a={float(a):.12f} |a-a_quad|={abs(float(a)-aq):.1e} "
                  f"|c1-c1_quad|={abs(float(c1)-cq):.1e}", flush=True)
        for beta in [0, 2, 10, 20, 100, 500, 1000]:
            for smax in ([1] if beta == 0 else [41, 201]):
                rq, single, den = su2_nested(beta, smax)
                print(f"beta={beta:5g} smax={smax:3d}  nested RQ (explicit vector, Decimal)={rq:.4f}  "
                      f"best single={single:.4f}  vCv={den:.3e}", flush=True)
