"""Independent single-plaquette checks (no import from src/).

Gauge-invariant states of one open plaquette = class functions of U_p, basis chi_j (orthonormal, Haar).
4 links each carry j, so H_jj = (g^2/2)*4 j(j+1); Tr U chi_j = chi_{j-1/2} + chi_{j+1/2} gives
off-diagonal -(2/g^2)*1.
1. exact Rayleigh-Schrodinger PT to 4th order in x=2/g^2 (sympy, rational)
2. dense diagonalization at large j_max
3. finite-difference solve of the theta equation as a second, basis-free check
4. weak-coupling: harmonic limit 4 and first anharmonic coefficient
"""
import numpy as np
from fractions import Fraction as Fr
from scipy.linalg import eigh_tridiagonal

G2 = [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]


def spin_gap(g2, n=200):
    j = np.arange(n) / 2
    d = 2 * g2 * j * (j + 1)
    e = -(2 / g2) * np.ones(n - 1)
    w = eigh_tridiagonal(d, e, select='i', select_range=(0, 1))[0]
    return w[1] - w[0], w


def fd_gap(g2, N):
    # H phi = (g^2/2)(-phi'' - phi) - (4/g^2) cos(theta) phi on (0, pi), Dirichlet
    h = np.pi / (N + 1)
    th = h * np.arange(1, N + 1)
    d = g2 / h**2 - g2 / 2 - (4 / g2) * np.cos(th)
    e = -g2 / (2 * h**2) * np.ones(N - 1)
    w = eigh_tridiagonal(d, e, select='i', select_range=(0, 1))[0]
    return w[1] - w[0]


def pt_series():
    """Exact RS PT (rational arithmetic) for D + lam V, D_kk = 2 (k/2)(k/2+1), V_{k,k+1} = -1.
    H = g^2 (D + lam V) with lam = (2/g^2)/g^2. Returns (e0, e2, e4) per state; odd orders vanish
    because V has zero diagonal and the graph is bipartite."""
    n = 10
    E0 = [Fr(2) * Fr(k, 2) * (Fr(k, 2) + 1) for k in range(n)]
    V = [[Fr(0)] * n for _ in range(n)]
    for k in range(n - 1):
        V[k][k + 1] = V[k + 1][k] = Fr(-1)
    res = {}
    for s in (0, 1):
        o = [k for k in range(n) if k != s]
        R = lambda k: 1 / (E0[s] - E0[k])
        e2 = sum(V[s][k] ** 2 * R(k) for k in o)
        e4 = sum(V[s][k] * V[k][m] * V[m][l] * V[l][s] * R(k) * R(m) * R(l) for k in o for m in o for l in o)             - e2 * sum(V[s][k] ** 2 * R(k) ** 2 for k in o)
        res[s] = (E0[s], e2, e4)
    return res


if __name__ == "__main__":
    res = pt_series()
    for k, (e0, e2, e4) in res.items():
        print(f"state j={k}/2: E/g^2 = {e0} + {e2} lam^2 + {e4} lam^4,  lam = 2/g^4")
    d0 = res[1][0] - res[0][0]; d2 = res[1][1] - res[0][1]; d4 = res[1][2] - res[0][2]
    # g^2 lam^2 = 4/g^6 ; g^2 lam^4 = 16/g^14
    print(f"Gap = {d0} g^2 + {4*d2}/g^6 + {16*d4}/g^14 + O(g^-22)")
    print()
    print(" g2     spin(n=200)      spin(n=81)       FD(Richardson)   PT2            PT4")
    for g2 in G2:
        a, _ = spin_gap(g2, 200)
        b, _ = spin_gap(g2, 81)
        f1, f2 = fd_gap(g2, 4000), fd_gap(g2, 8000)
        f = (4 * f2 - f1) / 3
        pt2 = 1.5 * g2 + float(4 * d2) / g2**3
        pt4 = pt2 + float(16 * d4) / g2**7
        print(f"{g2:5.1f} {a:16.10f} {b:16.10f} {f:16.10f} {pt2:14.8f} {pt4:14.8f}")
    print()
    print("weak coupling: harmonic omega=2 -> Delta=4; quartic term -theta^4/(6g^2), <theta^4>_n=(3/4)l^4(2n^2+2n+1),")
    print("l^2=g^2/2, n=1->3 => Delta = 4 - (5/8) g^2 + O(g^4)")
    for g2 in [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]:
        a, _ = spin_gap(g2, 600)
        print(f"  g2={g2:5.2f} Delta={a:.10f} (4-Delta)/g2={(4-a)/g2:.6f}")
