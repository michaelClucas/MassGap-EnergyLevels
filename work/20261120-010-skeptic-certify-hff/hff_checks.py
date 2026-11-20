"""Skeptic checks for C6 (H_FF = sum_l Q_l^* Q_l, psi = exp(-beta S_W / 2), SU(2)).

1. Single plaquette, gauge-invariant sector (class functions of U_p, Peter-Weyl basis chi_j):
   H_FF = 4 E^2 - (3 beta/2) cos(theta) + (beta^2/4) sin^2(theta)   (from the 007 explicit form)
   must have lowest eigenvalue exactly 0 (frustration-free) for every beta; this tests the
   coefficients Delta cos = -3 cos, |grad cos|^2 = sin^2 and E^2 = -Delta/4 together.
   Also: gauge-invariant gap >= Bakry-Emery bound 1/2 - beta (one plaquette per link).
2. Hessian of the Wilson action S = sum_p (1 - 1/2 tr U_p) in normal (exponential) coordinates
   U_l exp(sum_a t_a i sigma_a) on periodic lattices (d=2, L=3 and L=2; d=3, L=2): lambda_min(Hess S)
   must be >= -8(d-1) (007 bound). Random configurations plus the pi-flux configuration (all U_p = -1).
"""
import itertools

import numpy as np

rng = np.random.default_rng(1)
SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
       np.array([[1, 0], [0, -1]], complex)]
GEN = [1j * s for s in SIG]  # unit-speed generators: |i sigma_a|_op = 1, orthonormal on unit S^3


# ------------------------------------------------------------------ 1. single plaquette
def plaquette_hff(beta, nj=121):
    j = np.arange(nj) / 2
    C = np.zeros((nj, nj))  # multiplication by cos(theta) = 1/2 tr U = (1/2)(chi_{j+1/2}+chi_{j-1/2})
    for i in range(nj - 1):
        C[i, i + 1] = C[i + 1, i] = 0.5
    H = np.diag(4 * j * (j + 1)) - 1.5 * beta * C + 0.25 * beta**2 * (np.eye(nj) - C @ C)
    # truncation: C @ C is wrong only in the last row/col; nj large and the ground state is
    # concentrated at small j, so the top corner is irrelevant to E0, E1 at the printed accuracy.
    ev = np.linalg.eigvalsh(H)
    return ev[0], ev[1] - ev[0]


# ------------------------------------------------------------------ 2. Hessian of S
def expm_su2(t):
    """exp(sum_a t_a i sigma_a) = cos|t| + i sin|t| (t.sigma)/|t|."""
    n = np.linalg.norm(t)
    if n < 1e-15:
        return np.eye(2, dtype=complex)
    return np.cos(n) * np.eye(2) + 1j * np.sin(n) * sum(t[a] * SIG[a] for a in range(3)) / n


def lattice(d, L):
    sites = list(itertools.product(range(L), repeat=d))
    idx = {s: i for i, s in enumerate(sites)}
    link = lambda s, mu: idx[s] * d + mu
    plaqs = []
    for s in sites:
        for mu in range(d):
            for nu in range(mu + 1, d):
                spm = tuple((s[k] + (k == mu)) % L for k in range(d))
                spn = tuple((s[k] + (k == nu)) % L for k in range(d))
                plaqs.append((link(s, mu), link(spm, nu), link(spn, mu), link(s, nu)))
    return len(sites) * d, plaqs


def action(U, plaqs):
    S = 0.0
    for a, b, c, e in plaqs:
        P = U[a] @ U[b] @ U[c].conj().T @ U[e].conj().T
        S += 1 - 0.5 * np.trace(P).real
    return S


def hessian(U, plaqs, h=1e-3):
    nl = len(U)
    n = 3 * nl

    def S_at(t):
        V = [U[l] @ expm_su2(t[3 * l:3 * l + 3]) for l in range(nl)]
        return action(V, plaqs)

    Hm = np.zeros((n, n))
    s0 = S_at(np.zeros(n))
    for i in range(n):
        for k in range(i, n):
            if i == k:
                e = np.zeros(n); e[i] = h
                Hm[i, i] = (S_at(e) - 2 * s0 + S_at(-e)) / h**2
            else:
                if i // 3 != k // 3 and not any(i // 3 in p and k // 3 in p for p in plaqs):
                    continue  # links not sharing a plaquette: mixed derivative is 0
                e = np.zeros(n); e[i] = h; e[k] = h
                f = np.zeros(n); f[i] = h; f[k] = -h
                Hm[i, k] = Hm[k, i] = (S_at(e) - S_at(f) - S_at(-f) + S_at(-e)) / (4 * h * h)
    return Hm


def random_su2():
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    return q[0] * np.eye(2) + 1j * (q[1] * SIG[0] + q[2] * SIG[1] + q[3] * SIG[2])


def pi_flux(d, L):
    """U_mu(x) = exp(i sigma_3 pi sum_{nu<mu} x_nu): every plaquette holonomy -1 (L even)."""
    sites = list(itertools.product(range(L), repeat=d))
    U = []
    for s in sites:
        for mu in range(d):
            th = np.pi * sum(s[:mu])
            U.append(np.cos(th) * np.eye(2) + 1j * np.sin(th) * SIG[2])
    return U


def main():
    print("== 1. single plaquette H_FF (gauge-invariant sector)")
    for beta in [0.0, 0.1, 0.5, 1.0, 3.0, 10.0]:
        e0, gap = plaquette_hff(beta)
        print(f"  beta={beta:5}: E0 = {e0:+.2e} (must be 0), gauge-inv gap = {gap:.6f}, "
              f"BE bound 1/2-beta = {0.5-beta:+.3f}  ok={gap >= 0.5 - beta - 1e-12}")
    print("== 2. lambda_min(Hess S) vs -8(d-1)")
    for d, L in [(2, 2), (2, 3), (3, 2)]:
        nl, plaqs = lattice(d, L)
        cnt = np.zeros(nl, int)
        for p in plaqs:
            for l in p:
                cnt[l] += 1
        worst = np.inf
        for trial in range(4):
            U = [random_su2() for _ in range(nl)]
            worst = min(worst, np.linalg.eigvalsh(hessian(U, plaqs))[0])
        res = f"random min {worst:.4f}"
        if L % 2 == 0:
            U = pi_flux(d, L)
            P = U[plaqs[0][0]] @ U[plaqs[0][1]] @ U[plaqs[0][2]].conj().T @ U[plaqs[0][3]].conj().T
            ev = np.linalg.eigvalsh(hessian(U, plaqs))
            res += f"; pi-flux (1/2 tr U_p = {0.5*np.trace(P).real:+.0f}) min {ev[0]:.4f}"
        print(f"  d={d} L={L}: links {nl}, plaquettes/link {set(cnt.tolist())}, bound "
              f"-8(d-1) = {-8*(d-1)}: {res}")


if __name__ == "__main__":
    main()
