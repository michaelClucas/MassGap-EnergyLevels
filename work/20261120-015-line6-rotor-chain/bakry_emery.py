"""RK / Langevin parent Hamiltonian of the O(3) rotor chain and its Bakry-Emery bound.
Setting: 1D spatial lattice of rotors n_i in S^2 (unit sphere, Ric = 1), finite ring or open chain, any N.
S = -sum_b n_i.n_j,  pi_beta ~ exp(-beta S),  H_FF = sum_i Q_i^* Q_i,  Q_i = grad_i + (beta/2) grad_i S.
gap(H_FF) = Poincare constant of pi_beta for E(f) = sum_i int |grad_i f|^2.
(a) min over configurations of the per-bond Hessian ratio  [ (n1.n2)(|u1|^2+|u2|^2) - 2 u1.u2 ] / (|u1|^2+|u2|^2)
    (analytic value -2, attained at n2 = -n1, u2 = u1): checked by random search.
(b) Bakry-Emery: Ric + beta Hess S >= 1 - 4 beta (every site in 2 bonds), so gap >= 1 - 4 beta for beta < 1/4.
(c) Upper bound from f = sum_i n_i^z on the infinite chain: gap <= (2/3) / chi,
    chi = (1/3) (1+u)/(1-u),  u = coth(beta) - 1/beta  (exact 1D classical Heisenberg chain correlations).
Writes bakry_emery.csv."""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1)


def tangent(n, v):
    return v - (v @ n) * n


def hess_ratio(n1, n2, u1, u2):
    return ((n1 @ n2) * (u1 @ u1 + u2 @ u2) - 2 * u1 @ u2) / (u1 @ u1 + u2 @ u2)


def main():
    best = 0
    for _ in range(200000):
        n1 = rng.normal(size=3); n1 /= np.linalg.norm(n1)
        n2 = rng.normal(size=3); n2 /= np.linalg.norm(n2)
        u1 = tangent(n1, rng.normal(size=3)); u2 = tangent(n2, rng.normal(size=3))
        best = min(best, hess_ratio(n1, n2, u1, u2))
    n1 = np.array([0, 0, 1.0]); u = np.array([1.0, 0, 0])
    print("random-search min per-bond Hessian ratio:", best, " at antiparallel config:", hess_ratio(n1, -n1, u, u))
    rows = []
    for beta in [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 1, 2, 5, 10, 20, 50, 100]:
        u = 1 / np.tanh(beta) - 1 / beta
        chi = (1 + u) / (1 - u) / 3
        ub = (2 / 3) / chi
        lb = max(0.0, 1 - 4 * beta)
        xi = -1 / np.log(u)
        g = 1 / np.sqrt(beta)
        rows.append([beta, f"{g:.4f}", f"{lb:.4f}", f"{ub:.6g}", f"{xi:.4g}"])
        print(*rows[-1])
    with open(os.path.join(HERE, "bakry_emery.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["beta", "g_matching_beta=1/g^2", "BE_lower", "upper_f=sum_nz", "classical_xi"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
