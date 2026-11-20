"""Validation of the spin-network plaquette matrix elements and conventions.

1. Clebsch-Gordan series D^{1/2}_{ab} D^j_{mn} = sum_J C C D^J at random group elements.
2. Gauge invariance of spin-network wavefunctions under random vertex gauge transformations.
3. Pointwise identity (Tr U_p psi_G)(U) = sum_G' <G'|Tr U_p|G> psi_G'(U) at random U, using a
   truncation one step larger than the states tested, so no truncation error enters.
4. W hermitian; single plaquette W equals the analytic class-function tridiagonal matrix.
Run: python work/20261120-003-su2-ks-ed/validate.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src import ks_su2 as ks  # noqa: E402

rng = np.random.default_rng(1)


def check_cg_series(tjmax=5):
    err = 0.0
    for _ in range(5):
        v = rng.normal(size=3) * 2
        Dh = ks.D_matrix(1, v)
        for tj in range(tjmax + 1):
            Dj = ks.D_matrix(tj, v)
            for ta in (1, -1):
                for tb in (1, -1):
                    lhs = Dh[ks.m_index(1, ta), ks.m_index(1, tb)] * Dj
                    rhs = 0
                    for tJ in (tj - 1, tj + 1):
                        if tJ < 0:
                            continue
                        rhs = rhs + ks.cg_half(tj, tJ, ta).T @ ks.D_matrix(tJ, v) @ ks.cg_half(tj, tJ, tb)
                    err = max(err, np.max(np.abs(lhs - rhs)))
    # conj(D^j_{mn}) = (-1)^{m-n} D^j_{-m,-n}
    v = rng.normal(size=3)
    e2 = 0.0
    for tj in range(tjmax + 1):
        D = ks.D_matrix(tj, v)
        ms = np.arange(tj, -tj - 1, -2)
        ph = np.array([[(-1.0) ** ((int(a) - int(b)) // 2) for b in ms] for a in ms])
        e2 = max(e2, np.max(np.abs(D.conj() - ph * D[::-1, ::-1])))
    return err, e2


def random_gauge(lat, Us):
    gs = [rng.normal(size=3) * 2 for _ in range(lat.nv)]
    out = []
    for l, (s, e) in enumerate(lat.links):
        U = ks.D_matrix(1, Us[l])
        G = ks.D_matrix(1, gs[s]) @ U @ ks.D_matrix(1, gs[e]).conj().T
        out.append(su2_to_vec(G))
    return out


def su2_to_vec(U):
    # U = exp(-i v.sigma/2) = cos(|v|/2) - i sin(|v|/2) n.sigma
    c = np.clip(np.real(np.trace(U)) / 2, -1, 1)
    th = 2 * np.arccos(c)
    s = np.sin(th / 2)
    if abs(s) < 1e-14:
        return np.zeros(3)
    nx = -np.imag(U[0, 1] + U[1, 0]) / (2 * s)
    ny = np.real(U[1, 0] - U[0, 1]) / (2 * s)
    nz = -np.imag(U[0, 0] - U[1, 1]) / (2 * s)
    return th * np.array([nx, ny, nz])


def plaq_trace(lat, p, Us):
    M = np.eye(2, dtype=complex)
    for l, d in lat.plaqs[p]:
        U = ks.D_matrix(1, Us[l])
        M = M @ (U if d > 0 else U.conj().T)
    return np.trace(M).real


def check_lattice(lat, tjmax_states, n_states=25, n_points=3):
    basis, cas, W = ks.build(lat, tjmax_states + 1)
    herm = abs(W - W.conj().T).max()
    idx = {st: i for i, st in enumerate(basis)}
    small = [i for i, st in enumerate(basis) if max(st[0]) <= tjmax_states]
    pick = rng.choice(small, size=min(n_states, len(small)), replace=False)
    Wc = W.tocsc()
    gerr, perr = 0.0, 0.0
    for _ in range(n_points):
        Us = [rng.normal(size=3) * 1.5 for _ in lat.links]
        Ug = random_gauge(lat, Us)
        for i in pick:
            st = basis[i]
            psi = ks.wavefunction(lat, st, Us)
            gerr = max(gerr, abs(psi - ks.wavefunction(lat, st, Ug)))
            for p in range(len(lat.plaqs)):
                # W column i holds sum_p; rebuild single-plaquette column from moves
                col = {}
                for nst, a in ks.plaquette_moves(lat, st, p, tjmax_states + 1):
                    col[nst] = col.get(nst, 0) + a
                rhs = sum(a * ks.wavefunction(lat, nst, Us) for nst, a in col.items())
                perr = max(perr, abs(plaq_trace(lat, p, Us) * psi - rhs))
    return len(basis), herm, gerr, perr


def check_single_plaquette(tjmax=20):
    lat = ks.single_plaquette()
    basis, cas, W = ks.build(lat, tjmax)
    order = np.argsort([st[0][0] for st in basis])
    Wd = W.toarray()[np.ix_(order, order)]
    # analytic: Tr U chi_j = chi_{j-1/2} + chi_{j+1/2}, off-diagonal magnitude 1
    n = len(basis)
    A = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
    return np.max(np.abs(np.abs(Wd) - A)), np.max(np.abs(cas[order] - [4 * t * (t + 2) / 4 for t in range(n)]))


if __name__ == "__main__":
    e1, e2 = check_cg_series()
    print(f"CG series max error {e1:.2e}; conj(D) identity max error {e2:.2e}")
    a, b = check_single_plaquette()
    print(f"single plaquette: |W| vs analytic tridiagonal max dev {a:.2e}; Casimir dev {b:.2e}")
    for lat, tjm in [(ks.single_plaquette(), 3), (ks.ladder2(), 2), (ks.periodic(2, 2), 1)]:
        n, h, g, p = check_lattice(lat, tjm)
        print(f"{lat.name}: dim(j_max={(tjm + 1) / 2}) = {n}, |W-W^dag| = {h:.2e}, "
              f"gauge-invariance err = {g:.2e}, pointwise plaquette err = {p:.2e}")
