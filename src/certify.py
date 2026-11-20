"""Certified two-sided bounds on low eigenvalues of the UNTRUNCATED SU(2) Kogut-Susskind Hamiltonian
on a finite lattice at fixed coupling, from the truncated spin-network matrices of `src/ks_su2.py`.

Setting: 2+1D Hamiltonian lattice, SU(2), a = 1, finite lattice, fixed g^2. The Hilbert space is a
gauge-invariant flux sector S (defined by `sector_states`), split as S = P (+) Q with
P = span of spin networks in S with every link spin <= j_max and Q its orthogonal complement.

Bounds (statement and proof in work/20261120-006-certified-gap-bounds/README.md):
  upper:  lambda_k(H|S) <= lambda_k(A),  A = P H P                       (min-max / Rayleigh-Ritz)
  lower:  lambda_k(H|S) >= L_k := lambda_k(A - X / (d - lambda_k(A)))    (Schur-complement bound)
          where d <= inf spec(Q H Q) and X >= P H Q H P as forms on P, provided lambda_k(A) < d.

Two choices of X:
  'exact' : X = (4/g^4) W_{P,shell} W_{shell,P}, shell = spin networks in S with max spin j_max+1/2
            (Q W P lands only there because Tr U_p changes each loop spin by +-1/2).
  'crude' : X = (4/g^4) * 4 N_plaq * diag(n_i), n_i = # plaquettes of basis state i containing a link
            at spin j_max (Cauchy-Schwarz with ||Tr U_p|| <= 2; needs no shell basis).
d = (g^2/2) c_Q - 4 N_plaq / g^2, with c_Q a proved lower bound on sum_l j_l(j_l+1) over gauge-invariant
spin configurations with some spin > j_max (`casimir_floor`).

Rounding: dense LAPACK eigenvalues are backward stable; a margin of 10 n eps ||A||_inf is subtracted
from lower and added to upper bounds. This is not interval arithmetic.
"""
from __future__ import annotations

import itertools
from collections import deque

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

from src import ks_su2 as ks


# ---------------------------------------------------------------- sectors
def flux_cuts(lat):
    """Z2 electric-flux cuts. Open plaquette/ladder: none (every Z2 cycle is a sum of plaquette
    boundaries). periodicLxLy: x-links at x=0 (all y) and y-links at y=0 (all x); the parity of
    sum 2j over a cut is conserved by every Tr U_p and by Gauss's law independent of the cut chosen."""
    if not lat.name.startswith("periodic"):
        return []
    cx, cy = [], []
    for l, (s, e) in enumerate(lat.links):
        # periodic(): link 2*(x+Lx*y) is the x-link at (x,y), 2*(...)+1 the y-link
        site = l // 2
        if l % 2 == 0 and site % _Lx(lat) == 0:
            cx.append(l)
        if l % 2 == 1 and site // _Lx(lat) == 0:
            cy.append(l)
    return [cx, cy]


def _Lx(lat):
    return int(lat.name[len("periodic"):].split("x")[0])


def _vertex_ok(tjs_at_v):
    """SU(2) invariant exists iff sum of 2j even and 2*max <= sum (iterated triangle rule)."""
    if not tjs_at_v:
        return True
    s = sum(tjs_at_v)
    return s % 2 == 0 and 2 * max(tjs_at_v) <= s


def spin_configs(lat, tjmax):
    """All link-spin configurations (tj = 2j <= tjmax) admitting an invariant at every vertex.
    Depth-first, links ordered so vertices close early; prunes with the necessary condition
    2*max(assigned) <= sum(assigned) + (#unassigned) * tjmax."""
    L = len(lat.links)
    vlinks = [[l for l, _ in lat.legs[v]] for v in range(lat.nv)]
    order, done = [], set()
    while len(order) < L:  # greedy: pick the link that leaves the fewest open legs at its ends
        best = min((l for l in range(L) if l not in done),
                   key=lambda l: sum(len([m for m in vlinks[v] if m not in done and m != l])
                                     for v in set(lat.links[l])))
        order.append(best); done.add(best)
    pos = {l: i for i, l in enumerate(order)}
    tjs = [None] * L
    out = []

    def check(depth):
        l = order[depth]
        for v in set(lat.links[l]):
            assigned = [tjs[m] for m in vlinks[v] if pos[m] <= depth]
            r = len(vlinks[v]) - len(assigned)
            if r == 0:
                if not _vertex_ok(assigned):
                    return False
            elif 2 * max(assigned) > sum(assigned) + r * tjmax:
                return False
        return True

    def rec(depth):
        if depth == L:
            out.append(tuple(tjs)); return
        for t in range(tjmax + 1):
            tjs[order[depth]] = t
            if check(depth):
                rec(depth + 1)
        tjs[order[depth]] = None

    rec(0)
    return out


def sector_states(lat, tjmax, flux=None):
    """Orthonormal spin-network states (tjs, ints) of the gauge-invariant sector with Z2 flux
    parities `flux` (tuple, default all zero), all spins <= tjmax/2."""
    cuts = flux_cuts(lat)
    flux = tuple(flux) if flux is not None else tuple([0] * len(cuts))
    states = []
    for tjs in spin_configs(lat, tjmax):
        if tuple(sum(tjs[l] for l in c) % 2 for c in cuts) != flux:
            continue
        nint = [ks.intertwiners(lat.vertex_legs(v, tjs)).shape[0] for v in range(lat.nv)]
        if min(nint) == 0:
            continue
        for ints in itertools.product(*[range(n) for n in nint]):
            states.append((tjs, tuple(ints)))
    return states


def build_sector(lat, tjmax, flux=None):
    """Basis = full sector at spins <= tjmax/2 (not just the part reachable from the vacuum),
    Casimir sums, and W = sum_p Tr U_p on it. Raises if Tr U_p leaves the enumerated set (that
    would mean the enumeration is incomplete)."""
    basis = sector_states(lat, tjmax, flux)
    index = {st: i for i, st in enumerate(basis)}
    rows, cols, vals = [], [], []
    for i, st in enumerate(basis):
        for p in range(len(lat.plaqs)):
            for nst, a in ks.plaquette_moves(lat, st, p, tjmax):
                rows.append(index[nst]); cols.append(i); vals.append(a)  # KeyError = bug
    n = len(basis)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    if np.iscomplexobj(W.data):
        assert np.abs(W.data.imag).max(initial=0) < 1e-12, "W not real"
        W = sp.csr_matrix(W.real)
    W.eliminate_zeros()
    cas = np.array([sum(t * (t + 2) / 4 for t in st[0]) for st in basis])
    return basis, cas, W


def vacuum_component_size(basis, W, mask):
    """Number of states in the truncated set `mask` connected to the strong-coupling vacuum by
    W restricted to that set (= dimension of the iteration-1 BFS basis)."""
    idx = np.where(mask)[0]
    Wm = W[idx][:, idx]
    _, lab = connected_components(abs(Wm) > 1e-13, directed=False)
    vac = next(i for i, st in enumerate(basis) if max(st[0]) == 0)
    k = int(np.searchsorted(idx, vac))
    return int(np.sum(lab == lab[k]))


# ---------------------------------------------------------------- electric floor on Q
def casimir_lb_analytic(lat, J):
    """Proved lower bound on sum_l j_l(j_l+1) for any gauge-invariant configuration containing a link
    of spin J. At an endpoint with k other legs, Gauss's law gives sum(other j) >= J, and
    Cauchy-Schwarz gives sum j(j+1) >= (sum j)^2/k + sum j >= J^2/k + J. Both endpoints count when
    their other-leg sets are disjoint; otherwise only the better one."""
    best = np.inf
    for l, (s, e) in enumerate(lat.links):
        os_ = {m for m, _ in lat.legs[s] if m != l}
        oe = {m for m, _ in lat.legs[e] if m != l}
        if s == e or not os_ or not oe:
            continue
        fs = J * J / len(os_) + J
        fe = J * J / len(oe) + J
        extra = fs + fe if os_.isdisjoint(oe) else max(fs, fe)
        best = min(best, J * (J + 1) + extra)
    return best


def casimir_floor(lat, tjmax, budget=2e6):
    """c_Q: lower bound on sum_l j_l(j_l+1) over gauge-invariant configurations (any flux sector)
    with max spin >= (tjmax+1)/2. Exact minimum by vectorized enumeration of all configurations with
    spins <= tjcap, raising tjcap until the analytic bound for max spin > tjcap/2 exceeds the
    minimum found; falls back to the analytic bound if the enumeration exceeds `budget` configs."""
    tj0 = tjmax + 1
    lb0 = casimir_lb_analytic(lat, tj0 / 2)
    L = len(lat.links)
    tjcap = tj0
    while (tjcap + 1) ** L <= budget:
        grids = np.indices([tjcap + 1] * L, dtype=np.int16).reshape(L, -1)
        ok = grids.max(axis=0) >= tj0
        for v in range(lat.nv):
            ls = [l for l, _ in lat.legs[v]]
            sub = grids[ls]
            s = sub.sum(axis=0)
            ok &= (s % 2 == 0) & (2 * sub.max(axis=0) <= s)
        if ok.any():
            c = (grids[:, ok].astype(float) * (grids[:, ok] + 2) / 4).sum(axis=0).min()
            if casimir_lb_analytic(lat, (tjcap + 1) / 2) >= c:
                return float(c), "enumerated"
        tjcap += 1
    return float(lb0), "analytic"


# ---------------------------------------------------------------- the bound
DENSE_MAX = 2500
DENSE_LDL_MAX = 12000


def _eps_margin(M):
    n = M.shape[0]
    rowsum = np.asarray(abs(M).sum(axis=1)).ravel().max()
    return 10 * n * np.finfo(float).eps * max(1.0, rowsum)


def upper_bounds(M, k):
    """lambda_0..lambda_{k-1}(M) from above. Dense: LAPACK. Sparse: Lanczos vectors, then exact
    Rayleigh-Ritz on their span (min-max makes Ritz values upper bounds for any subspace)."""
    if M.shape[0] <= DENSE_MAX:
        Md = M.toarray() if sp.issparse(M) else M
        return sla.eigvalsh(Md, subset_by_index=[0, k - 1], driver="evr")
    _, V = spla.eigsh(M, k=k + 2, which="SA", tol=1e-12)
    V, _ = np.linalg.qr(V)
    return np.sort(np.linalg.eigvalsh(V.T @ (M @ V)))[:k]


def ldl_count(M, E):
    """Sylvester inertia of M - E via sparse LU with diagonal pivoting only (SuperLU, symmetric
    mode). With K' = L D L^T (symmetric by construction, D = diag U) and R = K' - P(M-E)P^T computed
    exactly, In(K') = In(D) and |lambda_k(M-E) - lambda_k(K')| <= ||R||_F (Weyl).
    Returns (#negative pivots, ||R||_F); (None, inf) if SuperLU permuted rows asymmetrically."""
    n = M.shape[0]
    K = (sp.csc_matrix(M) - E * sp.identity(n, format="csc")).tocsc()
    lu = spla.splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                   options=dict(SymmetricMode=True))
    if not np.array_equal(lu.perm_r, lu.perm_c):
        return None, np.inf
    D = lu.U.diagonal()
    Pr = sp.csc_matrix((np.ones(n), (lu.perm_r, np.arange(n))), shape=(n, n))
    Pc = sp.csc_matrix((np.ones(n), (np.arange(n), lu.perm_c)), shape=(n, n))
    R = lu.L @ sp.diags(D) @ lu.L.T - Pr @ K @ Pc
    return int(np.sum(D < 0)), float(sp.linalg.norm(R))


def dense_ldl_count(M, E):
    """Sylvester inertia of M - E from a dense Bunch-Kaufman factorization K = L D L^T (LAPACK
    sytrf via scipy.linalg.ldl; D has 1x1/2x2 blocks). R = L D L^T - K is formed explicitly, so
    In(K + R) = In(D) and Weyl gives |lambda_k(K) - lambda_k(K+R)| <= ||R||_F."""
    K = (M.toarray() if sp.issparse(M) else np.array(M)) - E * np.eye(M.shape[0])
    L, D, _ = sla.ldl(K, lower=True, hermitian=True)
    ev = sla.eigvalsh_tridiagonal(np.diag(D).copy(), np.diag(D, -1).copy())
    R = L @ (D @ L.T) - K
    return int(np.sum(ev < 0)), float(np.linalg.norm(R))


def lower_bound_k(M, k):
    """Proved (up to rounding) lower bound on lambda_k(M) for a finite symmetric matrix M."""
    if M.shape[0] <= DENSE_MAX:
        Md = M.toarray() if sp.issparse(M) else M
        return sla.eigvalsh(Md, subset_by_index=[0, k], driver="evr")[k] - _eps_margin(M)
    est = np.sort(spla.eigsh(M, k=k + 2, which="SA", tol=1e-12, return_eigenvectors=False))[k]
    delta = 1e-9 * max(1.0, abs(est))
    for _ in range(6):
        neg, err = (dense_ldl_count if M.shape[0] <= DENSE_LDL_MAX else ldl_count)(M, est - delta)
        if neg is not None and neg <= k and err < 1e-3 * delta + 1e-12:
            return est - delta - err
        delta *= 10
    return -np.inf


def schur_lower(A, X, d, E, k):
    """lambda_k(H) >= lambda_k(A - X/(d-E)) when E >= lambda_k(H) and E < d (see README)."""
    if not E < d:
        return -np.inf
    return lower_bound_k(sp.csr_matrix(A) - sp.csr_matrix(X) / (d - E), k)


def block_ops(cas, W, g2, rows, n_plaq, at_max_counts=None, shell=None):
    """A = H restricted to `rows`; X >= (P H Q)(Q H P) for the complement Q of `rows` in the sector.
    shell (index array) given: exact X from W[shell, rows]. Otherwise crude diagonal X."""
    A = (sp.diags(0.5 * g2 * cas[rows]) - (2.0 / g2) * W[rows][:, rows]).tocsr()
    if shell is not None:
        Bs = W[shell][:, rows]
        X = (4.0 / g2**2) * (Bs.T @ Bs)
    else:
        X = sp.diags((4.0 / g2**2) * 4.0 * n_plaq * np.asarray(at_max_counts, float))
    return A, sp.csr_matrix(X)


def certify(cas, W, g2, n_plaq, iP, c_Q, kmax=2, shell=None, at_max_P=None,
            iQ1=None, at_top_Q1=None, c_Q_top=None, lower=True):
    """Two-sided bounds on lambda_0..lambda_{kmax-1} of H on the sector.

    iP: basis indices of P (spins <= j_max). shell: indices of spin-(j_max+1/2) states (exact X);
    if None, crude X with at_max_P. Floor d on Q H Q: crude d0 = (g^2/2) c_Q - 4 N_p/g^2, improved
    when iQ1 (all basis states above j_max, the basis being complete up to j_top) is given, by the
    same Schur bound applied to Q = Q1 (+) Q2 with Q2 crude (c_Q_top, at_top_Q1)."""
    A, X = block_ops(cas, W, g2, iP, n_plaq, at_max_P, shell)
    d0 = 0.5 * g2 * c_Q - 4.0 * n_plaq / g2
    d, dsrc = d0, "crude"
    if lower and iQ1 is not None and 0 < len(iQ1) <= DENSE_MAX:
        A1, X1 = block_ops(cas, W, g2, iQ1, n_plaq, at_top_Q1, None)
        d2 = 0.5 * g2 * c_Q_top - 4.0 * n_plaq / g2
        e1 = upper_bounds(A1, 1)[0]
        dn = schur_lower(A1, X1, d2, e1, 0)
        if dn > d:
            d, dsrc = dn, "nested"
    hi = upper_bounds(A, kmax) + _eps_margin(A)
    lo = np.array([schur_lower(A, X, d, hi[k], k) if lower else -np.inf for k in range(kmax)])
    return dict(lo=lo, hi=hi, d=d, dsrc=dsrc, dim=A.shape[0],
                method="schur-exact-shell" if shell is not None else "schur-crude-shell")


def plaq_count_at_max(lat, basis, mask, tjmax):
    """n_i = number of plaquettes of state i containing a link with spin exactly tjmax/2."""
    out = []
    for st, m in zip(basis, mask):
        if m:
            out.append(sum(any(st[0][l] == tjmax for l, _ in p) for p in lat.plaqs))
    return np.array(out)
