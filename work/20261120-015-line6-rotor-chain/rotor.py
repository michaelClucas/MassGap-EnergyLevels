"""O(3) quantum rotor chain: ED in the angular-momentum basis + certified Feshbach/Schur lower bounds.

Setting: 1+1D Hamiltonian lattice (spatial chain, continuous time), a = 1, finite chain of L sites,
open ('open') or periodic ('per') boundary, fixed coupling g:

    H = gk * sum_i L_i^2 - gv * sum_<ij> n_i . n_j ,   gk = g/2, gv = 1/g  (the target model)

Single-site basis |l, m>, l <= lmax. All work is in the total-M = 0 sector: H conserves M, every
J multiplet has an M = 0 member, so lambda_0, lambda_1 of H|_{M=0} are those of H.

Bounds on the untruncated H (proof in README.md, section 'Certified ED'):
  upper  lambda_k(H) <= lambda_k(A), A = PHP, P = all sites l <= lmax          (min-max)
  lower  lambda_k(H) >= lambda_k(A - X/(d - U_k)) if U_k = lambda_k(A) < d     (Schur complement)
  X = 2 gv^2 sum_b K_b,  K_b = P_b h_b (1 - P_b) h_b P_b  (two-site, exact),  h_b = n_i . n_j
  d = max_alpha [ alpha gk (lmax+1)(lmax+2) + LB E0(H with gk -> (1-alpha) gk) ]
  LB E0 from the Anderson bond bound (sum of certified two-rotor ground energies) or, one level
  deep, from this same certified ED.
Rounding: eigsh residual norm + 10 n eps ||A||_1 margin. Lanczos completeness (no eigenvalue missed
below the computed ones) is assumed, not proved: 'numerical-certified'.
"""
from __future__ import annotations

import functools

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ---------------------------------------------------------------- single site
@functools.lru_cache(maxsize=None)
def site_ops(lmax):
    lm = [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]
    idx = {s: i for i, s in enumerate(lm)}
    d = len(lm)
    Z = np.zeros((d, d)); Pp = np.zeros((d, d))
    for (l, m), i in idx.items():
        if (l + 1, m) in idx:
            v = np.sqrt(((l + 1) ** 2 - m ** 2) / ((2 * l + 1) * (2 * l + 3)))
            Z[idx[(l + 1, m)], i] = v; Z[i, idx[(l + 1, m)]] = v
        if (l + 1, m + 1) in idx:  # n+ = sin(theta) e^{i phi}, Condon-Shortley
            Pp[idx[(l + 1, m + 1)], i] = np.sqrt((l + m + 1) * (l + m + 2) / ((2 * l + 1) * (2 * l + 3)))
        if (l - 1, m + 1) in idx:
            Pp[idx[(l - 1, m + 1)], i] = -np.sqrt((l - m) * (l - m - 1) / ((2 * l - 1) * (2 * l + 1)))
    Lsq = np.array([l * (l + 1) for l, m in lm], float)
    mm = np.array([m for l, m in lm]); ll = np.array([l for l, m in lm])
    return dict(d=d, Z=Z, Pp=Pp, Pm=Pp.T.copy(), Lsq=Lsq, m=mm, l=ll)


def check_site(lmax=6):
    """x^2+y^2+z^2 = 1 on the block l <= lmax-1 (validates the matrix elements)."""
    o = site_ops(lmax)
    S = o["Z"] @ o["Z"] + 0.5 * (o["Pp"] @ o["Pm"] + o["Pm"] @ o["Pp"])
    keep = o["l"] <= lmax - 1
    return float(np.abs(S[np.ix_(keep, keep)] - np.eye(keep.sum())).max())


@functools.lru_cache(maxsize=None)
def bond_ops(lmax):
    """h = n_1 . n_2 on two sites (sparse d^2 x d^2) and K = P h (1-P) h P with the exact shell (lmax+1)."""
    def nn(o):
        Z, Pp, Pm = (sp.csr_matrix(o[k]) for k in ("Z", "Pp", "Pm"))
        return (sp.kron(Z, Z) + 0.5 * (sp.kron(Pp, Pm) + sp.kron(Pm, Pp))).tocsr()
    o = site_ops(lmax); d = o["d"]
    ob = site_ops(lmax + 1); db = ob["d"]
    h, hb = nn(o), nn(ob)
    inP = ((ob["l"][:, None] <= lmax) & (ob["l"][None, :] <= lmax)).ravel()
    pidx = np.array([a * db + b for a in range(d) for b in range(d)])
    assert inP[pidx].all() and inP.sum() == d * d
    out = np.nonzero(~inP)[0]
    K = (hb[pidx][:, out] @ hb[out][:, pidx]).tocsr()
    assert abs(hb[pidx][:, pidx] - h).max() < 1e-14
    return h, K


# ---------------------------------------------------------------- chain
def bonds(L, bc):
    b = [(i, i + 1) for i in range(L - 1)]
    if bc == "per" and L >= 3:
        b.append((L - 1, 0))
    return b


def _embed2(op2, i, j, L, d):
    """Two-site operator on sites i<j (adjacent or (0, L-1)) in the full d^L space, sparse."""
    if j == i + 1:
        return sp.kron(sp.kron(sp.identity(d ** i, format="csr"), op2 if sp.issparse(op2) else sp.csr_matrix(op2)),
                       sp.identity(d ** (L - j - 1), format="csr"), format="csr")
    # wrap bond (L-1, 0): op2 acts on (site L-1, site 0); write as sum of products
    op2 = op2.toarray() if sp.issparse(op2) else op2
    out = None
    U, s, Vt = np.linalg.svd(op2.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d))
    for k in range(len(s)):
        if s[k] < 1e-14:
            continue
        Xa = (U[:, k] * s[k]).reshape(d, d)  # on site L-1
        Xb = Vt[k].reshape(d, d)              # on site 0
        t = sp.kron(sp.kron(sp.csr_matrix(Xb), sp.identity(d ** (L - 2), format="csr")),
                    sp.csr_matrix(Xa), format="csr")
        out = t if out is None else out + t
    return out


@functools.lru_cache(maxsize=64)
def m_sector(L, lmax, M=0):
    o = site_ops(lmax)
    d = o["d"]
    tot = np.zeros(1, dtype=np.int64)
    for _ in range(L):
        tot = (tot[:, None] + o["m"][None, :]).ravel()
    return np.nonzero(tot == M)[0]


def chain_parts(L, bc, lmax, M=0):
    """Returns (T, V, X) restricted to the M sector: T = sum L^2, V = sum n.n, X = sum_b K_b."""
    o = site_ops(lmax)
    d = o["d"]
    idx = m_sector(L, lmax, M)
    h, K = bond_ops(lmax)
    diag = np.zeros(1)
    for _ in range(L):
        diag = (diag[:, None] + o["Lsq"][None, :]).ravel()
    T = sp.diags(diag[idx]).tocsr()
    V = None; X = None
    for (i, j) in bonds(L, bc):
        a, b = (i, j) if j == i + 1 else (i, j)
        hv = _embed2(h, a, b, L, d)[idx][:, idx]
        kv = _embed2(K, a, b, L, d)[idx][:, idx]
        V = hv if V is None else V + hv
        X = kv if X is None else X + kv
    if V is None:
        V = sp.csr_matrix(T.shape); X = sp.csr_matrix(T.shape)
    return T, V.tocsr(), X.tocsr()


def lowest(Hm, k=2):
    n = Hm.shape[0]
    if n <= 1500:
        w, v = np.linalg.eigh(Hm.toarray() if sp.issparse(Hm) else Hm)
        w, v = w[:k], v[:, :k]
    else:
        w, v = spla.eigsh(Hm, k=k, which="SA", tol=1e-12, ncv=max(20, 3 * k + 5), maxiter=20000)
        o = np.argsort(w); w, v = w[o], v[:, o]
    res = np.linalg.norm(Hm @ v - v * w, axis=0)
    nrm = abs(Hm).sum(axis=0).max()
    margin = res + 10 * n * np.finfo(float).eps * nrm
    return w, margin


# ---------------------------------------------------------------- two-rotor (Anderson bond) bound
@functools.lru_cache(maxsize=None)
def two_rotor_lower(a, b, gv, l2=10):
    """Certified lower bound on inf spec(a L1^2 + b L2^2 - gv n1.n2) (M = 0 sector suffices)."""
    o = site_ops(l2); d = o["d"]
    idx = m_sector(2, l2, 0)
    h, K = bond_ops(l2)
    D = (np.kron(a * o["Lsq"], np.ones(d)) + np.kron(np.ones(d), b * o["Lsq"]))[idx]
    A = (np.diag(D) - gv * h[idx][:, idx].toarray())
    Xm = (gv ** 2) * K[idx][:, idx].toarray()
    floor = min(a, b) * (l2 + 1) * (l2 + 2) - gv
    U = np.linalg.eigvalsh(A)[0] + 1e-12
    if U >= floor:
        return -np.inf
    lo = np.linalg.eigvalsh(A - Xm / (floor - U))[0]
    return lo - 10 * A.shape[0] * np.finfo(float).eps * np.abs(A).sum(1).max()


def anderson_E0_lower(L, bc, gk, gv):
    bl = bonds(L, bc)
    deg = np.zeros(L)
    for i, j in bl:
        deg[i] += 1; deg[j] += 1
    if not bl:
        return 0.0
    # kinetic >= 0 and n.n <= 1 give the trivial bond bound -gv
    return sum(max(-gv, two_rotor_lower(round(gk / deg[i], 12), round(gk / deg[j], 12), gv)) for i, j in bl)


# ---------------------------------------------------------------- certified bounds
_parts_cache = {}


def parts(L, bc, lmax):
    key = (L, bc, lmax)
    if key not in _parts_cache:
        if len(_parts_cache) > 4:
            _parts_cache.clear()
        _parts_cache[key] = chain_parts(L, bc, lmax)
    return _parts_cache[key]


def certify(L, bc, gk, gv, lmax, depth=1, alphas=(0.3, 0.6, 0.85), want_k=2):
    """Returns dict with U0, U1 (upper), L0, L1 (lower, -inf if none), floor d, dim."""
    T, V, X = parts(L, bc, lmax)
    A = gk * T - gv * V
    nb = 2.0
    Xs = nb * gv ** 2 * X
    w, mar = lowest(A, k=want_k)
    U = w + mar
    c = (lmax + 1) * (lmax + 2)
    best = -np.inf
    for al in alphas:
        lb = anderson_E0_lower(L, bc, gk * (1 - al), gv)
        if depth > 0:
            sub = certify(L, bc, gk * (1 - al), gv, lmax, depth - 1, alphas, want_k=1)
            lb = max(lb, sub["L0"])
        best = max(best, al * gk * c + lb)
    d = best
    lows = []
    for k in range(want_k):
        if U[k] < d:
            wk, mk = lowest(A - Xs / (d - U[k]), k=k + 1)
            lows.append(wk[k] - mk[k])
        else:
            lows.append(-np.inf)
    out = dict(dim=A.shape[0], floor=d, U0=U[0], L0=lows[0])
    if want_k > 1:
        out.update(U1=U[1], L1=lows[1])
    return out
