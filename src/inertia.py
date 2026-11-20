"""Sylvester inertia counts: number of eigenvalues of a real symmetric matrix below a shift s,
from an LDL^T factorization of A - sI, with an explicit a posteriori residual.

Why a factorization gives the count (Sylvester's law of inertia): if K = A - sI and
K' = P L D L^T P^T with L unit lower triangular (so nonsingular) and P a permutation, then
K' is congruent to D and In(K') = In(D). D is block diagonal with 1x1 and (Bunch-Kaufman) 2x2
blocks; a 2x2 symmetric block has one negative and one positive eigenvalue if det < 0, and two of
the sign of its diagonal if det > 0. So #neg(K') is read off D exactly.
Plain LU with row pivoting (K = P L U, P acting on rows only) is NOT a congruence, so the signs of
diag(U) do not in general give In(K). Two factorizations are used here, both congruences:
  dense : LAPACK sytrf (Bunch-Kaufman, symmetric pivoting) via scipy.linalg.ldl.
  sparse: SuperLU with diag_pivot_thresh = 0 and SymmetricMode (no row pivoting, symmetric fill
          ordering MMD on A^T + A). The routine checks perm_r == perm_c, i.e. the row and column
          permutations agree, so K' = L diag(U) L^T is a symmetric permutation congruence. Without
          pivoting this is not backward stable for indefinite K, which is why the residual is
          computed explicitly instead of assumed.
Floating point: the computed L, D are exact floating-point numbers, so In(L D L^T) = In(D) holds
exactly. R = P L D L^T P^T - K is evaluated in floating point; ||R||_2 <= ||R_computed||_F
+ 2 n eps || |L| |D| |L|^T ||_F (the second term bounds the rounding in forming the product).
Weyl's inequality |lambda_i(K) - lambda_i(K+R)| <= ||R||_2 then gives, with r = that bound,
  #eig(A < s - r) <= #neg(D) <= #eig(A < s + r).
`inertia_below` factors at s + delta and s - delta (delta >= the residual at each), so
  n_lo <= #eig(A < s) <= n_hi,
and the count is determined when n_lo = n_hi. This is a floating-point margin computed from the
actual residual, not interval arithmetic: rounding inside the residual evaluation is bounded by
the standard model, not enclosed.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EPS = np.finfo(float).eps
DENSE_MAX = 5000


def _gamma(m):
    """Standard rounding constant gamma_m = m eps / (1 - m eps) (Higham, Accuracy and Stability, 3.1)."""
    return m * EPS / (1 - m * EPS)


def _neg_from_blocks(D):
    """Number of negative eigenvalues and zero pivots of a block-diagonal D (1x1, 2x2 blocks)."""
    n = D.shape[0]
    neg = zero = 0
    i = 0
    while i < n:
        if i + 1 < n and D[i + 1, i] != 0.0:
            a, b, c = D[i, i], D[i + 1, i], D[i + 1, i + 1]
            det = a * c - b * b
            if det < 0:
                neg += 1
            elif det > 0:
                neg += 2 if a < 0 else 0
            else:
                zero += 1
            i += 2
        else:
            if D[i, i] < 0:
                neg += 1
            elif D[i, i] == 0:
                zero += 1
            i += 1
    return neg, zero


def dense_count(A, s):
    """(#neg, #zero pivots, residual bound r) for A - sI, Bunch-Kaufman."""
    K = (A.toarray() if sp.issparse(A) else np.array(A, dtype=float)) - s * np.eye(A.shape[0])
    lu, D, _ = sla.ldl(K, lower=True, hermitian=True)
    neg, zero = _neg_from_blocks(D)
    d0, d1 = np.diag(D).copy(), np.diag(D, -1).copy()  # D is tridiagonal (1x1 and 2x2 blocks)

    def tri_times(dd0, dd1, M):  # (tridiagonal symmetric) @ M without forming a dense product
        out = dd0[:, None] * M
        out[1:] += dd1[:, None] * M[:-1]
        out[:-1] += dd1[:, None] * M[1:]
        return out

    LT = lu.T
    R = lu @ tri_times(d0, d1, LT) - K
    n = K.shape[0]
    # rounding in fl(L D L^T): |E| <= gamma_n |L||D||L^T| <= gamma_n B B^T entrywise with
    # B = |L| diag(dmaj)^(1/2), dmaj a diagonal majorant of |D| (2x2 blocks: |a|+|b|, |b|+|c|);
    # ||E||_2 <= gamma_n ||B||_2^2 <= gamma_n ||B||_1 ||B||_inf
    dmaj = np.abs(d0).copy(); dmaj[1:] += np.abs(d1); dmaj[:-1] += np.abs(d1)
    B = np.abs(lu) * np.sqrt(dmaj)[None, :]
    r = float(np.linalg.norm(R) + _gamma(n) * B.sum(axis=0).max() * B.sum(axis=1).max())
    return neg, zero, r, "dense-bunch-kaufman"


def sparse_count(A, s, chunk=2000):
    """(#neg, #zero, residual bound r) for A - sI via SuperLU without row pivoting.
    Returns (None, None, inf, msg) if SuperLU permuted rows and columns differently."""
    n = A.shape[0]
    K = (sp.csc_matrix(A) - s * sp.identity(n, format="csc")).tocsc()
    lu = spla.splu(K, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                   options=dict(SymmetricMode=True))
    if not np.array_equal(lu.perm_r, lu.perm_c):
        return None, None, np.inf, "sparse: asymmetric permutation"
    D = lu.U.diagonal()
    neg = int(np.sum(D < 0)); zero = int(np.sum(D == 0))
    L = lu.L.tocsc()
    LT = L.T.tocsc()
    # K permuted: Kp[i, j] = K[perm_r^-1 ...]; SuperLU: Pr K Pc = L U with Pr[perm_r[i], i] = 1
    Pr = sp.csc_matrix((np.ones(n), (lu.perm_r, np.arange(n))), shape=(n, n))
    Pc = sp.csc_matrix((np.ones(n), (np.arange(n), lu.perm_c)), shape=(n, n))
    Kp = (Pr @ K @ Pc).tocsc()
    Dg = sp.diags(D)
    r2 = 0.0
    for j0 in range(0, n, chunk):
        j1 = min(n, j0 + chunk)
        Rc = L @ (Dg @ LT[:, j0:j1]) - Kp[:, j0:j1]
        r2 += float(np.sum(np.asarray(Rc.data) ** 2))
    # rounding in fl(L D L^T), inner products of length <= m = max row count of L:
    # ||E||_2 <= gamma_m ||B||_1 ||B||_inf, B = |L| |D|^(1/2)
    B = abs(L) @ sp.diags(np.sqrt(np.abs(D)))
    m = int(np.diff(L.tocsr().indptr).max())
    r = float(np.sqrt(r2) + _gamma(m) * abs(B).sum(axis=0).max() * abs(B).sum(axis=1).max())
    return neg, zero, r, f"sparse-superlu-nopivot(nnzL={L.nnz})"


def count_at(A, s, dense_max=DENSE_MAX):
    if A.shape[0] <= dense_max:
        return dense_count(A, s)
    return sparse_count(A, s)


def inertia_below(A, s, delta=None, dense_max=DENSE_MAX, tries=5, need_lo=True):
    """Bracket N(s) = #{eigenvalues of symmetric A below s}.

    Factors A - (s + delta) I and A - (s - delta) I and accepts a side only if its residual bound
    r <= delta. Then n_lo <= N(s) <= n_hi (see module docstring). delta starts at
    1e-14 * max(1, |s|) (or the given value) and is reset to 2 r whenever a residual bound r
    exceeds it (up to `tries` rounds), so delta is about twice the actual backward error.
    need_lo=False skips the s - delta side (only the upper count is needed to verify a lower
    eigenvalue bound). Returns dict(n_lo, n_hi, delta, r_lo, r_hi, method, exact)."""
    if delta is None:  # guess close to the typical residual bound, so usually one factorization
        fro = float(sp.linalg.norm(A)) if sp.issparse(A) else float(np.linalg.norm(A))
        delta = max(1e-14 * max(1.0, abs(s)), 4 * A.shape[0] * EPS * (fro + abs(s) * np.sqrt(A.shape[0])))
    hi = lo = None
    for _ in range(tries):
        if hi is None or hi[2] > delta:
            hi = count_at(A, s + delta, dense_max)
            if hi[2] > delta:
                delta = 2 * hi[2]
                hi = None
                continue
        if need_lo:
            lo = count_at(A, s - delta, dense_max)
            if lo[2] > delta:
                delta = 2 * lo[2]
                hi = None
                continue
        break
    ok_hi = hi is not None and hi[2] <= delta
    ok_lo = lo is not None and lo[2] <= delta
    out = dict(n_lo=lo[0] if ok_lo else None, n_hi=hi[0] if ok_hi else None, delta=delta,
               r_lo=lo[2] if lo else np.nan, r_hi=hi[2] if hi else np.nan,
               method=(hi or lo or (0, 0, 0, "-"))[3])
    out["exact"] = out["n_lo"] is not None and out["n_lo"] == out["n_hi"]
    return out


def at_most(A, s, k, weaken=True, **kw):
    """Verdict for the statement '#eig(A < s) <= k' (equivalently lambda_k(A) >= s, 0-based):
    'verified' if n_hi <= k, 'failed' if n_lo > k, else 'inconclusive'. If inconclusive and
    weaken=True, retries at s - 3 delta and reports 'verified at s - w' with w stated."""
    res = inertia_below(A, s, need_lo=False, **kw)  # the upper count alone can verify
    if not (res["n_hi"] is not None and res["n_hi"] <= k) and "delta" not in kw             and np.isfinite(res["r_hi"]) and res["r_hi"] < 0.5 * res["delta"]:
        # the starting delta was a pessimistic guess: retry with delta = 2 x measured residual
        res = inertia_below(A, s, need_lo=False, delta=2 * res["r_hi"], **kw)
    if not (res["n_hi"] is not None and res["n_hi"] <= k):
        res = inertia_below(A, s, delta=res["delta"], **kw)
    if res["n_hi"] is not None and res["n_hi"] <= k:
        v = "verified"
    elif res["n_lo"] is not None and res["n_lo"] > k:
        v = "failed"
    else:
        v = "inconclusive"
        if weaken:
            w = 3 * res["delta"]
            r2 = inertia_below(A, s - w, need_lo=False, delta=res["delta"], **kw)
            if r2["n_hi"] is not None and r2["n_hi"] <= k:
                v = f"verified at claimed - {w:.1e}"
    res["verdict"] = v
    return res
