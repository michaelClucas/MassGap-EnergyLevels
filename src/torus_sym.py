"""Symmetry-blocked SU(2) Kogut-Susskind Hamiltonian on the periodic 2x2 torus (2+1D), all Z2
electric-flux sectors. Extends `src/ks_su2.py` and `src/certify.py` without changing them.

Setting: 2+1D Hamiltonian lattice, SU(2), a = 1, periodic 2x2, finite volume, fixed g^2,
H = (g^2/2) sum_l j_l(j_l+1) - (2/g^2) sum_p Tr U_p (same as `src/ks_su2.py`).

Key device: `torus2x2()` reorders the legs at every vertex by direction (E-out, N-out, W-in, S-in).
`ks.intertwiners` is a deterministic function of the leg tuple, so every vertex then uses the same
intertwiner basis for the same four spins, and lattice translations act on spin networks as pure
permutations of (link spins, intertwiner labels), with no phases. This is a basis change at each
vertex relative to `ks.periodic(2, 2)`, so spectra are identical (checked in the work script).

Point-group elements (reflections) reverse links in general; `apply_op` implements a general
lattice automorphism, using D^j_{mn}(U^-1) = (-1)^{n-m} D^j_{-n,-m}(U) on reversed links and the
overlap of the transformed vertex tensor with the image vertex's intertwiner basis.
"""
from __future__ import annotations

import itertools

import numpy as np
import scipy.sparse as sp

from src import ks_su2 as ks
from src import certify as C

L = 2


def vid(x, y):
    return (x % L) + L * (y % L)


def xl(x, y):
    return 2 * vid(x, y)


def yl(x, y):
    return 2 * vid(x, y) + 1


def torus2x2():
    lat = ks.periodic(L, L)
    for y in range(L):
        for x in range(L):
            lat.legs[vid(x, y)] = [(xl(x, y), 0), (yl(x, y), 0), (xl(x - 1, y), 1), (yl(x, y - 1), 1)]
    return lat


# ---------------------------------------------------------------- lattice automorphisms
def make_op(fpt):
    """Lattice automorphism from a point map fpt(x, y) -> (x', y') that is an affine symmetry of the
    torus. Returns (vperm, lperm, lflip): vertex v -> vperm[v], link l -> lperm[l], reversed iff
    lflip[l]."""
    vperm = [None] * (L * L)
    lperm = [None] * (2 * L * L)
    lflip = [False] * (2 * L * L)
    ends = {}
    for y in range(L):
        for x in range(L):
            ends[(vid(x, y), vid(x + 1, y), 0)] = xl(x, y)
            ends[(vid(x, y), vid(x, y + 1), 1)] = yl(x, y)
    for y in range(L):
        for x in range(L):
            vperm[vid(x, y)] = vid(*fpt(x, y))
            for (l, nb) in ((xl(x, y), (x + 1, y)), (yl(x, y), (x, y + 1))):
                a, b = fpt(x, y), fpt(*nb)
                # image displacement (mod L, as +-1 step); 2x2: +1 == -1 mod 2, so use the
                # unwrapped affine image to decide orientation
                dx, dy = b[0] - a[0], b[1] - a[1]
                if (dx, dy) == (1, 0):
                    lperm[l] = xl(*a)
                elif (dx, dy) == (-1, 0):
                    lperm[l] = xl(b[0], b[1]); lflip[l] = True
                elif (dx, dy) == (0, 1):
                    lperm[l] = yl(*a)
                elif (dx, dy) == (0, -1):
                    lperm[l] = yl(b[0], b[1]); lflip[l] = True
                else:
                    raise ValueError((x, y, dx, dy))
    return vperm, lperm, lflip


def translation(a, b):
    return make_op(lambda x, y: (x + a, y + b))


def sigma_x():  # x -> -x
    return make_op(lambda x, y: (-x, y))


def sigma_d():  # x <-> y
    return make_op(lambda x, y: (y, x))


def rot90():  # (x, y) -> (-y, x)
    return make_op(lambda x, y: (-y, x))


TRANS = [(0, 0), (1, 0), (0, 1), (1, 1)]
_TOPS = {g: translation(*g) for g in TRANS}


def translate(state, g):
    """Pure permutation action of translation g on (tjs, ints) for `torus2x2()`."""
    vperm, lperm, _ = _TOPS[g]
    tjs, ints = state
    nt = [0] * len(tjs)
    ni = [0] * len(ints)
    for l, t in enumerate(tjs):
        nt[lperm[l]] = t
    for v, i in enumerate(ints):
        ni[vperm[v]] = i
    return tuple(nt), tuple(ni)


def canon(state):
    """(rep, h, stab): rep = min over translations, state = h . rep (Z2xZ2: h = h^-1), |stabilizer|."""
    imgs = [(translate(state, g), g) for g in TRANS]
    rep = min(s for s, _ in imgs)
    hs = [g for s, g in imgs if s == rep]
    return rep, hs[0], len(hs)


def _flip_leg(t, axis, tj, kind):
    """Reverse one leg of a vertex tensor: m -> -m with phase (-1)^(j-m) (row leg) or (-1)^(j+m)
    (column leg); index i = j - m in the m-descending basis."""
    i = np.arange(tj + 1)
    ph = (-1.0) ** i if kind == 0 else (-1.0) ** (tj - i)
    t = np.moveaxis(t, axis, 0)
    t = (ph[:, None] * t.reshape(tj + 1, -1))[::-1].reshape(t.shape)
    return np.moveaxis(t, 0, axis)


def apply_op(lat, op, state):
    """General automorphism acting on a spin network: returns list of (new_state, amplitude)."""
    vperm, lperm, lflip = op
    tjs, ints = state
    nt = [0] * len(tjs)
    for l, t in enumerate(tjs):
        nt[lperm[l]] = t
    nt = tuple(nt)
    sign = 1.0
    for l, t in enumerate(tjs):
        if lflip[l] and t % 2:
            sign = -sign
    cols = [None] * lat.nv
    for v in range(lat.nv):
        old_legs = lat.vertex_legs(v, tjs)
        I = ks.intertwiners(old_legs)[ints[v]].reshape([tj + 1 for tj, _ in old_legs])
        w = vperm[v]
        new_order = []
        for k, (l, kind) in enumerate(lat.legs[v]):
            if lflip[l]:
                I = _flip_leg(I, k, tjs[l], kind)
            nk = 1 - kind if lflip[l] else kind
            new_order.append(lat.legs[w].index((lperm[l], nk)))
        # axis k of I goes to position new_order[k]
        I = np.transpose(I, np.argsort(new_order))
        Inew = ks.intertwiners(lat.vertex_legs(w, nt))
        cols[w] = Inew.conj() @ I.reshape(-1)
    out = []
    for choice in itertools.product(*[range(len(c)) for c in cols]):
        a = sign
        for v, c in enumerate(choice):
            a *= cols[v][c]
        if abs(a) > 1e-13:
            out.append(((nt, tuple(choice)), a))
    return out


def op_matrix(lat, op, basis):
    index = {s: i for i, s in enumerate(basis)}
    rows, cols, vals = [], [], []
    for i, s in enumerate(basis):
        for ns, a in apply_op(lat, op, s):
            rows.append(index[ns]); cols.append(i); vals.append(a)
    M = sp.csr_matrix((vals, (rows, cols)), shape=(len(basis),) * 2)
    if np.iscomplexobj(M.data):
        assert np.abs(M.data.imag).max(initial=0) < 1e-12
        M = sp.csr_matrix(M.real)
    return M


# ---------------------------------------------------------------- momentum blocks
KS = [(0, 0), (1, 0), (0, 1), (1, 1)]  # momentum (kx, ky) in units of pi


def chi(k, g):
    return -1.0 if (k[0] * g[0] + k[1] * g[1]) % 2 else 1.0


_LAT = None


def _moves_of_rep(args):
    rep, tjmax = args
    global _LAT
    if _LAT is None:
        _LAT = torus2x2()
    acc = {}
    for p in range(len(_LAT.plaqs)):
        for ns, a in ks.plaquette_moves(_LAT, rep, p, tjmax):
            r2, h, _ = canon(ns)
            key = (r2, h)
            acc[key] = acc.get(key, 0.0) + a
    return [(r2, h, float(np.real(a))) for (r2, h), a in acc.items() if abs(a) > 1e-14]


def build_blocks(flux, tjmax, nproc=1, log=None):
    """Momentum-blocked H pieces for the given Z2 flux sector, complete basis up to spin tjmax/2.
    Returns dict k -> (reps_k, cas_k, W_k) with W_k real symmetric (scipy CSR)."""
    lat = torus2x2()
    states = C.sector_states(lat, tjmax, flux)
    reps, stab = {}, {}
    for s in states:
        r, _, ns = canon(s)
        if r not in reps:
            reps[r] = len(reps); stab[r] = ns
    rep_list = list(reps)
    jobs = [(r, tjmax) for r in rep_list]
    if nproc > 1:
        import multiprocessing as mp
        res = []
        with mp.get_context("spawn").Pool(nproc) as pool:
            for i, r in enumerate(pool.imap(_moves_of_rep, jobs, chunksize=32)):
                res.append(r)
                if log and (i + 1) % 2000 == 0:
                    log(f"  {i + 1}/{len(jobs)} reps")
    else:
        res = [_moves_of_rep(j) for j in jobs]
    out = {}
    nrep = len(rep_list)
    S = np.array([stab[r] for r in rep_list], float)
    tr = [(i, reps[r2], h, a) for i, mv in enumerate(res) for (r2, h, a) in mv]
    ci = np.array([t[0] for t in tr]); ri = np.array([t[1] for t in tr])
    hx = np.array([t[2][0] for t in tr]); hy = np.array([t[2][1] for t in tr])
    av = np.array([t[3] for t in tr])
    del tr, res
    # stabilizer subgroup of each rep (for the momentum selection rule)
    stabg = [[g for g in TRANS if translate(r, g) == r] for r in rep_list]
    for k in KS:
        ok = np.array([all(chi(k, g) == 1.0 for g in sg) for sg in stabg])
        idx = np.where(ok)[0]
        pos = -np.ones(nrep, int); pos[idx] = np.arange(len(idx))
        m = ok[ci] & ok[ri]
        ph = np.where((k[0] * hx[m] + k[1] * hy[m]) % 2 == 1, -1.0, 1.0)
        vals = av[m] * ph * np.sqrt(S[ri[m]] / S[ci[m]])
        W = sp.csr_matrix((vals, (pos[ri[m]], pos[ci[m]])), shape=(len(idx), len(idx)))
        W.sum_duplicates(); W.eliminate_zeros()
        cas = np.array([sum(t * (t + 2) / 4 for t in rep_list[i][0]) for i in idx])
        out[k] = ([rep_list[i] for i in idx], cas, W, S[idx])
    return out, len(states)


def lift(reps_k, stab_k, k, coeffs, basis_index, n):
    """Momentum-basis vector -> full spin-network basis vector."""
    v = np.zeros(n)
    G = len(TRANS)
    for r, s, c in zip(reps_k, stab_k, coeffs):
        nrm = 1.0 / np.sqrt(G * s)
        for g in TRANS:
            v[basis_index[translate(r, g)]] += nrm * chi(k, g) * c
    return v


# ---------------------------------------------------------------- improved electric floor
def plaquette_e0_lower(kin, g2, J2max=400):
    """Certified lower bound on inf spec of h = kin * j(j+1) - (2/g2) chi_{1/2}(U) on class functions
    of one plaquette loop (character basis, tridiagonal: diagonal kin*j(j+1), off-diagonal -2/g2). Truncated at 2j <= J2max; the tail is
    handled by the Schur Lemma of src/certify.py with floor kin*J(J+1) - 4/g2 on the tail and the
    single coupling entry, so the result is a lower bound for the untruncated operator."""
    tj = np.arange(J2max + 1)
    diag = kin * tj * (tj + 2) / 4
    A = sp.diags([diag, -2.0 / g2 * np.ones(J2max), -2.0 / g2 * np.ones(J2max)], [0, 1, -1]).tocsr()
    Jt = (J2max + 1) / 2
    d = kin * Jt * (Jt + 1) - 4.0 / g2
    X = sp.csr_matrix(([(2.0 / g2) ** 2], ([J2max], [J2max])), shape=A.shape)
    U = C.upper_bounds(A, 1)[0] + C._eps_margin(A)
    if not U < d:
        return -np.inf
    Ad = (A - X / (d - U)).toarray()
    import scipy.linalg as sla
    return sla.eigvalsh_tridiagonal(np.diag(Ad).copy(), np.diag(Ad, 1).copy(),
                                    select="i", select_range=(0, 0))[0] - C._eps_margin(A)


def split_floor(g2, c_Q, n_plaq=4, ts=np.linspace(0, 1, 41)):
    """Proved floor on Q H Q (Q: gauge-invariant states with sum_l j(j+1) >= c_Q on the 2x2 torus):
    H = (1-t)(g^2/2) C + sum_p h_p(t), h_p(t) = t (g^2/4) sum_{l in p} C_l - (2/g^2) Tr U_p
    (each link lies in exactly 2 plaquettes). inf spec h_p(t) on L^2(SU(2)^4) is attained in the
    corner-gauge-invariant sector (unique positive ground state), i.e. class functions of U_p, where
    sum_{l in p} C_l acts as 4 j(j+1). Returns (best floor, t)."""
    best, tb = -np.inf, None
    for t in ts:
        e = -4.0 / g2 if t == 0 else plaquette_e0_lower(t * g2, g2)  # 4 * t g^2/4 = t g^2
        f = (1 - t) * 0.5 * g2 * c_Q + n_plaq * e
        if f > best:
            best, tb = f, t
    return best, tb
