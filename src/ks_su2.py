"""SU(2) Kogut-Susskind Hamiltonian in the gauge-invariant spin-network basis.

Lattice units (a = 1):
    H = (g^2/2) sum_links E_l^2  -  (1/g^2) sum_plaq Tr(U_p + U_p^dag)
      = (g^2/2) sum_links j_l(j_l+1)  -  (2/g^2) sum_plaq Tr U_p        (SU(2): Tr U real)
Tr is the trace in the fundamental (j=1/2) representation.

Basis: orthonormal spin networks. Link states |j m n> <-> sqrt(2j+1) D^j_{mn}(U) (Condon-Shortley),
an orthonormal basis of L^2(SU(2), Haar). At every vertex the link indices are contracted with an
orthonormal basis of invariant tensors (intertwiners), computed numerically as the common null
space of the total su(2) generators acting on the leg indices. Gauss's law holds exactly.

Plaquette matrix elements: Tr U_p multiplies each loop link by D^{1/2}; the Clebsch-Gordan series
D^{1/2}_{ab} D^j_{mn} = sum_J C^{J,a+m}_{1/2 a; j m} C^{J,b+n}_{1/2 b; j n} D^J_{a+m,b+n}
splits into one CG factor per link end; the fundamental index contracted at each plaquette corner
is local to that vertex, so <G'|Tr U_p|G> = prod_links sqrt((2j+1)/(2J+1)) * prod_corners M_v,
where M_v is the overlap of the new intertwiner with the old one dressed by the two CG maps.

Spins are stored as integers tj = 2j throughout.
"""
from __future__ import annotations

import itertools
from collections import deque
from functools import lru_cache

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp


# ---------------------------------------------------------------- su(2) basics
@lru_cache(maxsize=None)
def spin_ops(tj: int):
    """Condon-Shortley Jx, Jy, Jz for spin tj/2, basis ordered m = j, j-1, ..., -j."""
    j = tj / 2
    m = j - np.arange(tj + 1)
    jz = np.diag(m).astype(complex)
    jp = np.zeros((tj + 1, tj + 1), complex)
    for k in range(1, tj + 1):  # <m+1|J+|m>, m = m[k]
        jp[k - 1, k] = np.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    jx = (jp + jp.conj().T) / 2
    jy = (jp - jp.conj().T) / (2j)
    return jx, jy, jz


def m_index(tj: int, tm: int) -> int:
    """Row index of 2m = tm in the spin-tj/2 basis (ordered from m = j down)."""
    return (tj - tm) // 2


@lru_cache(maxsize=None)
def cg_half(tj: int, tJ: int, ts: int) -> np.ndarray:
    """Matrix (2J+1)x(2j+1) with entries C^{J, m+s}_{1/2 s; j m}, s = ts/2, J = j +- 1/2.

    Standard Condon-Shortley values for coupling 1/2 (first) with j (second).
    """
    out = np.zeros((tJ + 1, tj + 1))
    j = tj / 2
    s = ts / 2
    for tm in range(-tj, tj + 1, 2):
        tM = tm + ts
        if abs(tM) > tJ:
            continue
        M = tM / 2
        if tJ == tj + 1:
            # <1/2 s; j m | j+1/2 M>
            val = np.sqrt((j + 0.5 + 2 * s * M) / (2 * j + 1))
        elif tJ == tj - 1:
            # <1/2 s; j m | j-1/2 M> = -<j m; 1/2 s | j-1/2 M>  (order-swap phase)
            val = 2 * s * np.sqrt((j + 0.5 - 2 * s * M) / (2 * j + 1))
        else:
            raise ValueError
        out[m_index(tJ, tM), m_index(tj, tm)] = val
    return out


def sigma(ts: int) -> int:
    """sigma(+1/2)=+1, sigma(-1/2)=-1; (-1)^{b-a} = sigma(b) sigma(a) for a, b in {+-1/2}."""
    return 1 if ts > 0 else -1


# ---------------------------------------------------------------- intertwiners
@lru_cache(maxsize=None)
def intertwiners(legs: tuple) -> np.ndarray:
    """Orthonormal invariant tensors for legs = ((tj, kind), ...), kind 0 = row (link starts here),
    1 = column (link ends here). Returns array (n_int, prod dims).

    Coefficient tensors c of sum c_{..m..n..} prod D_{mn}(U) transform at vertex v as
    row legs: c -> D(g)^T c, column legs: c -> D(g^-1) c. Generators: J^T (row) and -J (column).
    """
    dims = [tj + 1 for tj, _ in legs]
    D = int(np.prod(dims)) if dims else 1
    if D == 1:
        return np.ones((1, 1), complex)
    if len(legs) <= 3:  # quick triangle / parity reject
        tjs = sorted(tj for tj, _ in legs)
        if sum(tjs) % 2 or (len(tjs) == 2 and tjs[0] != tjs[1]) or \
                (len(tjs) == 3 and tjs[2] > tjs[0] + tjs[1]) or (len(tjs) == 1 and tjs[0] > 0):
            return np.zeros((0, D), complex)
    # Invariant vectors = weight-0 vectors (G_z = 0) annihilated by G_+ = G_x + i G_y.
    idx = np.indices(dims).reshape(len(dims), -1)  # multi-index of every flat position
    ops = []
    Gz = np.zeros(D)
    for k, (tj, kind) in enumerate(legs):
        jx, jy, jz = spin_ops(tj)
        if kind == 0:
            opz, opp = jz.T, jx.T + 1j * jy.T
        else:
            opz, opp = -jz, -(jx + 1j * jy)
        Gz += np.real(np.diag(opz))[idx[k]]
        ops.append(opp)
    zero = np.where(np.abs(Gz) < 1e-9)[0]
    if len(zero) == 0:
        return np.zeros((0, D), complex)
    strides = np.array([int(np.prod(dims[k + 1:])) for k in range(len(dims))])
    rows, cols, vals = [], [], []
    for k, opp in enumerate(ops):
        r_of_c = {c: (r, opp[r, c]) for r, c in zip(*np.nonzero(np.abs(opp) > 1e-14))}
        for ci, col in enumerate(zero):
            c = idx[k, col]
            if c in r_of_c:
                r, val = r_of_c[c]
                rows.append(col + (r - c) * strides[k]); cols.append(ci); vals.append(val)
    if rows:
        urows, inv = np.unique(rows, return_inverse=True)
        A = np.zeros((len(urows), len(zero)), complex)
        np.add.at(A, (inv, np.array(cols)), np.array(vals))
    else:
        A = np.zeros((1, len(zero)), complex)
    ns = sla.null_space(A, rcond=1e-10)
    out = np.zeros((ns.shape[1], D), complex)
    out[:, zero] = ns.T
    return out  # rows: orthonormal intertwiners


# ---------------------------------------------------------------- lattices
class Lattice:
    """vertices: int count; links: list of (start, end); plaquettes: list of [(link, +1/-1), ...]
    in loop order (+1 traversed start->end)."""

    def __init__(self, name, n_vertices, links, plaquettes):
        self.name = name
        self.nv = n_vertices
        self.links = links
        self.plaqs = plaquettes
        # legs at each vertex: list of (link, kind)
        self.legs = [[] for _ in range(n_vertices)]
        for l, (s, e) in enumerate(links):
            self.legs[s].append((l, 0))
            self.legs[e].append((l, 1))
        for p in plaquettes:  # sanity: consecutive links share a vertex
            for k in range(len(p)):
                l1, d1 = p[k]
                l2, d2 = p[(k + 1) % len(p)]
                end1 = links[l1][1] if d1 > 0 else links[l1][0]
                start2 = links[l2][0] if d2 > 0 else links[l2][1]
                assert end1 == start2, (p, k)

    def vertex_legs(self, v, tjs):
        return tuple((tjs[l], kind) for l, kind in self.legs[v])


def single_plaquette():
    # square 0-1-2-3, links 0:(0,1) 1:(1,2) 2:(3,2) 3:(0,3); loop U0 U1 U2^dag U3^dag
    return Lattice("plaq1", 4, [(0, 1), (1, 2), (3, 2), (0, 3)],
                   [[(0, 1), (1, 1), (2, -1), (3, -1)]])


def ladder2():
    # vertices bottom 0,1,2 ; top 3,4,5. Two plaquettes side by side, open boundary.
    links = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
    plaqs = [[(0, 1), (5, 1), (2, -1), (4, -1)],
             [(1, 1), (6, 1), (3, -1), (5, -1)]]
    return Lattice("ladder2", 6, links, plaqs)


def periodic(Lx, Ly):
    def vid(x, y):
        return (x % Lx) + Lx * (y % Ly)
    links, xl, yl = [], {}, {}
    for y in range(Ly):
        for x in range(Lx):
            xl[(x, y)] = len(links); links.append((vid(x, y), vid(x + 1, y)))
            yl[(x, y)] = len(links); links.append((vid(x, y), vid(x, y + 1)))
    plaqs = []
    for y in range(Ly):
        for x in range(Lx):
            plaqs.append([(xl[(x, y)], 1), (yl[((x + 1) % Lx, y)], 1),
                          (xl[(x, (y + 1) % Ly)], -1), (yl[(x, y)], -1)])
    return Lattice(f"periodic{Lx}x{Ly}", Lx * Ly, links, plaqs)


# ---------------------------------------------------------------- plaquette action
def _end_map(tj, tJ, direction, end, t):
    """CG map on the leg of a loop link at one end, for trace index t (=2*a).

    direction +1: incoming index a sits at the row end, outgoing b at the column end.
    direction -1 (U^dag_{ab} = (-1)^{b-a} D_{-b,-a}): incoming a at column end, outgoing b at row end.
    end: 'in' or 'out' (which loop index this end carries).
    """
    if direction > 0:
        return cg_half(tj, tJ, t)
    return sigma(t) * cg_half(tj, tJ, -t)


@lru_cache(maxsize=None)
def corner_matrix(old_legs, new_legs, posA, dirA, posB, dirB):
    """M[i', i] = <I'_{i'} | T | I_i> at a corner. posA: leg of the loop link arriving (its
    outgoing index sits here), posB: leg of the loop link leaving (its incoming index)."""
    Iold = intertwiners(old_legs)
    Inew = intertwiners(new_legs)
    if Iold.shape[0] == 0 or Inew.shape[0] == 0:
        return None
    dims_old = [tj + 1 for tj, _ in old_legs]
    n = len(old_legs)
    tjA, tJA = old_legs[posA][0], new_legs[posA][0]
    tjB, tJB = old_legs[posB][0], new_legs[posB][0]
    res = 0
    for t in (1, -1):
        RA = _end_map(tjA, tJA, dirA, 'out', t)
        RB = _end_map(tjB, tJB, dirB, 'in', t)
        X = Iold.reshape([-1] + dims_old)
        X = np.moveaxis(np.tensordot(RA, X, axes=([1], [1 + posA])), 0, 1 + posA)
        X = np.moveaxis(np.tensordot(RB, X, axes=([1], [1 + posB])), 0, 1 + posB)
        res = res + X.reshape(X.shape[0], -1)
    M = Inew.conj() @ res.T
    if np.max(np.abs(M)) < 1e-13:
        return None
    return M


def plaquette_moves(lat: Lattice, state, p, tjmax):
    """Yield (new_state, amplitude) for Tr U_p acting on basis state = (tjs, ints)."""
    tjs, ints = state
    loop = lat.plaqs[p]
    nl = len(loop)
    # corner k sits between loop[k-1] (arriving) and loop[k] (leaving)
    corners = []
    for k in range(nl):
        lB, dB = loop[k]
        v = lat.links[lB][0] if dB > 0 else lat.links[lB][1]
        lA, dA = loop[k - 1]
        legs = lat.legs[v]
        posA = next(i for i, (l, kind) in enumerate(legs)
                    if l == lA and kind == (1 if dA > 0 else 0))
        posB = next(i for i, (l, kind) in enumerate(legs)
                    if l == lB and kind == (0 if dB > 0 else 1))
        corners.append((v, posA, dA, posB, dB))
    for shifts in itertools.product((1, -1), repeat=nl):
        new = list(tjs)
        ok = True
        amp = 1.0
        for (l, _), s in zip(loop, shifts):
            new[l] = tjs[l] + s
            if new[l] < 0 or new[l] > tjmax:
                ok = False
                break
            amp *= np.sqrt((tjs[l] + 1) / (new[l] + 1))
        if not ok:
            continue
        new = tuple(new)
        mats = []
        for (v, posA, dA, posB, dB) in corners:
            M = corner_matrix(lat.vertex_legs(v, tjs), lat.vertex_legs(v, new), posA, dA, posB, dB)
            if M is None:
                ok = False
                break
            mats.append((v, M[:, ints[v]]))
        if not ok:
            continue
        cvs = [v for v, _ in mats]
        for choice in itertools.product(*[range(len(col)) for _, col in mats]):
            a = amp
            for (v, col), c in zip(mats, choice):
                a *= col[c]
            if abs(a) < 1e-13:
                continue
            nints = list(ints)
            for v, c in zip(cvs, choice):
                nints[v] = c
            yield (new, tuple(nints)), a


def build(lat: Lattice, tjmax: int, max_dim=2_000_000):
    """BFS from the strong-coupling vacuum (all j=0). Returns basis list, Casimir sums (sum_l j(j+1)),
    and sparse matrix W = sum_p Tr U_p restricted to the connected sector."""
    vac = (tuple([0] * len(lat.links)), tuple([0] * lat.nv))
    index = {vac: 0}
    basis = [vac]
    rows, cols, vals = [], [], []
    q = deque([vac])
    while q:
        st = q.popleft()
        i = index[st]
        for p in range(len(lat.plaqs)):
            for nst, a in plaquette_moves(lat, st, p, tjmax):
                if nst not in index:
                    index[nst] = len(basis)
                    basis.append(nst)
                    q.append(nst)
                    if len(basis) > max_dim:
                        raise MemoryError("basis too large")
                rows.append(index[nst]); cols.append(i); vals.append(a)
    n = len(basis)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    cas = np.array([sum(t * (t + 2) / 4 for t in st[0]) for st in basis])
    return basis, cas, W


def hamiltonian(cas, W, g2):
    return sp.diags(0.5 * g2 * cas) - (2.0 / g2) * W


def lowest(H, k=4):
    n = H.shape[0]
    if n <= 1500:
        w = np.linalg.eigvalsh(H.toarray())
        return np.sort(w.real)[:k]
    from scipy.sparse.linalg import eigsh
    w = eigsh(H, k=k, which='SA', tol=1e-12, return_eigenvectors=False)
    return np.sort(w.real)


# ---------------------------------------------------------------- pointwise validation
def D_matrix(tj, v):
    """D^j(U) for U = exp(-i v.sigma/2)."""
    jx, jy, jz = spin_ops(tj)
    return sla.expm(-1j * (v[0] * jx + v[1] * jy + v[2] * jz))


def wavefunction(lat, state, Us):
    """Evaluate the spin network psi_G at link group elements given as rotation vectors Us[l]."""
    tjs, ints = state
    letters = iter("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    lab = {}
    for l in range(len(lat.links)):
        lab[(l, 0)] = next(letters); lab[(l, 1)] = next(letters)
    ops, specs = [], []
    for l, tj in enumerate(tjs):
        ops.append(np.sqrt(tj + 1) * D_matrix(tj, Us[l]))
        specs.append(lab[(l, 0)] + lab[(l, 1)])
    for v in range(lat.nv):
        legs = lat.vertex_legs(v, tjs)
        I = intertwiners(legs)[ints[v]]
        if len(legs) == 0:
            continue
        ops.append(I.reshape([tj + 1 for tj, _ in legs]))
        specs.append("".join(lab[(l, kind)] for l, kind in lat.legs[v]))
    return np.einsum(",".join(specs) + "->", *ops, optimize=True)
