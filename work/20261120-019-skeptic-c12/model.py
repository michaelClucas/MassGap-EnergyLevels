"""Constituent-gluon bound-state model for glueballs (Result-B attempt, work 013).

Units: sqrt(sigma_fundamental) = 1 throughout.
Two-gluon states: H = K + sigma_A r - kappa/r (+ C0), kappa = N alpha_s, sigma_A = (C_A/C_F) sigma.
Kinetic K: spinless Salpeter 2 sqrt(p^2 + m^2) ("rel") or 2m + p^2/m ("nr").
Spin-dependent terms: Breit-Fermi form with V_V = -kappa/r, V_S = sigma_A r, first-order PT.
Three-gluon (C = -, d-coupled) states: hyperspherical K=0 approximation, pairwise Casimir potentials.
2+1D: 2D radial operator, log Coulomb with g^2 N fixed by lattice sqrt(sigma)/(g^2 N).
"""
import csv
import math
import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"  # 1 thread: multithreaded eigh on small matrices thrashes on the pinned cores
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.special import erf, sph_harm_y

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


# ---------------------------------------------------------------- data
def load(fname, N, dim, source=None):
    """Return {(JPC, level): (value, err)} for quantity m/sqrt_sigma."""
    out = {}
    with open(DATA / fname, newline="") as f:
        for row in csv.DictReader(f):
            if row["N"] != str(N) or row["dim"] != str(dim) or row["quantity"] != "m/sqrt_sigma":
                continue
            if source and row["source_id"] != source:
                continue
            key = (row["JPC"], int(row["level"]))
            if key in out:  # duplicate rows (AT2016 SU(3) 2-- typo): keep first
                continue
            out[key] = (float(row["value"]), float(row["err_stat"]))
    return out


def load_scalar(fname, N, quantity):
    with open(DATA / fname, newline="") as f:
        for row in csv.DictReader(f):
            if row["N"] == str(N) and row["quantity"] == quantity:
                return float(row["value"]), float(row["err_stat"])
    raise KeyError(quantity)


# ---------------------------------------------------------------- angular algebra
def cg(j1, m1, j2, m2, J, M):
    """Clebsch-Gordan <j1 m1 j2 m2 | J M>, integer spins, Racah formula."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J or not (abs(j1 - j2) <= J <= j1 + j2):
        return 0.0
    f = math.factorial
    pre = math.sqrt((2 * J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J) / f(j1 + j2 + J + 1))
    pre *= math.sqrt(f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2))
    s = 0.0
    for k in range(0, j1 + j2 + J + 1):
        d = [k, j1 + j2 - J - k, j1 - m1 - k, j2 + m2 - k, J - j2 + m1 + k, J - j1 - m2 + k]
        if min(d) < 0:
            continue
        s += (-1) ** k / np.prod([f(x) for x in d])
    return pre * s


def spin1_ops():
    """Cartesian spin-1 matrices in the |m=+1,0,-1> basis."""
    s2 = math.sqrt(2)
    Sp = np.array([[0, s2, 0], [0, 0, s2], [0, 0, 0]], dtype=complex)
    Sx = (Sp + Sp.conj().T) / 2
    Sy = (Sp - Sp.conj().T) / (2j)
    Sz = np.diag([1, 0, -1]).astype(complex)
    return [Sx, Sy, Sz]


@lru_cache(None)
def tensor_me(L, S, J):
    """<(L S) J, M=J | 3 (S1.n)(S2.n) - S1.S2 | same>, two spin-1 constituents, by sphere quadrature."""
    ms = [1, 0, -1]
    # coupled spin state |S mS> in 3x3 product basis
    def spin_state(mS):
        v = np.zeros(9, dtype=complex)
        for i, m1 in enumerate(ms):
            for k, m2 in enumerate(ms):
                v[3 * i + k] = cg(1, m1, 1, m2, S, mS)
        return v
    xg, wg = np.polynomial.legendre.leggauss(40)
    th = np.arccos(xg)
    ph = np.linspace(0, 2 * np.pi, 80, endpoint=False)
    Sop = spin1_ops()
    I3 = np.eye(3)
    S1 = [np.kron(s, I3) for s in Sop]
    S2 = [np.kron(I3, s) for s in Sop]
    S1S2 = sum(a @ b for a, b in zip(S1, S2))
    total = 0.0
    for t, w in zip(th, wg):
        for p in ph:
            n = [math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)]
            psi = np.zeros(9, dtype=complex)
            for mL in range(-L, L + 1):
                mS = J - mL
                if abs(mS) > S:
                    continue
                c = cg(L, mL, S, mS, J, J)
                if c == 0:
                    continue
                psi += c * sph_harm_y(L, mL, t, p) * spin_state(mS)
            A = sum(n[i] * S1[i] for i in range(3))
            B = sum(n[i] * S2[i] for i in range(3))
            op = 3 * A @ B - S1S2
            total += w * (2 * np.pi / len(ph)) * np.real(psi.conj() @ op @ psi)
    return total


def s1s2(S, ncon=2):
    """sum over pairs S_i.S_j for spin-1 constituents with total spin S."""
    return (S * (S + 1) - 2 * ncon) / 2.0


def ls(L, S, J):
    return (J * (J + 1) - L * (L + 1) - S * (S + 1)) / 2.0


# ---------------------------------------------------------------- radial solvers
class Grid:
    def __init__(self, R=16.0, n=700):
        self.h = R / (n + 1)
        self.r = self.h * np.arange(1, n + 1)
        self.D2 = (np.diag(-2 * np.ones(n)) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)) / self.h ** 2


_cache = {}


def p2_eig(grid, cent):
    """Eigen-decomposition of -d^2/dr^2 + cent/r^2 (Dirichlet), cached."""
    key = (id(grid), cent)
    if key not in _cache:
        T = -grid.D2 + np.diag(cent / grid.r ** 2)
        _cache[key] = np.linalg.eigh(T)
    return _cache[key]


def solve(grid, cent, V, m, nconst, kin, nlev=3):
    """Lowest eigenpairs of nconst * sqrt(P^2/nconst' + m^2) + V on the radial grid.

    Two-body: P^2 is the relative p^2, kinetic 2 sqrt(p^2 + m^2).
    Three-body K=0: P^2 = p_x^2 + p_y^2 (Jacobi), kinetic ~ 3 sqrt(P^2/3 + m^2) (symmetric approximation).
    """
    t, U = p2_eig(grid, cent)
    t = np.clip(t, 0, None)
    if nconst == 2:
        p2 = t
    else:
        p2 = t / 3.0
    if kin == "rel":
        k = nconst * np.sqrt(p2 + m ** 2)
    else:
        k = nconst * (m + p2 / (2 * m))
    H = (U * k) @ U.T + np.diag(V)
    e, W = np.linalg.eigh(H)
    W = W[:, :nlev] / math.sqrt(grid.h)  # normalize so sum u^2 h = 1
    return e[:nlev], W


def avg(grid, u, f):
    return float(np.sum(u ** 2 * f) * grid.h)


# ---------------------------------------------------------------- 3+1D spectrum
def casimir_ratio(N):
    return 2.0 * N ** 2 / (N ** 2 - 1)  # C_A / C_F


TWO_GLUON_4D = []  # (L, S, J, JPC) with L+S even (Bose, symmetric colour singlet), P=(-1)^L, C=+
for L in range(0, 4):
    for S in range(0, 3):
        if (L + S) % 2:
            continue
        for J in range(abs(L - S), L + S + 1):
            P = "+" if L % 2 == 0 else "-"
            TWO_GLUON_4D.append((L, S, J, f"{J}{P}+"))


def spectrum4d(params, N=3, kin="rel", grid=None, spin=True, smear="compton", Lmax=3, three=True):
    """Return list of (JPC, mass, label) for 3+1D, params = (m, kappa, C0, delta3)."""
    m, kappa, C0 = params[:3]
    grid = grid or GRID4
    sA = casimir_ratio(N)
    r = grid.r
    rs = 1.0 / m if smear == "compton" else smear
    # one-gluon-exchange potential of a Gaussian-smeared colour charge (width rs, Godfrey-Isgur style):
    # V_V = -kappa erf(r/rs)/r, so grad^2 V_V = 4 pi kappa delta_s(r) with a normalized Gaussian delta_s.
    x = r / rs
    g = np.exp(-x ** 2)
    VV = -kappa * erf(x) / r
    dVV = kappa * (erf(x) / r ** 2 - (2 / math.sqrt(math.pi)) * g / (rs * r))
    d2VV = kappa * ((4 / math.sqrt(math.pi)) * g / (rs * r ** 2) - 2 * erf(x) / r ** 3
                    + (4 / math.sqrt(math.pi)) * g / rs ** 3)
    delta_s = g / (math.pi ** 1.5 * rs ** 3)
    V = sA * r + VV
    out = []
    for L in range(0, Lmax + 1):
        e, W = solve(grid, L * (L + 1), V, m, 2, kin, nlev=3)
        for n in range(3):
            u = W[:, n]
            d0 = avg(grid, u, delta_s)  # <delta_s> = int delta_s(r) u^2 dr
            fls = avg(grid, u, (3 * dVV - sA) / r)
            ften = avg(grid, u, dVV / r - d2VV)
            for (LL, S, J, jpc) in TWO_GLUON_4D:
                if LL != L:
                    continue
                dE = 0.0
                if spin:
                    # V_SS = (2/3m^2) grad^2 V_V S1.S2
                    dE += (2 / (3 * m ** 2)) * 4 * math.pi * kappa * d0 * s1s2(S)
                    if L > 0 and S > 0:
                        # V_LS = (1/2m^2)(1/r)(3 V_V' - V_S') L.S
                        dE += (1 / (2 * m ** 2)) * fls * ls(L, S, J)
                        # V_T = (1/3m^2)(V_V'/r - V_V'') [3 (S1.n)(S2.n) - S1.S2]  (spin-operator form)
                        dE += (1 / (3 * m ** 2)) * ften * tensor_me(L, S, J)
                out.append((jpc, e[n] + dE + C0, f"gg L={L} S={S} n={n}"))
    if three:
        out += three_gluon_4d(params, N, kin, spin)
    return out


def three_gluon_4d(params, N, kin, spin):
    """d-coupled (C = -) S-wave three-gluon states: spin symmetric S = 1, 3 -> 1--, 3--."""
    m, kappa, C0 = params[:3]
    grid = GRID4_3
    sA = casimir_ratio(N)
    rho = grid.r
    c1 = 32 / (15 * math.pi)  # <cos a>, K=0, 3+3 dims
    c2 = 16 / (3 * math.pi)   # <1/cos a>
    # pairwise: sigma_A/2 * r_ij - (kappa/2)/r_ij, r_ij = sqrt2 rho cos a; 3 pairs
    # hyperangle average of the smeared Coulomb approximated by smearing the averaged one (rs = 1/m)
    V = 3 * ((sA / 2) * math.sqrt(2) * c1 * rho
             - (kappa / 2) * c2 / (math.sqrt(2) * rho) * erf(math.sqrt(2) * rho * m))
    e, W = solve(grid, 15 / 4, V, m, 3, kin, nlev=2)
    out = []
    for n in range(2):
        u = W[:, n]
        # <delta^3(r_12)> = <delta^3(x)>/2^{3/2}, <delta^3(x)> = (4/pi^2) int u^2/rho^3 (K=0, 6D)
        # <delta_s(r_12)>, r_12 = sqrt2 rho cos a, K=0 hyperangle measure (16/pi) cos^2 a sin^2 a da
        a, wa = np.polynomial.legendre.leggauss(48)
        a = (a + 1) * math.pi / 4
        wa = wa * math.pi / 4 * (16 / math.pi) * np.cos(a) ** 2 * np.sin(a) ** 2
        rs = 1.0 / m
        r12 = math.sqrt(2) * np.outer(rho, np.cos(a))
        fd = (np.exp(-(r12 / rs) ** 2) / (math.pi ** 1.5 * rs ** 3)) @ wa
        d12 = avg(grid, u, fd)
        for S, jpc in ((1, "1--"), (3, "3--")):
            dE = 0.0
            if spin:
                # pair (colour factor kappa/2), summed over pairs: sum S_i.S_j = s1s2(S, 3)
                dE = (2 / (3 * m ** 2)) * 4 * math.pi * (kappa / 2) * d12 * s1s2(S, 3)
            out.append((jpc, e[n] + dE + 1.5 * C0, f"ggg K=0 S={S} n={n}"))
    return out


# ---------------------------------------------------------------- 2+1D spectrum
def spectrum3d(params, N, g2N, kin="rel"):
    """2+1D: spinless constituents (one polarization), intrinsic parity -, L even (Bose).
    V = sigma_A r + (g^2 N / 2 pi) ln r + C0. J = |L|; J>0 parity doubled (J^{+C} = J^{-C})."""
    m, C0 = params
    grid = GRID3
    sA = casimir_ratio(N)
    r = grid.r
    V = sA * r + (g2N / (2 * math.pi)) * np.log(r)
    out = []
    for L in (0, 2, 4):
        e, _ = solve(grid, L ** 2 - 0.25, V, m, 2, kin, nlev=4)
        for n in range(4):
            jpcs = ["0++"] if L == 0 else [f"{L}++", f"{L}-+"]
            for j in jpcs:
                out.append((j, e[n] + C0, f"gg |L|={L} n={n}"))
    # d-coupled three gluons, K=0 in 4D Jacobi space: intrinsic parity (-)^3, C = - -> 0--
    rho = GRID3_3.r
    c1, cl = 2 / 3, -0.5  # <cos a>, <ln cos a> for 2+2 dims
    # pair: (sA/2) r_ij + (g2N/2)/(2 pi) ln r_ij,  r_ij = sqrt2 rho cos a, 3 pairs
    V3 = 3 * ((sA / 2) * math.sqrt(2) * c1 * rho + (g2N / 2) / (2 * math.pi) * (np.log(math.sqrt(2) * rho) + cl))
    e, _ = solve(GRID3_3, 0.75, V3, m, 3, kin, nlev=3)
    for n in range(3):
        out.append(("0--", e[n] + 1.5 * C0, f"ggg K=0 n={n}"))
    return out


def by_channel(spec):
    """Sort within J^PC: level 0, 1, ... (ordering rule: energy order within a channel)."""
    d = {}
    for jpc, E, lab in spec:
        d.setdefault(jpc, []).append((E, lab))
    return {k: sorted(v) for k, v in d.items()}


GRID4 = Grid(R=14.0, n=700)
GRID4_3 = Grid(R=14.0, n=700)
GRID3 = Grid(R=16.0, n=900)
GRID3_3 = Grid(R=16.0, n=700)
