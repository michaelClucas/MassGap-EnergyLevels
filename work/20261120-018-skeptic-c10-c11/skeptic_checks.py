"""Skeptic re-derivations for C10 (work 014) and C11 (work 015).

Independent code: nothing is imported from work 014/015. Sections:
  P  Peierls gain: closed form via complete elliptic integral E(m) vs (1/pi) d^2 ln(1/d)
  G  Gribov Schwinger function: prefactor and sign change
  D  anharmonic oscillator: exact gaps by sinc-DVR (spectrally accurate) vs Gaussian Hartree gap
  W  2D lattice Watson integral asymptotics W_2(m) ~ (1/4pi) ln(32/m^2)
  A  AKLT open-chain gaps eps_n (projector normalization), Knabe and Gosset-Mozgunov constants,
     dense (no Lanczos) certificate at n = 10, Lanczos at n = 12, ring gap n = 12
  R  rotor chain: own matrix elements; Rayleigh-Schroedinger PT for E0 and k = 0 magnon
     through x^4; spot checks of certified intervals by large-l_max ED
  B  Bakry-Emery: Riemannian Hessian of -sum n_i.n_{i+1} on (S^2)^N, minimum per site
Run: python tools/pin.py .venv/Scripts/python.exe -u work/20261120-018-skeptic-c10-c11/skeptic_checks.py
"""
import itertools
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy import integrate, special, optimize

T0 = time.time()


def hdr(s):
    print(f"\n== {s} ==  ({time.time() - T0:.0f}s)", flush=True)


# ---------------------------------------------------------------- P
hdr("P. Peierls: E/N = -(2t/pi) E(m = 1 - d^2) for half-filled spinless dimerized ring")
# band: e(k) = +-2t sqrt(cos^2 k + d^2 sin^2 k), k in (-pi/2, pi/2], lower band filled, 2 sites per cell
for d in [0.1, 0.01, 0.001]:
    num = -integrate.quad(lambda k: 2 * np.sqrt(np.cos(k) ** 2 + d * d * np.sin(k) ** 2), -np.pi / 2, np.pi / 2,
                          epsabs=1e-14, epsrel=1e-14)[0] / (2 * np.pi)  # per site
    closed = -(2 / np.pi) * special.ellipe(1 - d * d)
    gain = -(2 / np.pi) - closed
    asym = (d * d / np.pi) * (np.log(1 / d) + np.log(4) - 0.5)
    print(f"d={d:6g}  E/N quad {num:.12f}  closed {closed:.12f}  gain {gain:.6e}  "
          f"(1/pi)d^2[ln(1/d)+ln4-1/2] {asym:.6e}  ratio {gain / asym:.6f}")
print("=> gain = (1/pi) d^2 ln(1/d) + ((ln 4 - 1/2)/pi) d^2 + o(d^2); b = %.5f" % ((np.log(4) - 0.5) / np.pi))

# ---------------------------------------------------------------- G
hdr("G. Gribov D(p) = p^2/(p^4+gamma^4), C(t) = int dp/(2pi) D(p) e^{ipt}")
for gam in [1.0, 2.0]:
    for t in [0.0, 1.0, 2.0, 4.0]:
        f = lambda p: p * p / (p ** 4 + gam ** 4) / np.pi
        if t == 0:
            num = integrate.quad(f, 0, np.inf)[0]
        else:
            num = integrate.quad(f, 0, np.inf, weight="cos", wvar=t)[0]
        a = gam * t / np.sqrt(2)
        closed = np.exp(-a) * (np.cos(a) - np.sin(a)) / (2 * np.sqrt(2) * gam)
        print(f"gamma={gam} t={t}: numeric {num:+.8e}  e^-a(cos a - sin a)/(2 sqrt2 gamma) {closed:+.8e}")
print("=> prefactor is 1/(2 sqrt2 gamma) (dimension 1/mass); work-014 README writes gamma/(2 sqrt2): typo,"
      " sign change at a = pi/4 unaffected")

# ---------------------------------------------------------------- D
hdr("D. anharmonic oscillator, sinc-DVR exact vs Gaussian (Hartree) symmetric gap Omega")


def dvr_levels(mu, lam, xmax=9.0, n=900, k=2):
    x = np.linspace(-xmax, xmax, n)
    h = x[1] - x[0]
    i = np.arange(n)
    dif = i[:, None] - i[None, :]
    with np.errstate(divide="ignore"):
        T = np.where(dif == 0, np.pi ** 2 / 3, 2.0 * (-1.0) ** dif / np.where(dif == 0, 1, dif) ** 2)
    T = T / (2 * h * h)
    H = T + np.diag(0.5 * mu * x * x + lam * x ** 4)
    from scipy.linalg import eigh
    return eigh(H, eigvals_only=True, subset_by_index=[0, k - 1])


def hartree_sym(mu, lam):
    # Omega^2 = mu + 12 lam <x^2>, <x^2> = 1/(2 Omega)
    return optimize.brentq(lambda W: W * W - mu - 6 * lam / W, 1e-6, 100)


for mu, lam in [(1, .1), (1, 1), (0, 1), (-1, .1), (-2, .1), (-3, .1), (-4, .1)]:
    e1 = dvr_levels(mu, lam, n=700)
    e2 = dvr_levels(mu, lam, n=1000, xmax=10)
    gap = e2[1] - e2[0]
    W = hartree_sym(mu, lam)
    print(f"mu={mu:+.0f} lam={lam}: E0 {e2[0]:.10f}  gap {gap:.6e} (grid change {abs((e1[1]-e1[0])-gap):.1e})"
          f"  Omega_sym {W:.4f}  ratio {W / gap:.4g}")

# ---------------------------------------------------------------- W
hdr("W. 2D Watson integral W_2(m) = int d^2k/(2pi)^2 1/(m^2 + 4 - 2cos k1 - 2cos k2)")
for m in [1e-2, 1e-3, 1e-4]:
    # integrate over k1 analytically: int dk1/(2pi) 1/(a - 2cos k1) = 1/sqrt(a^2-4)
    f = lambda k2: 1 / np.sqrt((m * m + 4 - 2 * np.cos(k2)) ** 2 - 4) / np.pi
    w = integrate.quad(f, 0, np.pi, limit=400, points=[1e-3, 1e-2, 1e-1])[0]
    print(f"m={m:g}: W_2 {w:.8f}  (1/4pi) ln(32/m^2) {np.log(32 / m / m) / (4 * np.pi):.8f}")

# ---------------------------------------------------------------- A
hdr("A. AKLT: h = P_2 (projector), open-chain eps_n, Knabe / Gosset-Mozgunov")
Sz = np.diag([1.0, 0.0, -1.0])
Sp = np.sqrt(2) * np.diag([1.0, 1.0], 1)
Sm = Sp.T
SS = np.kron(Sz, Sz) + 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp))
P2 = 0.5 * SS + SS @ SS / 6 + np.eye(9) / 3
assert np.allclose(P2 @ P2, P2) and abs(np.trace(P2) - 5) < 1e-12
P2s = sp.csr_matrix(P2)


def chain(n, periodic):
    D = 3 ** n
    H = sp.csr_matrix((D, D))
    bonds = [(i, i + 1) for i in range(n - 1)] + ([(n - 1, 0)] if periodic else [])
    for (i, j) in bonds:
        if j == i + 1:
            H = H + sp.kron(sp.kron(sp.identity(3 ** i), P2s), sp.identity(3 ** (n - i - 2)), format="csr")
        else:  # wrap bond (n-1, 0): permute via explicit construction
            ops = []
            for a in range(9):
                for b in range(9):
                    if P2[a, b] != 0:
                        a0, a1 = divmod(a, 3)  # (site n-1, site 0)
                        b0, b1 = divmod(b, 3)
                        e_last = sp.csr_matrix(([1.0], ([a0], [b0])), shape=(3, 3))
                        e_first = sp.csr_matrix(([1.0], ([a1], [b1])), shape=(3, 3))
                        ops.append(P2[a, b] * sp.kron(sp.kron(e_first, sp.identity(3 ** (n - 2))), e_last))
            H = H + sum(ops)
    mz = np.zeros(1, dtype=int)
    for _ in range(n):
        mz = (mz[:, None] + np.array([1, 0, -1])[None, :]).ravel()
    idx = np.where(mz == 0)[0]
    return H.tocsr()[idx][:, idx]


def knabe(n, e):
    return (n - 1) / (n - 2) * (e - 1 / (n - 1))


def gm(n, e):
    return 5 / 6 * (n * n + n) / (n * n - 4) * (e - 6 / (n * (n + 1)))


for n in [6, 8, 10]:
    H = chain(n, False).toarray()
    ev = np.linalg.eigvalsh(H)
    nz = ev[ev > 1e-9][0]
    print(f"open n={n} dense (dim {H.shape[0]}): #zero modes in Sz=0 {np.sum(ev < 1e-9)}, eps_n {nz:.8f}, "
          f"Knabe {knabe(n, nz):.5f}, GM {gm(n, nz):.5f}")
for n in [12]:
    H = chain(n, False)
    ev = np.sort(sla.eigsh(H, k=6, which="SA", tol=1e-12)[0])
    nz = ev[ev > 1e-9][0]
    print(f"open n={n} Lanczos (dim {H.shape[0]}): lowest {np.round(ev, 8)}  eps_n {nz:.8f}, "
          f"Knabe {knabe(n, nz):.5f}, GM {gm(n, nz):.5f}")
    R = chain(n, True)
    evr = np.sort(sla.eigsh(R, k=3, which="SA", tol=1e-12)[0])
    print(f"ring n={n}: lowest {np.round(evr, 8)}  gap {evr[1] - evr[0]:.6f}")

# ---------------------------------------------------------------- R
hdr("R. O(3) rotor chain, H/g = (1/2) sum L^2 + x V, V = -sum n.n, x = 1/g^2")


def site_ops(lmax):
    st = [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]
    ix = {s: k for k, s in enumerate(st)}
    d = len(st)
    nz = np.zeros((d, d))
    npl = np.zeros((d, d))
    for (l, m), k in ix.items():
        # cos(th) Y_lm = A(l+1,m) Y_{l+1,m} + A(l,m) Y_{l-1,m}
        for (l2, c) in [(l + 1, np.sqrt(((l + 1) ** 2 - m * m) / ((2 * l + 1) * (2 * l + 3)))),
                        (l - 1, np.sqrt((l * l - m * m) / ((2 * l - 1) * (2 * l + 1))) if l > 0 else 0)]:
            if (l2, m) in ix:
                nz[ix[(l2, m)], k] += c
        # sin(th) e^{i phi} Y_lm = -sqrt((l+m+1)(l+m+2)/((2l+1)(2l+3))) Y_{l+1,m+1}
        #                          + sqrt((l-m)(l-m-1)/((2l-1)(2l+1))) Y_{l-1,m+1}
        c1 = -np.sqrt((l + m + 1) * (l + m + 2) / ((2 * l + 1) * (2 * l + 3)))
        if (l + 1, m + 1) in ix:
            npl[ix[(l + 1, m + 1)], k] += c1
        if l > 0 and (l - 1, m + 1) in ix:
            npl[ix[(l - 1, m + 1)], k] += np.sqrt((l - m) * (l - m - 1) / ((2 * l - 1) * (2 * l + 1)))
    L2 = np.array([l * (l + 1) for (l, m) in st], float)
    M = np.array([m for (l, m) in st])
    return nz, npl, L2, M, st


nz_, np_, _, _, st_ = site_ops(4)
nn1 = nz_ @ nz_ + 0.5 * (np_ @ np_.T + np_.T @ np_)
low = [k for k, (l, m) in enumerate(st_) if l <= 3]
print("site check: max |n.n - 1| on l <= lmax-1:", np.abs(nn1[np.ix_(low, low)] - np.eye(len(low))).max())


def rotor_H(L, lmax, periodic, Mtot=0):
    nz, npl, L2, Ms, st = site_ops(lmax)
    d = len(st)
    Mt = np.zeros(1, dtype=int)
    for _ in range(L):
        Mt = (Mt[:, None] + Ms[None, :]).ravel()
    idx = np.where(Mt == Mtot)[0]
    K = np.zeros(1)
    for _ in range(L):
        K = (K[:, None] + L2[None, :]).ravel()
    H0 = 0.5 * K[idx]
    bond2 = sp.csr_matrix(np.kron(nz, nz) + 0.5 * (np.kron(npl, npl.T) + np.kron(npl.T, npl)))
    bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if periodic and L > 2 else [])
    V = sp.csr_matrix((len(idx), len(idx)))
    for (i, j) in bonds:
        if j == i + 1:
            op = sp.kron(sp.kron(sp.identity(d ** i), bond2), sp.identity(d ** (L - i - 2)), format="csr")
        else:  # n_{L-1}.n_0
            Z = sp.csr_matrix(nz)
            P = sp.csr_matrix(npl)
            mid = sp.identity(d ** (L - 2))
            op = (sp.kron(sp.kron(Z, mid), Z) + 0.5 * (sp.kron(sp.kron(P, mid), P.T) +
                                                       sp.kron(sp.kron(P.T, mid), P))).tocsr()
        V = V - op[idx][:, idx]
    return H0, V.tocsr(), idx, st, d


def rspt(H0, V, psi0, E00, nmax):
    deg = np.abs(H0 - E00) < 1e-12
    Rd = np.where(deg, 0.0, 1.0 / np.where(deg, 1.0, E00 - H0))
    psis = [psi0]
    E = [E00]
    for n in range(1, nmax + 1):
        E.append(psi0 @ (V @ psis[n - 1]))
        rhs = V @ psis[n - 1] - sum(E[k] * psis[n - k] for k in range(1, n + 1))
        psis.append(Rd * rhs)
    return E


for (L, lmax, nmax) in [(5, 2, 3), (5, 3, 4)]:
    H0, V, idx, st, d = rotor_H(L, lmax, True)
    pos = {s: k for k, s in enumerate(st)}
    full_to_sec = {f: k for k, f in enumerate(idx)}
    vac = 0
    psi0 = np.zeros(len(idx)); psi0[full_to_sec[vac]] = 1
    mag = np.zeros(len(idx))
    for i in range(L):
        f = pos[(1, 0)] * d ** (L - 1 - i)
        mag[full_to_sec[f]] = 1 / np.sqrt(L)
    Eg = rspt(H0, V, psi0, 0.0, nmax)
    Em = rspt(H0, V, mag, 1.0, nmax)
    print(f"ring L={L} lmax={lmax} (M=0 dim {len(idx)}): E0/g coeffs {np.round(Eg, 10)}")
    print(f"   gap/g coeffs {np.round(np.array(Em) - np.array(Eg), 10)}")
print("exact rationals: -2/3 = %.10f, 2/27 = %.10f, 7/405 = %.10f; E0 x^2 per site -1/6" % (-2 / 3, 2 / 27, 7 / 405))

hdr("R2. spot checks of certified rotor intervals (own ED, larger l_max, upper bounds by min-max)")


def gap_ed(L, lmax, periodic, g):
    H0, V, idx, st, d = rotor_H(L, lmax, periodic)
    H = sp.diags(g * H0) + (1 / g) * V
    ev = np.sort(sla.eigsh(H.tocsc(), k=4, which="SA", tol=1e-12)[0])
    return ev, len(idx)


for (L, per, g, lm, interval) in [(2, False, 0.25, 12, (0.1318, 0.1318)),
                                  (2, False, 0.5, 12, (0.2879914092, 0.288024011)),
                                  (3, False, 0.35, 7, (0.13109, 0.13204)),
                                  (3, True, 0.35, 7, (0.12446, 0.12824))]:
    ev, dim = gap_ed(L, lm, per, g)
    print(f"L={L} {'ring' if per else 'open'} g={g} lmax={lm} dim={dim}: E0 {ev[0]:.8f}  levels {np.round(ev, 6)}"
          f"  gap {ev[1] - ev[0]:.6f}  015 interval {interval}")

# ---------------------------------------------------------------- B
hdr("B. Bakry-Emery: Riemannian Hessian of S = -sum n_i.n_{i+1} on (S^2)^N")
rng = np.random.default_rng(1)


def hess_min(ns, periodic=True):
    N = len(ns)
    # tangent basis per site
    basis = []
    for n in ns:
        a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = a - (a @ n) * n; e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        basis.append((e1, e2))
    Hm = np.zeros((2 * N, 2 * N))
    bonds = [(i, (i + 1) % N) for i in range(N if periodic else N - 1)]
    for (i, j) in bonds:
        c = ns[i] @ ns[j]
        for a in range(2):
            Hm[2 * i + a, 2 * i + a] += c
            Hm[2 * j + a, 2 * j + a] += c
            for b in range(2):
                v = -basis[i][a] @ basis[j][b]
                Hm[2 * i + a, 2 * j + b] += v
                Hm[2 * j + b, 2 * i + a] += v
    return np.linalg.eigvalsh(Hm)[0]


neel = [np.array([0, 0, 1.0]) * (-1) ** i for i in range(6)]
print("N=6 ring, Neel config: min Hess eigenvalue per |u|^2 =", round(hess_min(neel), 12))
best = 0
for _ in range(3000):
    ns = [v / np.linalg.norm(v) for v in rng.normal(size=(6, 3))]
    best = min(best, hess_min(ns))
print("N=6 ring, random configs: min found", round(best, 4), "(bound -4)")
print(f"\ntotal {time.time() - T0:.0f}s")
