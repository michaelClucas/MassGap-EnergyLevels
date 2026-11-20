"""Line 5 (interaction-driven gaps): explicit calculations for README.md.

Sections (run all by default, or pass letters, e.g. `ABD`):
  A  Large-N (spherical) gap equation on Z^d, d = 2, 3, 4; applied to N = 2 (2D XY).
  B  Fermi-surface log (BCS) vs bosonic Hartree gap equation in d = 4.
  C  Coupling at the gap scale for an RG-invariant gap m = kappa * Lambda (one loop).
  D  Bosonic Gaussian (Hartree / BCS-type) variational gap vs exact gap, 1D anharmonic QM.
  E  Gauge-parameter dependence of the one-loop Hartree (tadpole) gluon mass.
  F  Positivity violation of the Gribov propagator (zero-momentum Schwinger function).
  G  Peierls: fermions vs bosons in a dimerized ring; Landau levels of spin-1 vs Dirac.

No measured numbers are used (none needed); the only external number quoted is
the 2D XY transition point in the README, cited there, not used here.
Run: python tools/pin.py .venv/Scripts/python.exe -u work/20261120-014-line5-interaction-gaps/interaction_gaps.py
"""
import sys

import numpy as np
from scipy import integrate, optimize, special
from scipy.linalg import eigh_tridiagonal


# ---------------------------------------------------------------- A
def W(d, m):
    """W_d(m) = int_{[-pi,pi]^d} d^dk/(2pi)^d 1/(m^2 + sum_i (2 - 2cos k_i)).

    Uses 1/a = int_0^inf e^{-a t} dt and int dk/2pi e^{-t(2-2cos k)} = e^{-2t} I0(2t),
    with t = e^u to handle the slow tail at small m."""
    def f(u):
        t = np.exp(u)
        if t < 1e6:
            b = special.ive(0, 2 * t)
        else:  # scipy ive returns nan near 2t ~ 1e9; asymptotic series, rel. error < 1e-17 here
            b = (1 + 1 / (16 * t) + 9 / (512 * t * t)) / np.sqrt(4 * np.pi * t)
        return t * np.exp(-m * m * t) * b ** d
    umax = np.log(60.0 / max(m * m, 1e-300)) if m > 0 else 60.0
    val, err = integrate.quad(f, -40, umax, limit=800, epsabs=1e-12, epsrel=1e-10)
    return val


def solve_gap(d, beta):
    """Spherical-model (N = inf O(N)) gap equation beta = W_d(m); None if no solution."""
    if d >= 3 and beta >= W(d, 0.0):
        return None
    g = lambda lm: W(d, np.exp(lm)) - beta
    return float(np.exp(optimize.brentq(g, -40, 5, xtol=1e-12)))


def section_A():
    print("== A. Large-N gap equation beta = W_d(m) on Z^d (spherical model) ==")
    for d in (3, 4):
        print(f"d = {d}: W_d(0) = {W(d, 0.0):.6f}  -> no massive solution for beta >= W_d(0) (gapless/ordered phase)")
    print("d = 2: W_2(m) ~ (1/4pi) ln(32/m^2) diverges as m -> 0: massive solution for every beta")
    print(f"{'beta':>6} {'m (gap eq)':>14} {'sqrt(32) exp(-2 pi beta)':>26}")
    for beta in (0.3, 0.5, 0.8, 1.0, 1.5):
        m = solve_gap(2, beta)
        print(f"{beta:6.2f} {m:14.6e} {np.sqrt(32) * np.exp(-2 * np.pi * beta):26.6e}")
    print("Applied to N = 2 (2D XY, s in S^1, weight exp(beta_XY sum s.s')): beta_sph = beta_XY / 2")
    print(f"{'beta_XY':>8} {'m_gapeq':>12} {'xi_gapeq = 1/m':>16}")
    for bxy in (0.6, 0.9, 1.2, 1.5, 2.0, 3.0):
        m = solve_gap(2, bxy / 2)
        print(f"{bxy:8.2f} {m:12.4e} {1 / m:16.4e}")
    print("Leading-order gap equation predicts m > 0 at every beta_XY; Froehlich-Spencer 1981 prove")
    print("power-law decay (m = 0) for beta_XY large. Gap equation predicts a gap where none exists.\n")


# ---------------------------------------------------------------- B
def section_B():
    print("== B. BCS log vs bosonic Hartree ==")
    print("BCS (constant DOS, cutoff w): 1 = lam * asinh(w/Delta)  =>  Delta/w = 1/sinh(1/lam)")
    print("Boson Hartree, d = 4 lattice (cutoff = 1/a): m^2 = lam * W_4(m), massless bare")
    print(f"{'lam':>6} {'Delta/w (BCS)':>14} {'m (boson d=4)':>14} {'sqrt(lam W4(0))':>16}")
    W40 = W(4, 0.0)
    for lam in (0.05, 0.1, 0.2, 0.5, 1.0):
        dbcs = 1 / np.sinh(1 / lam)
        m = optimize.brentq(lambda m: m * m - lam * W(4, m), 1e-8, 10)
        print(f"{lam:6.2f} {dbcs:14.4e} {m:14.4e} {np.sqrt(lam * W40):16.4e}")
    print("BCS: Delta ~ 2w exp(-1/lam): nonperturbative in lam, RG-invariant scale, from the 1D (Fermi-surface) log.")
    print("Boson d = 4: m ~ sqrt(lam) * cutoff: perturbative tadpole, removed by a mass counterterm;")
    print("in a gauge theory it is forbidden by BRST and vanishes in dim. reg.\n")


# ---------------------------------------------------------------- C
def section_C():
    print("== C. Coupling at the gap scale, one loop, SU(N) pure YM ==")
    print("alpha(mu) = 4 pi / (beta0 ln(mu^2/Lambda^2)), beta0 = 11N/3; a = alpha N /(4 pi) = 3/(11 ln kappa^2)")
    print(f"{'kappa=m/Lambda':>15} {'a(mu=m)':>10}")
    for kappa in (2, 4, 7, 15, 100, 1e4):
        print(f"{kappa:15.4g} {3 / (11 * np.log(kappa ** 2)):10.4f}")
    print("kappa is a pure number fixed by the theory, not a tunable parameter (contrast BCS: lam free).\n")


# ---------------------------------------------------------------- D
def exact_levels(mu, lam, n=6000, L=None):
    if L is None:
        L = 8.0 + 4.0 * np.sqrt(max(-mu, 0) / max(lam, 1e-12))
    x = np.linspace(-L, L, n)
    h = x[1] - x[0]
    V = 0.5 * mu * x ** 2 + lam * x ** 4
    diag = 1.0 / h ** 2 + V
    off = -0.5 / h ** 2 * np.ones(n - 1)
    return eigh_tridiagonal(diag, off, select="i", select_range=(0, 2))[0]


def gaussian_mf(mu, lam):
    """Gaussian trial state centered at x0, width s^2 = 1/(2w) (squeezed coherent state:
    the bosonic analogue of a BCS / Hartree-Fock-Bogoliubov vacuum)."""
    def E(p):
        x0, lw = p
        w = np.exp(lw)
        s2 = 1 / (2 * w)
        return w / 4 + 0.5 * mu * (x0 ** 2 + s2) + lam * (x0 ** 4 + 6 * x0 ** 2 * s2 + 3 * s2 ** 2)
    sym = optimize.minimize_scalar(lambda lw: E((0.0, lw)), bounds=(-10, 10), method="bounded")
    w_sym = np.exp(sym.x)
    best = sym.fun, 0.0, w_sym
    for x0 in (0.5, 1.0, 2.0, 4.0):
        r = optimize.minimize(E, (x0, 0.0), method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-14, maxiter=5000))
        if r.fun < best[0] - 1e-10:
            best = r.fun, abs(r.x[0]), np.exp(r.x[1])
    # Hartree frequency at the optimum: w^2 = <V''> (stationarity), = MF excitation energy
    return w_sym, sym.fun, best


def section_D():
    print("== D. Gaussian (Hartree/HFB) variational gap vs exact gap, H = p^2/2 + mu x^2/2 + lam x^4 ==")
    print(f"{'mu':>5} {'lam':>5} {'E0 exact':>10} {'E0 MF':>10} {'gap exact':>11} {'gap MF(sym)':>12} {'MF broken?':>10} {'gap MF(brk)':>12} {'MF/exact':>9} {'symMF/exact':>12}")
    rows = [(1, 0.1), (1, 1), (0, 1), (-1, 0.1), (-2, 0.1), (-3, 0.1), (-4, 0.1)]
    for mu, lam in rows:
        e = exact_levels(mu, lam)
        w_sym, e_sym, (e_b, x0, w_b) = gaussian_mf(mu, lam)
        gap = e[1] - e[0]
        broken = x0 > 1e-4
        g_mf = 0.0 if broken else w_sym   # broken MF: two degenerate vacua, predicted splitting 0
        print(f"{mu:5.1f} {lam:5.2f} {e[0]:10.5f} {min(e_sym, e_b):10.5f} {gap:11.4e} {w_sym:12.4f} {str(broken):>10} {g_mf:12.4f} {(g_mf / gap):9.3g} {(w_sym / gap):12.4g}")
    print("E0 MF >= E0 exact always (Rayleigh-Ritz). Symmetric MF gap/exact ranges 1.01 to ~1e6 (no constant);")
    print("broken MF predicts 0 where exact > 0. Neither MF gap gives two-sided control of E1 - E0.\n")


# ---------------------------------------------------------------- E
def section_E():
    print("== E. Covariant-gauge one-loop tadpole (4-gluon vertex) at zero momentum ==")
    print("D_{rs}(k) = [(delta - k k/k^2) + xi k k/k^2]/k^2; angular average of k_r k_s/k^2 = delta_rs/d")
    d = 4
    rng = np.random.default_rng(0)
    k = rng.normal(size=(200000, d))
    kk = np.einsum("ni,nj->ij", k / np.linalg.norm(k, axis=1)[:, None], k / np.linalg.norm(k, axis=1)[:, None]) / len(k)
    print(f"{'xi':>5} {'Pi_00 / (g^2 C_A I)':>22} {'(d-1+xi)(d-1)/d':>18}")
    for xi in (0.0, 0.5, 1.0, 3.0):
        Dbar = (np.eye(d) - kk) + xi * kk        # times I = int d^dk/(2pi)^d 1/k^2
        Pi = np.eye(d) * np.trace(Dbar) - Dbar   # structure g_{mn} D_rr - D_{mn}
        print(f"{xi:5.1f} {Pi[0, 0]:22.4f} {(d - 1 + xi) * (d - 1) / d:18.4f}")
    print("Hartree mass^2 from the quartic vertex is linear in xi: Feynman/Landau = 4/3 in d = 4.")
    print("A physical pole mass is xi-independent (Nielsen identities), so this 'mass' is not one.\n")


# ---------------------------------------------------------------- F
def section_F():
    print("== F. Gribov propagator D(p) = p^2/(p^4 + gamma^4), Schwinger function at zero spatial momentum ==")
    g = 1.0
    print(f"{'t*gamma':>8} {'C(t) numeric':>14} {'C(t) closed form':>17}")
    for t in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
        num = integrate.quad(lambda p: p * p / (p ** 4 + g ** 4) / np.pi, 0, np.inf, weight="cos", wvar=t)[0]
        a = g * t / np.sqrt(2)
        cf = g / (2 * np.sqrt(2)) * np.exp(-a) * (np.cos(a) - np.sin(a))
        print(f"{t:8.2f} {num:14.6e} {cf:17.6e}")
    print("C(t) < 0 for gamma t/sqrt2 in (pi/4, 5pi/4): no positive Kallen-Lehmann measure, so gamma is not")
    print("an energy gap of any positive-metric Hamiltonian.\n")


# ---------------------------------------------------------------- G
def section_G():
    print("== G1. Peierls: dimerized ring, hoppings t(1 +- delta), per-site band energy, t = 1 ==")
    M = 200000
    kgrid = 2 * np.pi * (np.arange(M) + 0.5) / M
    def e_ferm(dl):   # spinless, half filling: lower band full, one particle per 2-site cell
        t1, t2 = 1 + dl, 1 - dl
        return -np.mean(np.sqrt(t1 * t1 + t2 * t2 + 2 * t1 * t2 * np.cos(kgrid))) / 2
    def e_bos(dl):    # n bosons per site, all in the band bottom: -(t1 + t2) per particle
        return -(2.0)
    e0 = e_ferm(0)
    print(f"{'delta':>7} {'fermion gain/site':>18} {'delta^2 ln(1/delta)':>20} {'boson gain/particle':>20}")
    for dl in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005):
        print(f"{dl:7.3f} {e0 - e_ferm(dl):18.6e} {dl * dl * np.log(1 / dl):20.6e} {e_bos(0) - e_bos(dl):20.6e}")
    ds = np.array([0.02, 0.01, 0.005, 0.0025, 0.00125])
    gains = np.array([e0 - e_ferm(x) for x in ds])
    A = np.vstack([ds ** 2 * np.log(1 / ds), ds ** 2]).T
    a, b = np.linalg.lstsq(A, gains, rcond=None)[0]
    print(f"fit gain = a delta^2 ln(1/delta) + b delta^2, delta in [0.00125, 0.02]: a = {a:.5f}, b = {b:.5f}  (1/pi = {1 / np.pi:.5f})")
    print("Fermions: gain ~ c delta^2 ln(1/delta) beats elastic cost K delta^2/2 for every K: dimerization at")
    print("any coupling (Kennedy-Lieb). Free bosons: gain exactly 0, so no instability at any K.\n")

    print("== G2. Landau levels in field B (charge 1), lowest E^2 at k_z = 0 ==")
    print("spin-1 charged vector (W boson / off-diagonal gluon, g-factor 2): E^2 = (2n+1)B - 2 B s_z, s_z = -1, 0(unphys.), +1")
    print("Dirac electron (g = 2):                                         E^2 = (2n+1)B - 2 B s_z, s_z = +-1/2")
    for n in range(3):
        ev = [(2 * n + 1) - 2 * s for s in (-1, 1)]
        ed = [(2 * n + 1) - 2 * s for s in (-0.5, 0.5)]
        print(f"n={n}: vector E^2/B = {ev}   Dirac E^2/B = {ed}")
    print("Vector n=0, s_z=+1: E^2 = -B < 0 (Nielsen-Olesen tachyon). Dirac minimum is 0, never negative.\n")


SECTIONS = dict(A=section_A, B=section_B, C=section_C, D=section_D, E=section_E, F=section_F, G=section_G)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCDEFG"
    for s in which:
        SECTIONS[s]()
