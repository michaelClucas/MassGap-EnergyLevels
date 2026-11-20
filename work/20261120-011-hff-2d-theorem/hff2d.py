"""Gap of the frustration-free RK Hamiltonian H_FF(beta) in 2+1D (2D spatial lattice), free boundary.

gap(H_FF) = (1/4) * Poincare constant of the Wilson measure for the link Dirichlet form
E(f) = sum_l int |grad_l f|^2 dpi_beta  (unit S^3 per SU(2) link; unit circle per U(1) link).
All quantities below are for the gauge-invariant sector.

Sections
  A  single plaquette, SU(2): Sturm-Liouville on class functions (Ritz + finite differences)
  B  1x2 strip, SU(2): check of the reduced Dirichlet form E = 3|r1|^2+3|r2|^2+|l1-r2|^2 against
     direct differentiation of the 7-link configuration, then Ritz (upper bounds) on invariant polynomials
  C  L x L box, U(1): exact variational upper bound (lifted-flux test function)
  D  large box, SU(2): exact variational upper bound from nested concentric Wilson loops
Run: python tools/pin.py .venv/Scripts/python.exe work/20261120-011-hff-2d-theorem/hff2d.py
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.special import ive

rng = np.random.default_rng(11)
_leggauss = lru_cache(maxsize=None)(leggauss)


def ritz_min(K, C, tol=1e-11):
    """Smallest generalized eigenvalue of (K, C) on the numerically nonsingular part of C."""
    w, V = np.linalg.eigh(C)
    keep = w > tol * w.max()
    T = V[:, keep] / np.sqrt(w[keep])
    return np.linalg.eigvalsh(T.T @ K @ T)[0]


# ---------------- A: single plaquette ----------------
def plaquette_gap_ritz(beta, deg=40, nq=400):
    # class function f(theta); measure sin^2 e^{beta cos}; |grad f|^2 = f'(theta)^2 on unit S^3
    tmax = np.pi if beta < 30 else min(np.pi, 14 / np.sqrt(beta))
    x, w = leggauss(nq)
    t = 0.5 * tmax * (x + 1)
    w = 0.5 * tmax * w * np.sin(t) ** 2 * np.exp(beta * (np.cos(t) - 1))
    w /= w.sum()
    s = (1 - np.cos(t)) * max(beta, 1.0) / 4  # scaled variable, O(1) on the support
    P = np.polynomial.legendre.legvander(2 * s / s.max() - 1, deg)[:, 1:]
    dsdt = np.sin(t) * max(beta, 1.0) / 4 * 2 / s.max()
    dP = np.stack([np.polynomial.legendre.legval(2 * s / s.max() - 1,
                   np.polynomial.legendre.legder(np.eye(deg + 1)[k])) for k in range(1, deg + 1)], 1) * dsdt[:, None]
    Pc = P - (w @ P)
    C = Pc.T @ (w[:, None] * Pc)
    K = dP.T @ (w[:, None] * dP)
    lam1 = ritz_min(K, C)
    return lam1  # Poincare const of class functions; H_FF gap = (1/4)*4*lam1


def plaquette_gap_fd(beta, n=4000):
    tmax = np.pi if beta < 30 else min(np.pi, 14 / np.sqrt(beta))
    h = tmax / n
    tm = (np.arange(n) + 0.5) * h          # cell centres
    te = np.arange(n + 1) * h              # edges
    wc = np.sin(tm) ** 2 * np.exp(beta * (np.cos(tm) - 1))
    we = np.sin(te) ** 2 * np.exp(beta * (np.cos(te) - 1))
    # -(1/wc) d/dt (we df/dt), Neumann (flux zero at ends: we=0 at 0, truncation at tmax)
    main = np.zeros(n); off = np.zeros(n - 1)
    main[:-1] += we[1:-1] / h**2; main[1:] += we[1:-1] / h**2
    off[:] = -we[1:-1] / h**2
    from scipy.linalg import eigh_tridiagonal
    d = main / wc; e = off / np.sqrt(wc[:-1] * wc[1:])
    ev = eigh_tridiagonal(d, e, select='i', select_range=(0, 1), eigvals_only=True)
    return ev[1]


# ---------------- B: 1x2 strip ----------------
def qmul(a, b):
    a0, a1, a2, a3 = np.moveaxis(a, -1, 0); b0, b1, b2, b3 = np.moveaxis(b, -1, 0)
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)


def qinv(a):
    return a * np.array([1, -1, -1, -1])


def qexp(v):
    # v: imaginary quaternion as a 4-vector (v[0] = 0)
    th = np.linalg.norm(v)
    return np.array([1., 0, 0, 0]) * np.cos(th) + np.sin(th) * v / th if th > 0 else np.array([1., 0, 0, 0])


E_IMAG = np.eye(4)[1:]


def strip_W(U):
    b1, b2, t1, t2, v0, v1, v2 = U
    W1 = qmul(qmul(qmul(b1, v1), qinv(t1)), qinv(v0))
    loop2 = qmul(qmul(qmul(b2, v2), qinv(t2)), qinv(v1))
    W2 = qmul(qmul(b1, loop2), qinv(b1))
    return W1, W2


def check_strip_form(ntrial=5, eps=1e-6):
    """Compare sum_l |grad_l f|^2 by central differences with 3|r1|^2+3|r2|^2+|l1-r2|^2."""
    def f(W1, W2):
        s = W1[1:] @ W2[1:]
        return np.sin(1.3 * W1[0]) + W2[0] ** 3 + 0.7 * s + W1[0] * s ** 2
    def grads(W1, W2):
        def fx(x): return f(x, W2)
        def fy(y): return f(W1, y)
        r1 = np.array([(fx(qmul(W1, qexp(eps*e))) - fx(qmul(W1, qexp(-eps*e)))) / (2*eps) for e in E_IMAG])
        l1 = np.array([(fx(qmul(qexp(eps*e), W1)) - fx(qmul(qexp(-eps*e), W1))) / (2*eps) for e in E_IMAG])
        r2 = np.array([(fy(qmul(W2, qexp(eps*e))) - fy(qmul(W2, qexp(-eps*e)))) / (2*eps) for e in E_IMAG])
        l2 = np.array([(fy(qmul(qexp(eps*e), W2)) - fy(qmul(qexp(-eps*e), W2))) / (2*eps) for e in E_IMAG])
        return r1, l1, r2, l2
    worst = 0.0
    for _ in range(ntrial):
        U = rng.normal(size=(7, 4)); U /= np.linalg.norm(U, axis=1, keepdims=True)
        direct = 0.0
        for k in range(7):
            for e in E_IMAG:
                Up = U.copy(); Um = U.copy()
                Up[k] = qmul(U[k], qexp(eps*e)); Um[k] = qmul(U[k], qexp(-eps*e))
                direct += ((f(*strip_W(Up)) - f(*strip_W(Um))) / (2*eps)) ** 2
        r1, l1, r2, l2 = grads(*strip_W(U))
        formula = 3*r1@r1 + 3*r2@r2 + (l1-r2)@(l1-r2)
        worst = max(worst, abs(direct - formula) / abs(direct), np.abs(l1 + l2 - r1 - r2).max())
    return worst


def strip_gap_ritz(beta, D=7, n1=40, nphi=20):
    """Ritz upper bound on the Poincare constant of the 1x2 strip (gauge-invariant sector).
    Basis: monomials u^a v^b s^c (a+b+c<=D), u=beta(1-x0), v=beta(1-y0), s=beta x.y (vector parts)."""
    tmax = np.pi if beta < 30 else min(np.pi, 14 / np.sqrt(beta))
    xt, wt = leggauss(n1); t = 0.5*tmax*(xt+1); wt = 0.5*tmax*wt
    xp, wp = leggauss(nphi); ph = 0.5*np.pi*(xp+1); wp = 0.5*np.pi*wp
    T1, T2, PH = np.meshgrid(t, t, ph, indexing='ij')
    Wt = (wt[:, None, None]*wt[None, :, None]*wp[None, None, :]
          * np.sin(T1)**2*np.sin(T2)**2*np.sin(PH)*np.exp(beta*(np.cos(T1)+np.cos(T2)-2))).ravel()
    Wt /= Wt.sum()
    T1, T2, PH = T1.ravel(), T2.ravel(), PH.ravel()
    x = np.stack([np.cos(T1), 0*T1, 0*T1, np.sin(T1)], -1)
    y = np.stack([np.cos(T2), np.sin(T2)*np.sin(PH), 0*T2, np.sin(T2)*np.cos(PH)], -1)
    bs = max(beta, 1.0)
    u = bs*(1-x[:, 0]); v = bs*(1-y[:, 0]); s = bs*np.einsum('ij,ij->i', x[:, 1:], y[:, 1:])
    exps = [(a, b, c) for a in range(D+1) for b in range(D+1-a) for c in range(D+1-a-b) if a+b+c > 0]
    nb = len(exps)
    # tangent directions: x e_a (right), e_a x (left), y e_a (right)
    xr = np.stack([qmul(x, np.broadcast_to(e, x.shape)) for e in E_IMAG], 1)   # (N,3,4)
    xl = np.stack([qmul(np.broadcast_to(e, x.shape), x) for e in E_IMAG], 1)
    yr = np.stack([qmul(y, np.broadcast_to(e, y.shape)) for e in E_IMAG], 1)
    # d(u,v,s)/dx and /dy in R^4
    zero = np.zeros_like(x)
    du_dx = zero.copy(); du_dx[:, 0] = -bs
    ds_dx = zero.copy(); ds_dx[:, 1:] = bs*y[:, 1:]
    dv_dy = zero.copy(); dv_dy[:, 0] = -bs
    ds_dy = zero.copy(); ds_dy[:, 1:] = bs*x[:, 1:]
    proj = lambda g, tang: np.einsum('nk,nak->na', g, tang)
    U_r1, S_r1 = proj(du_dx, xr), proj(ds_dx, xr)
    U_l1, S_l1 = proj(du_dx, xl), proj(ds_dx, xl)
    V_r2, S_r2 = proj(dv_dy, yr), proj(ds_dy, yr)
    F = np.empty((len(Wt), nb)); r1 = np.empty((len(Wt), nb, 3)); l1 = np.empty_like(r1); r2 = np.empty_like(r1)
    pw = lambda z, k: z**k if k >= 0 else 0*z
    for i, (a, b, c) in enumerate(exps):
        F[:, i] = pw(u, a)*pw(v, b)*pw(s, c)
        fu = a*pw(u, a-1)*pw(v, b)*pw(s, c); fv = b*pw(u, a)*pw(v, b-1)*pw(s, c); fs = c*pw(u, a)*pw(v, b)*pw(s, c-1)
        r1[:, i] = fu[:, None]*U_r1 + fs[:, None]*S_r1
        l1[:, i] = fu[:, None]*U_l1 + fs[:, None]*S_l1
        r2[:, i] = fv[:, None]*V_r2 + fs[:, None]*S_r2
    Fc = F - Wt @ F
    C = Fc.T @ (Wt[:, None]*Fc)
    K = np.zeros((nb, nb))
    for comp in range(3):
        for A, cf in ((r1[:, :, comp], 3.0), (r2[:, :, comp], 3.0), (l1[:, :, comp]-r2[:, :, comp], 1.0)):
            K += cf*(A.T @ (Wt[:, None]*A))
    return ritz_min(K, C)


# ---------------- C: U(1), L x L box ----------------
def u1_bound(beta, L, n=200001):
    F = np.linspace(-np.pi, np.pi, n)[:-1]; h = F[1]-F[0]
    logrho = beta*np.cos(F) - np.log(2*np.pi*ive(0, beta)) - beta
    rho = np.exp(logrho)
    lam_M = 4 - 4*np.cos(np.pi/(L+1))
    out = {}
    # (i) Kramers-optimal lifted function: g' = 1 - lam/rho, lam = 2pi/int(1/rho)
    inv_int = np.sum(np.exp(-logrho))*h
    lam = 2*np.pi/inv_int
    gp = 1 - lam*np.exp(-logrho)
    g = np.concatenate([[0], np.cumsum(0.5*(gp[1:]+gp[:-1])*h)])
    Eg = np.sum(rho*g)*h; Varg = np.sum(rho*(g-Eg)**2)*h
    Egp = np.sum(rho*gp)*h; Vargp = np.sum(rho*(gp-Egp)**2)*h
    out['lifted'] = 0.25*(Egp**2*lam_M + 4*Vargp)/Varg
    out['lifted_Linf'] = 0.25*4*Vargp/Varg
    # (ii) g = sin F for comparison
    Es = np.sum(rho*np.sin(F))*h; Vs = np.sum(rho*np.sin(F)**2)*h - Es**2
    Ec = np.sum(rho*np.cos(F))*h; Vc = np.sum(rho*np.cos(F)**2)*h - Ec**2
    out['sin'] = 0.25*(Ec**2*lam_M + 4*Vc)/Vs
    out['lam_M'] = lam_M
    return out


# ---------------- D: SU(2), nested Wilson loops ----------------



@lru_cache(maxsize=None)
def su2_moments(beta, nq=800):
    tmax = np.pi if beta < 30 else min(np.pi, 14/np.sqrt(beta))
    x, w = _leggauss(nq); t = 0.5*tmax*(x+1)
    w = 0.5*tmax*w*np.sin(t)**2*np.exp(beta*(np.cos(t)-1)); w /= w.sum()
    a = w @ np.cos(t)                      # c_{1/2} = E cos theta
    c1 = (4*(w @ np.cos(t)**2) - 1)/3      # c_1 = E chi_1 / 3
    return a, c1


def su2_nested_bound(beta, smax):
    """Concentric squares of odd side s=1,3,..,smax (disjoint boundaries); f = sum_k a_k cos theta_k.
    Returns (1/4) min Rayleigh quotient = rigorous upper bound on gap(H_FF) for any box containing them."""
    a, c1 = su2_moments(beta)
    sides = np.arange(1, smax+1, 2); n = sides**2
    E2 = (1 + 3*c1**n)/4
    D = np.diag(4*sides*(1 - E2))
    N = len(sides); C = np.empty((N, N))
    for i in range(N):
        for j in range(N):
            k, m = min(i, j), max(i, j)
            C[i, j] = a**(n[m]-n[k])*E2[k] - a**(n[k]+n[m])
    single = min(D[i, i]/C[i, i] for i in range(N))
    return 0.25*ritz_min(D, C, tol=1e-14), 0.25*single


if __name__ == '__main__':
    import sys, time
    secs = sys.argv[1] if len(sys.argv) > 1 else 'ABCD'
    t0 = time.time()
    if 'A' in secs:
        print('# A: single plaquette (gauge-invariant), gap(H_FF) = lambda_1(class fns); compare 2*beta', flush=True)
        for beta in [0.0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]:
            print(f'beta={beta:7.2f}  ritz={plaquette_gap_ritz(beta):.8f}  fd={plaquette_gap_fd(beta):.6f}  2beta={2*beta:.1f}', flush=True)
    if 'C' in secs:
        print('# C: U(1), L x L box: upper bounds on gap(H_FF); compare 1/beta', flush=True)
        for beta in [1, 2, 4, 6, 8, 10, 15, 20]:
            rows = [u1_bound(beta, L) for L in (4, 16, 64, 10**4)]
            print(f'beta={beta:3d}  lifted L=4,16,64,1e4: ' + ' '.join(f'{r["lifted"]:.3e}' for r in rows)
                  + f'  L=inf {rows[0]["lifted_Linf"]:.3e}  sin(L=1e4) {rows[-1]["sin"]:.3e}  1/beta {1/beta:.3e}', flush=True)
    if 'D' in secs:
        print('# D: SU(2) nested concentric Wilson loops, upper bound on gap(H_FF) uniform over boxes of side >= smax', flush=True)
        for beta in [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]:
            res = [su2_nested_bound(beta, smax) for smax in (1, 9, 41, 201)]
            print(f'beta={beta:7.1f}  smax=1,9,41,201: ' + ' '.join(f'{r[0]:.4f}' for r in res)
                  + f'  best single loop {res[-1][1]:.4f}  1/beta {1/beta:.4f}', flush=True)
    if 'B' in secs:
        print('# B: strip reduced-form check, max rel. deviation / conj-invariance residual:', f'{check_strip_form():.2e}', flush=True)
        print('# B: 1x2 strip Ritz upper bound on gap(H_FF)=(1/4)lambda; linearized large-beta value 1.5*beta', flush=True)
        for beta in [0.0, 1, 5, 20, 50]:
            g4 = 0.25*strip_gap_ritz(beta, D=4, n1=30, nphi=16); g6 = 0.25*strip_gap_ritz(beta, D=6, n1=30, nphi=16)
            print(f'beta={beta:6.2f}  D4={g4:.6f}  D6={g6:.6f}  1.5beta={1.5*beta:.2f}  [{time.time()-t0:.0f}s]', flush=True)
