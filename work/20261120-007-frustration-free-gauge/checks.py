"""Light numerical checks for work/20261120-007 (frustration-free SU(2) gauge Hamiltonian).

1. Single-link Witten-Laplacian gap gamma(kappa) on SU(2) = S^3:
     h(kappa) = E^2 + (kappa^2/16) sin^2(theta) - (3 kappa/8) cos(theta),   E^2 = -Laplacian_{S^3}/4,
   exact ground state exp(kappa cos(theta)/2) with energy 0.  Sector L of the SO(3) that fixes the
   staple direction: u = sin(theta) psi,
     -(1/4)(u'' + u) + L(L+1)/(4 sin^2 theta) u + V u,   u(0) = u(pi) = 0.
   gamma(kappa) = min(E_1(L=0), E_0(L=1)); higher L are higher (centrifugal term, checked up to L=3).
2. Random-matrix test of the two-projector inequality used in README section 4:
     A B + B A >= -c (A + B),  A = 1-P, B = 1-Q,  c = ||P Q - P^Q||.
Run: python tools/pin.py .venv/Scripts/python.exe work/20261120-007-frustration-free-gauge/checks.py
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh, null_space, orth

N = 4000
th = np.linspace(0, np.pi, N + 2)[1:-1]
h = th[1] - th[0]


def levels(kappa, L, k=2):
    V = kappa**2 / 16 * np.sin(th) ** 2 - 3 * kappa / 8 * np.cos(th)
    diag = 0.25 * (2 / h**2 - 1) + 0.25 * L * (L + 1) / np.sin(th) ** 2 + V
    off = -0.25 / h**2 * np.ones(N - 1)
    return eigh_tridiagonal(diag, off, select="i", select_range=(0, k - 1))[0]


print("# single-link gap; e0 is the L=0 ground energy (exact value 0) as discretization check")
print("kappa   e0(L=0)     E1(L=0)   E0(L=1)   E0(L=2)   E0(L=3)   gamma")
for kappa in [0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32]:
    l0 = levels(kappa, 0)
    l1, l2, l3 = (levels(kappa, L, 1)[0] for L in (1, 2, 3))
    print(f"{kappa:5.2f} {l0[0]: .2e} {l0[1]:9.5f} {l1:9.5f} {l2:9.5f} {l3:9.5f} {min(l0[1], l1):9.5f}")

rng = np.random.default_rng(0)
worst = np.inf
for trial in range(2000):
    n = rng.integers(3, 9)
    common = rng.integers(0, 2)  # 0 or 1 dimensional common kernel
    def rand_proj(r):
        X = rng.normal(size=(n, r)) + 1j * rng.normal(size=(n, r))
        Qm = orth(X)
        return Qm @ Qm.conj().T
    v = orth(rng.normal(size=(n, common)) + 1j * rng.normal(size=(n, common))) if common else np.zeros((n, 0))
    r1, r2 = rng.integers(0, n - common), rng.integers(0, n - common)
    def proj_containing(r):
        X = np.hstack([v, rng.normal(size=(n, r)) + 1j * rng.normal(size=(n, r))])
        Qm = orth(X)
        return Qm @ Qm.conj().T
    P, Q = proj_containing(r1), proj_containing(r2)
    # projector onto range(P) cap range(Q)
    M = np.vstack([np.eye(n) - P, np.eye(n) - Q])
    K = null_space(M, rcond=1e-9)
    PQ = K @ K.conj().T if K.size else np.zeros((n, n))
    c = np.linalg.norm(P @ Q - PQ, 2)
    A, B = np.eye(n) - P, np.eye(n) - Q
    T = A @ B + B @ A + c * (A + B)
    worst = min(worst, eigh((T + T.conj().T) / 2, eigvals_only=True)[0])
print(f"# projector inequality: min eigenvalue of AB+BA+c(A+B) over 2000 random cases = {worst:.2e} (>= -1e-10 expected)")
