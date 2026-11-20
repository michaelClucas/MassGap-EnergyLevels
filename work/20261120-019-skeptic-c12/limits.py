"""Skeptic 019: test the 013 model implementation in limits with known answers.

Run from repo root: python tools/pin.py .venv/Scripts/python.exe work/20261120-019-skeptic-c12/limits.py
Imports the 013 model.py unchanged (copied here as model.py; identical file).
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import erf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as M  # noqa: E402

ok = True


def check(name, got, want, tol):
    global ok
    good = abs(got - want) <= tol
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}: got {got:.6f}, want {want:.6f} (tol {tol:g})")


g = M.Grid(R=14.0, n=700)
# Airy zeros |a_1|, and l=1, l=2 lowest eigenvalues of -u'' + l(l+1)/r^2 u + r u (KS Table 1 also lists 3.3613, 4.2482)
A1, L1, L2 = 2.338107410, 3.361254566, 4.248182160

# 1. NR: H = 2m + p^2/m + s r  ->  E = 2m + (s^2/m)^(1/3) eps_l
for m, s in ((1.0, 1.0), (3.1323, 9 / 4)):
    for cent, eps in ((0, A1), (2, L1), (6, L2)):
        e, _ = M.solve(g, cent, s * g.r, m, 2, "nr", nlev=1)
        check(f"nr linear m={m} s={s} l(l+1)={cent}", e[0], 2 * m + (s * s / m) ** (1 / 3) * eps, 2e-3)

# 2. Relativistic: 2 sqrt(p^2+m^2) + r, L = 0, against an independent sine-basis (exact |p| in the box) solver
def sine_basis(m, s, R=14.0, nb=500, nq=6000):
    k = np.arange(1, nb + 1) * math.pi / R
    r = (np.arange(nq) + 0.5) * R / nq
    B = np.sqrt(2 / R) * np.sin(np.outer(k, r))
    V = (B * (s * r)) @ B.T * (R / nq)
    H = np.diag(2 * np.sqrt(k ** 2 + m ** 2)) + V
    return np.linalg.eigvalsh(H)[:2]


for m in (0.0, 1.0, 3.1323):
    e, _ = M.solve(g, 0, g.r, m, 2, "rel", nlev=2)
    ref = sine_basis(m, 1.0)
    check(f"rel linear m={m} n=0 vs sine basis", e[0], ref[0], 5e-3)
    check(f"rel linear m={m} n=1 vs sine basis", e[1], ref[1], 5e-3)
# massless: H = 2|p| + r = sqrt(2) (|p'| + r'); literature E0(|p| + r, l=0) = 2.2322 (Lucha-Schoeberl, from memory)
e, _ = M.solve(g, 0, g.r, 0.0, 2, "rel", nlev=1)
check("rel massless 2|p|+r vs sqrt2*2.2322 (literature value quoted from memory)", e[0], math.sqrt(2) * 2.2322, 5e-3)

# 3. Tensor matrix elements: S=0 gives 0; sum_J (2J+1) <T>_{LSJ} = 0; S=1 J-pattern proportional to -4 : 2 : -2/5
for L, S in ((1, 1), (2, 2), (2, 0), (3, 1)):
    tot = sum((2 * J + 1) * M.tensor_me(L, S, J) for J in range(abs(L - S), L + S + 1))
    check(f"tensor trace L={L} S={S}", tot, 0.0, 1e-8)
t = [M.tensor_me(1, 1, J) for J in (0, 1, 2)]
check("tensor 3P J=1/J=0 ratio", t[1] / t[0], 2 / -4, 1e-8)
check("tensor 3P J=2/J=0 ratio", t[2] / t[0], (-2 / 5) / -4, 1e-8)
print(f"      tensor L=1 S=1 J=0,1,2: {[round(x, 5) for x in t]};  L=2 S=2 J=0..4: "
      f"{[round(M.tensor_me(2, 2, J), 5) for J in range(5)]}")
# independent value for L=1,S=1,J=0: for spin-1 x spin-1 -> S=1 the rank-2 spin tensor has the same reduced
# structure as spin-1/2 x spin-1/2 -> S=1 up to a factor; compute that factor directly from operators below.
Sop = M.spin1_ops()
I3 = np.eye(3)
S1 = [np.kron(a, I3) for a in Sop]
S2 = [np.kron(I3, a) for a in Sop]
Stot = [a + b for a, b in zip(S1, S2)]
S2tot = sum(a @ a for a in Stot)
w, U = np.linalg.eigh(S2tot)
P1 = U[:, np.isclose(w, 2)]  # S = 1 subspace
# <S=1, mS=1| 3 S1z S2z - S1.S2 |S=1, mS=1> for spin-1 pair; compare with spin-1/2 pair value (3/4 - 1/4 = 1/2)
Sz = P1.conj().T @ Stot[2] @ P1
wz, Uz = np.linalg.eigh(Sz)
v = P1 @ Uz[:, np.argmax(wz)]
T_zz = np.real(v.conj() @ (3 * S1[2] @ S2[2] - sum(a @ b for a, b in zip(S1, S2))) @ v)
print(f"      spin-1 pair S=1 mS=1: <3 S1z S2z - S1.S2> = {T_zz:.4f} (spin-1/2 pair: +0.5) -> "
      f"tensor force sign {'same' if T_zz > 0 else 'opposite'} to quarks for S=1")
# For spin-1/2, <3P0| 3(S1.n)(S2.n) - S1.S2 |3P0> = -1 (Pauli S12 = -4); scaled by T_zz/0.5
check("tensor 3P0 absolute (spin-1 pair, scaled from spin-1/2 result)", t[0], -1.0 * T_zz / 0.5, 1e-6)

# 4. smeared OGE: derivatives and Laplacian identity grad^2 V_V = 4 pi kappa delta_s
kappa, rs = 0.7187, 1 / 3.1323
r = np.linspace(0.05, 3, 3001)
x = r / rs
VV = -kappa * erf(x) / r
dVV = kappa * (erf(x) / r ** 2 - (2 / math.sqrt(math.pi)) * np.exp(-x ** 2) / (rs * r))
d2VV = kappa * ((4 / math.sqrt(math.pi)) * np.exp(-x ** 2) / (rs * r ** 2) - 2 * erf(x) / r ** 3
                + (4 / math.sqrt(math.pi)) * np.exp(-x ** 2) / rs ** 3)
fd1 = np.gradient(VV, r)
fd2 = np.gradient(fd1, r)
check("max |V_V' - finite diff|", float(np.max(np.abs(dVV - fd1)[5:-5])), 0.0, 1e-3)
check("max |V_V'' - finite diff|", float(np.max(np.abs(d2VV - fd2)[10:-10])), 0.0, 5e-3)
delta = np.exp(-x ** 2) / (math.pi ** 1.5 * rs ** 3)
check("max |grad^2 V_V - 4 pi kappa delta_s|", float(np.max(np.abs(d2VV + 2 * dVV / r - 4 * math.pi * kappa * delta))),
      0.0, 1e-8)
rr = np.linspace(0, 10 * rs, 20001)
check("int 4 pi r^2 delta_s", float(np.trapezoid(4 * math.pi * rr ** 2 * np.exp(-(rr / rs) ** 2), rr)
                                    / (math.pi ** 1.5 * rs ** 3)), 1.0, 1e-6)

# 5. colour factors
check("C_A/C_F, N=3", M.casimir_ratio(3), 9 / 4, 1e-12)
check("C_A/C_F, N=1e6", M.casimir_ratio(10 ** 6), 2.0, 1e-9)
# S1.S2 for two spin-1: S=0 -> -2, S=2 -> +1; three spin-1, S=1 -> -2, S=3 -> +3
for S, n, want in ((0, 2, -2), (1, 2, -1), (2, 2, 1), (1, 3, -2), (3, 3, 3)):
    check(f"sum S_i.S_j  S={S} n={n}", M.s1s2(S, n), want, 1e-12)

# 6. hyperspherical K=0 (6D) averages
a = np.linspace(0, math.pi / 2, 200001)
wgt = np.cos(a) ** 2 * np.sin(a) ** 2
Z = np.trapezoid(wgt, a)
check("<cos a> K=0 6D", float(np.trapezoid(wgt * np.cos(a), a) / Z), 32 / (15 * math.pi), 1e-8)
check("<1/cos a> K=0 6D", float(np.trapezoid(wgt / np.maximum(np.cos(a), 1e-300), a) / Z), 16 / (3 * math.pi), 1e-6)

print("\nALL PASS" if ok else "\nSOME FAILED")
