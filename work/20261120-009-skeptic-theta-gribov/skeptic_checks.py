"""Skeptic checks on C4 (theta/Bloch) and C7 (Payne-Weinberger on the Gribov region).

Reuses the 008 model code via an unmodified copy in rerun/ (rerun/gribov_scope.py, rerun/theta_bands.py).
Writes skeptic_summary.txt. Usage (repo root):
  python tools/pin.py .venv/Scripts/python.exe work/20261120-009-skeptic-theta-gribov/skeptic_checks.py
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / "rerun"))
import gribov_scope as G  # noqa: E402
import theta_bands as TB  # noqa: E402

out = []


def say(s=""):
    print(s, flush=True)
    out.append(s)


# ---------------- C4 -----------------
def pendulum_full(V0, theta, mmax=60):
    m = np.arange(-mmax, mmax + 1)
    q = (m + theta / (2 * np.pi))
    H = np.diag(0.5 * q**2 + V0) + np.diag(np.full(len(m) - 1, -V0 / 2), 1) + np.diag(np.full(len(m) - 1, -V0 / 2), -1)
    w, v = np.linalg.eigh(H)
    return w, v, q / (2 * np.pi)  # H' = dH/dtheta is diagonal q/(2 pi); H'' = 1/(4 pi^2)


def c4_checks():
    say("== C4-a: N=2 pendula, Bloch fiber of T_total = T1 T2 only (theta on total winding) ==")
    say("Fiber T_total = e^{i theta} decomposes further over T1 = e^{i phi}, phi continuous;")
    say("its spectrum is the union over phi of e_a(phi) + e_b(theta - phi).")
    phis = np.linspace(0, 2 * np.pi, 401)
    for V0 in [1.0, 4.0]:
        e = np.array([TB.bands(V0, p, k=2) for p in phis])  # e[i, n]
        for th in [0.0, np.pi]:
            j = np.round(((th - phis) % (2 * np.pi)) / (2 * np.pi) * 400).astype(int) % 400
            s00 = e[:, 0] + e[j, 0]
            s01 = e[:, 0] + e[j, 1]
            say(f" V0={V0} theta={th:.3f}: lowest band of fiber = [{s00.min():.8f}, {s00.max():.8f}] "
                f"(width {s00.max() - s00.min():.3e}, continuous => no isolated ground state, in-sector gap 0); "
                f"next band starts {s01.min():.6f}")
    say(" Per-copy fiber (T1 = T2 = e^{i theta}, the model used in 008 A4) instead has gap = single-well gap.")
    say("")
    say("== C4-b: theta-curvature sum rule gives only an UPPER bound on the fiber gap ==")
    say(" E0'' = <H''> - 2 sum_n |<n|H'|0>|^2/(E_n-E_0) >= <H''> - 2 Var0(H')/Delta  =>  Delta <= 2 Var0(H')/(<H''> - E0'')")
    for V0 in [0.05, 1.0, 4.0, 8.0]:
        for th in [0.0, np.pi / 2]:
            w, v, Hp = pendulum_full(V0, th)
            psi = v[:, 0]
            mean = np.sum(Hp * psi**2)
            var = np.sum(Hp**2 * psi**2) - mean**2
            Hpp = 1 / (4 * np.pi**2)
            me = v.T @ (Hp * psi)
            E0pp = Hpp + 2 * np.sum(me[1:] ** 2 / (w[0] - w[1:]))
            ub = 2 * var / (Hpp - E0pp)
            say(f" V0={V0:5.2f} theta={th:.3f}: Delta={w[1] - w[0]:.5f}  sum-rule upper bound={ub:.5f}  "
                f"holds={w[1] - w[0] <= ub + 1e-12}")
    say("")


# ---------------- C7 -----------------
def c7_pw_nonmonotone():
    say("== C7-a: adding V >= 0 can lower or raise the gap (Neumann -d^2/dx^2 on [-1,1], D = 2) ==")
    n = 2000
    h = 2.0 / n
    x = -1 + h * (np.arange(n) + 0.5)
    d0 = np.full(n, 2.0 / h**2); d0[0] = d0[-1] = 1.0 / h**2
    off = np.full(n - 1, -1.0 / h**2)
    for name, V in [("V=0", 0 * x), ("barrier 2000*1{|x|<0.05}", 2000.0 * (np.abs(x) < 0.05)),
                    ("harmonic 50 x^2", 50 * x**2)]:
        ev = eigh_tridiagonal(d0 + V, off, select="i", select_range=(0, 1), eigvals_only=True)
        say(f" {name:28s} gap = {ev[1] - ev[0]:.5f}   (PW pi^2/D^2 = {np.pi**2 / 4:.5f})")
    say("")


def c7_constant_field_hessian():
    say("== C7-b: Vmag non-convex at arbitrarily small |A| inside Omega, any L (constant fields) ==")
    for L in [3, 4]:
        lat = G.Lattice(L, 3)
        A = np.zeros((3, lat.V, 3)); A[0, :, 0] = 1.0; A[1, :, 1] = 1.0
        c = lat.Tsc.T @ A.reshape(3 * lat.V, 3); nrm = np.linalg.norm(c); c /= nrm
        assert np.allclose((lat.Tsc @ c).reshape(3, lat.V, 3) * nrm, A)  # constant field is transverse
        tp, tm, *_ = G.width(lat, c)
        for frac in [0.05, 0.2]:
            s = frac * min(tp, tm)
            x0 = (s * c).ravel()
            if L == 3:
                hm = G.hess_min_eig(lat, x0, h=1e-3)
            else:
                hm = float("nan")
            xsite = s / np.sqrt(2 * lat.V)  # per-link component
            say(f" L={L}: chord through 0 along (A1=e1, A2=e2 const) is [-{tm:.3f}, {tp:.3f}]; point |A|={s:.3f} "
                f"(inside Omega); analytic restricted-Hessian min = -x^2 = {-xsite**2:.3e}; full numeric Hessian min = {hm:.3e}")
    say(" Analytic: on constant fields Vmag = (V/2)|A1 x A2|^2; with A1 = x e1, A2 = y e2 in unit-norm coords"
        " the 2x2 Hessian at x=y has eigenvalues 3x^2, -x^2. A negative eigenvalue of a principal restriction"
        " implies one of the full Hessian.")
    say("")


def c7_diameters(rng):
    say("== C7-c: d=3 diameter lower bounds, more restarts, L = 3..6 ==")
    res = {}
    for L, rs in [(3, 12), (4, 8), (5, 4), (6, 3)]:
        t0 = time.time()
        lat = G.Lattice(L, 3)
        D = G.max_width(lat, rng, restarts=rs, iters=300)
        res[L] = D
        say(f" L={L} V={lat.V} restarts={rs} D_lower={D:.3f} D^2/V={D**2 / lat.V:.2f} D/L={D / L:.3f} [{time.time() - t0:.0f}s]")
    Ls = np.array(sorted(res)); Ds = np.array([res[L] for L in Ls])
    p = np.polyfit(np.log(Ls), np.log(Ds), 1)[0]
    say(f" fit D_lower ~ L^{p:.3f}  (D^2 ~ L^{2 * p:.2f}; 008 claimed D^2 between L^2 and L^3)")
    for a, b in zip(Ls[:-1], Ls[1:]):
        say(f"  local exponent L={a}->{b}: {np.log(res[b] / res[a]) / np.log(b / a):.3f}")
    say("")
    return res, p


def c7_ratio(res, p):
    m_s = None
    with open(ROOT / "data" / "glueball_sun_4d_athenodorou_teper_2021.csv") as f:
        for r in csv.DictReader(f):
            if (r["N"], r["JPC"], r["level"], r["quantity"]) == ("2", "0++", "0", "m/sqrt_sigma"):
                m_s = float(r["value"])
    say(f"== C7-d: kinetic PW quantity / m0++ (SU(2) 4D m/sqrt(sigma) = {m_s} from data/, AT2021), L sqrt(sigma) = 3 ==")
    say(" Measured lattices only (Omega diameters, noncompact discretisation):")
    for L, D in sorted(res.items()):
        for g2 in [1.0, 2.0, 4.0]:
            r = 0.5 * g2 * np.pi**2 / D**2 / (m_s * 3.0 / L)
            say(f"  N_s={L} g^2={g2}: <= {r:.3e}")
    say(f" Extrapolation with fitted D ~ L^{p:.2f} anchored at the largest L (not a bound):")
    Lmax = max(res); Dmax = res[Lmax]
    for Ns in [8, 16, 32]:
        D = Dmax * (Ns / Lmax) ** p
        say(f"  N_s={Ns} g^2=2: ~ {0.5 * 2 * np.pi**2 / D**2 / (m_s * 3.0 / Ns):.3e}")
    say(f" Continuum T^3, Lambda (constant-abelian chord, 008 (i)): Delta_PW <= g^2/(24 L) => ratio <= g^2/(24*3*{m_s}) "
        f"= {1 / (24 * 3 * m_s):.3e} g^2")
    say("")


def main():
    t0 = time.time()
    c4_checks()
    c7_pw_nonmonotone()
    c7_constant_field_hessian()
    res, p = c7_diameters(np.random.default_rng(2026))
    c7_ratio(res, p)
    say(f"total time {time.time() - t0:.0f}s")
    (HERE / "skeptic_summary.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
