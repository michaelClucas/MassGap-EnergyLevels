"""Task A: Bloch bands of the quantum pendulum with a theta term, and what they control.

H = p^2/(2I) + V0 (1 - cos x) on the line; T: x -> x + 2 pi commutes with H ("large gauge
transformation"). Fiber H_theta acts on psi(x + 2 pi) = e^{i theta} psi(x), i.e. Fourier modes
e^{i (m + theta/2pi) x}. Units I = 1. Writes theta_bands.csv and theta_summary.txt.
Usage: python work/20261120-008-theta-obstruction-gribov-scope/theta_bands.py
"""
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

OUT = Path(__file__).resolve().parent
MMAX = 60


def bands(V0, theta, k=3):
    m = np.arange(-MMAX, MMAX + 1)
    d = 0.5 * (m + theta / (2 * np.pi)) ** 2 + V0
    e = np.full(len(m) - 1, -V0 / 2)
    return eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1), eigvals_only=True)


def chi_exact(V0):
    """E0''(0) by second-order perturbation theory in theta (exact formula, no finite differences)."""
    m = np.arange(-MMAX, MMAX + 1)
    H = np.diag(0.5 * m**2 + V0) + np.diag(np.full(len(m) - 1, -V0 / 2), 1)         + np.diag(np.full(len(m) - 1, -V0 / 2), -1)
    w, v = np.linalg.eigh(H)
    dH = m / (2 * np.pi)  # dH/dtheta, diagonal
    me = v.T @ (dH * v[:, 0])
    return 1 / (4 * np.pi**2) + 2 * np.sum(me[1:] ** 2 / (w[0] - w[1:]))


def main():
    thetas = np.linspace(0, 2 * np.pi, 81)
    rows, summ = [], []
    summ.append("V0  s=sqrt(V0)  omega  W0=E0(pi)-E0(0)  chi=E0''(0)  W1  "
                "direct_gap(0)  direct_gap(pi)  indirect_gap  log(W0)+8s-1.5log(s)  N_overlap=ceil(gap/W0)")
    for V0 in [0.05, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
        E = np.array([bands(V0, t) for t in thetas])
        for t, ev in zip(thetas, E):
            rows.append((V0, t, *ev))
        chi = chi_exact(V0)
        W0 = E[:, 0].max() - E[:, 0].min()
        W1 = E[:, 1].max() - E[:, 1].min()
        dg = E[:, 1] - E[:, 0]
        ind = E[:, 1].min() - E[:, 0].max()
        s = np.sqrt(V0)
        summ.append(f"{V0:5.2f} {s:6.3f} {s:7.4f} {W0:.4e} {chi:.4e} {W1:.4e} "
                    f"{dg[0]:.5f} {dg[40]:.5f} {ind:.5f} {np.log(W0) + 8 * s - 1.5 * np.log(s):7.3f} "
                    f"{int(np.ceil(dg.min() / W0))}")
    # weak-potential Bragg check: gap at theta = pi between bands 0 and 1 ~ 2|V_1| = V0
    summ.append("")
    summ.append("Bragg (nearly-free) check at theta=pi: V0, E1-E0, ratio to V0")
    for V0 in [0.001, 0.01, 0.05]:
        ev = bands(V0, np.pi)
        summ.append(f"{V0:.3f} {ev[1] - ev[0]:.6e} {(ev[1] - ev[0]) / V0:.5f}")
    # Extensivity: N decoupled pendula, theta couples to total winding. Fiber gap = single gap,
    # n=0 band width = N*W0. Indirect gap < 0 once N*W0 > direct gap.
    summ.append("")
    summ.append("Extensivity (N independent wells, common theta): V0=4, N, n=0 bandwidth N*W0, fiber gap")
    E = np.array([bands(4.0, t) for t in thetas])
    W0 = E[:, 0].max() - E[:, 0].min()
    gmin = (E[:, 1] - E[:, 0]).min()
    for N in [1, 10**3, 10**5, 10**6, 10**7]:
        summ.append(f"N={N:6d} bandwidth={N * W0:.4e} fiber_gap_min={gmin:.5f} "
                    f"indirect={'>0' if N * W0 < gmin else '<0 (bands overlap)'}")
    np.savetxt(OUT / "theta_bands.csv", np.array(rows), delimiter=",",
               header="V0,theta,E0,E1,E2", comments="", fmt="%.10g")
    (OUT / "theta_summary.txt").write_text("\n".join(summ) + "\n")
    print("\n".join(summ))


if __name__ == "__main__":
    main()
