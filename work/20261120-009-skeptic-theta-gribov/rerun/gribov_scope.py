"""Task B: geometry of the Gribov region for a noncompact (linear-A) SU(2) lattice Coulomb gauge.

Setting: spatial lattice Z_L^d (periodic), d = 2 or 3, lattice units a = 1, real fields A_mu^a(x),
a = 1..3. Coulomb (spatial Landau) condition: backward divergence sum_mu A_mu(x) - A_mu(x - mu) = 0.
FP quadratic form on adjoint ghosts omega (constants projected out):
  q_A(omega) = sum_{x,mu} |omega(x+mu) - omega(x)|^2 + sum_{x,mu} A_mu(x) . (omega(x) x omega(x+mu))
(the linear-A limit of the compact SU(2) lattice Landau FP form, normalization of the A term is a
convention). Gribov region Omega = {A transverse : q_A >= 0 on constants-orthogonal omega}.
q_A is affine in A, so Omega is an intersection of half-spaces, hence convex.

Outputs gribov_summary.txt and gribov_widths.csv. Usage:
  python work/20261120-008-theta-obstruction-gribov-scope/gribov_scope.py
"""
import csv
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, null_space

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
EPS = np.zeros((3, 3, 3))
for _a, _b, _c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    EPS[_a, _b, _c], EPS[_a, _c, _b] = 1, -1


class Lattice:
    def __init__(self, L, d):
        self.L, self.d, self.V = L, d, L**d
        idx = np.arange(self.V).reshape((L,) * d)
        self.fwd = [np.roll(idx, -1, axis=mu).ravel() for mu in range(d)]  # site x+mu
        div = np.zeros((self.V, self.V * d))
        for mu in range(d):
            for x in range(self.V):
                div[x, mu * self.V + x] += 1
                div[self.fwd[mu][x], mu * self.V + x] -= 1
        self.Tsc = null_space(div)  # orthonormal basis of transverse scalar link fields
        self.nT = self.Tsc.shape[1]
        lap = np.zeros((self.V, self.V))
        for mu in range(d):
            for x in range(self.V):
                y = self.fwd[mu][x]
                lap[x, x] += 1; lap[y, y] += 1; lap[x, y] -= 1; lap[y, x] -= 1
        w, U = eigh(lap)
        # drop the constant mode; P maps reduced ghost coords to 3V with M0^{-1/2} folded in
        self.P = np.kron(U[:, 1:] / np.sqrt(w[1:]), np.eye(3))

    def A_from_coeffs(self, c):
        return (self.Tsc @ c.reshape(self.nT, 3)).reshape(self.d, self.V, 3)

    def M1(self, A):
        V = self.V
        M = np.zeros((V, 3, V, 3))
        for mu in range(self.d):
            y = self.fwd[mu]
            B = np.einsum("xa,abc->xbc", A[mu], EPS)
            for x in range(V):
                M[x, :, y[x], :] += 0.5 * B[x]
                M[y[x], :, x, :] += 0.5 * B[x].T
        return M.reshape(3 * V, 3 * V)

    def K(self, A):
        return self.P.T @ self.M1(A) @ self.P

    def grad_coeffs(self, u):
        """Gradient wrt transverse coefficients of u^T M1(A) u."""
        u = u.reshape(self.V, 3)
        G = np.zeros((self.d, self.V, 3))
        for mu in range(self.d):
            G[mu] = np.cross(u, u[self.fwd[mu]])
        return self.Tsc.T @ G.reshape(self.d * self.V, 3)


def width(lat, c):
    """Omega contains t*Ahat for -t_- <= t <= t_+  (|c| = 1). M0 + t M1 >= 0 <=> 1 + t K >= 0."""
    w, v = np.linalg.eigh(lat.K(lat.A_from_coeffs(c)))
    tp = 1 / -w[0] if w[0] < 0 else np.inf
    tm = 1 / w[-1] if w[-1] > 0 else np.inf
    return tp, tm, w, v


def max_width(lat, rng, restarts, iters=150):
    """Gradient ascent of t_+ + t_- over unit directions. Any value found is a lower bound on diam."""
    best = 0.0
    for _ in range(restarts):
        c = rng.standard_normal((lat.nT, 3)); c /= np.linalg.norm(c)
        tp, tm, w, v = width(lat, c)
        cur, step = tp + tm, 0.3
        for _ in range(iters):
            g = (lat.grad_coeffs(lat.P @ v[:, 0]) / w[0] ** 2
                 - lat.grad_coeffs(lat.P @ v[:, -1]) / w[-1] ** 2)
            g -= np.sum(g * c) * c
            cn = c + step * g / (np.linalg.norm(g) + 1e-30)
            cn /= np.linalg.norm(cn)
            tpn, tmn, wn, vn = width(lat, cn)
            if tpn + tmn > cur:
                c, w, v, cur = cn, wn, vn, tpn + tmn
                step *= 1.2
            else:
                step *= 0.5
                if step < 1e-5:
                    break
        best = max(best, cur)
    return best


def Vmag(lat, A):
    """Noncompact lattice magnetic energy 1/2 sum |F|^2, F = dA + A_mu x A_nu."""
    tot = 0.0
    for mu in range(lat.d):
        for nu in range(mu + 1, lat.d):
            F = ((A[nu][lat.fwd[mu]] - A[nu]) - (A[mu][lat.fwd[nu]] - A[mu])
                 + np.cross(A[mu], A[nu]))
            tot += 0.5 * np.sum(F**2)
    return tot


def hess_min_eig(lat, x0, h=1e-3):
    f = lambda c: Vmag(lat, lat.A_from_coeffs(c))
    n = x0.size
    E = np.eye(n) * h
    f0 = f(x0)
    fp = np.array([f(x0 + E[i]) for i in range(n)])
    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = (fp[i] - 2 * f0 + f(x0 - E[i])) / h**2
        for j in range(i):
            H[i, j] = H[j, i] = (f(x0 + E[i] + E[j]) - fp[i] - fp[j] + f0) / h**2
    return np.linalg.eigvalsh(H)[0]


def main():
    rng = np.random.default_rng(12345)
    lines, rows = [], []
    t0 = time.time()
    lines.append("d L V dim(transverse) D_lower median_width min_width D_lower^2/V lowmode_width")
    for d, Ls in [(2, [3, 4, 5, 6, 8]), (3, [2, 3, 4])]:
        for L in Ls:
            lat = Lattice(L, d)
            D = max_width(lat, rng, restarts=3 if lat.V > 30 else 6)
            ws = []
            for _ in range(20):
                c = rng.standard_normal((lat.nT, 3)); c /= np.linalg.norm(c)
                tp, tm, *_ = width(lat, c)
                ws.append(tp + tm)
            # softest transverse direction: colour-3, A_x varying as cos(2 pi y / L)
            A = np.zeros((d, lat.V, 3))
            coords = np.array(np.unravel_index(np.arange(lat.V), (L,) * d)).T
            A[0, :, 2] = np.cos(2 * np.pi * coords[:, 1] / L)
            cl = lat.Tsc.T @ A.reshape(d * lat.V, 3); cl /= np.linalg.norm(cl)
            tp, tm, *_ = width(lat, cl)
            r = dict(d=d, L=L, V=lat.V, dim=3 * lat.nT, D_lower=D, median_width=np.median(ws),
                     min_width=np.min(ws), D2_over_V=D**2 / lat.V, lowmode_width=tp + tm)
            rows.append(r)
            lines.append(" ".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in r.values()))
            print(lines[-1], f"[{time.time() - t0:.0f}s]", flush=True)

    lat = Lattice(3, 3)
    worst = -np.inf
    for _ in range(10):
        c1 = rng.standard_normal((lat.nT, 3)); c1 /= np.linalg.norm(c1)
        c2 = rng.standard_normal((lat.nT, 3)); c2 /= np.linalg.norm(c2)
        p1, p2 = 0.9 * width(lat, c1)[0] * c1, 0.9 * width(lat, c2)[0] * c2
        ld = [np.linalg.slogdet(np.eye(3 * (lat.V - 1)) + lat.K(lat.A_from_coeffs((1 - t) * p1 + t * p2)))[1]
              for t in np.linspace(0, 1, 21)]
        worst = max(worst, np.max(np.diff(ld, 2)))
    lines.append(f"log det(M) along 10 random chords in Omega (d=3, L=3): max second difference = "
                 f"{worst:.3e} (<= 0: concave)")

    lat = Lattice(2, 3)
    for frac in [0.0, 0.25, 0.5, 0.75, 0.95]:
        c = rng.standard_normal((lat.nT, 3)); c /= np.linalg.norm(c)
        tp = width(lat, c)[0]
        m = hess_min_eig(lat, (frac * tp * c).ravel())
        lines.append(f"Vmag Hessian min eigenvalue (d=3, L=2) at A = {frac:.2f} t_+ Ahat "
                     f"(|A| = {frac * tp:.3f}): {m:.4f}")

    m_s = None
    with open(ROOT / "data" / "glueball_sun_4d_athenodorou_teper_2021.csv") as f:
        for r in csv.DictReader(f):
            if (r["N"], r["JPC"], r["level"], r["quantity"]) == ("2", "0++", "0", "m/sqrt_sigma"):
                m_s = float(r["value"])
    c3 = np.mean([r["D2_over_V"] for r in rows if r["d"] == 3][-2:])
    lines.append(f"SU(2) 4D m0++/sqrt(sigma) = {m_s} (data/, AT2021). d=3: D_lower^2 ~ {c3:.3f} V.")
    lines.append("Kinetic-only PW bound a*Delta_PW <= (g^2/2) pi^2 / D^2 with D^2 = c V = c N_s^3;")
    lines.append("a*m0++ = (m/sqrt sigma)(L sqrt sigma)/N_s.  Ratio Delta_PW/m0++ at L sqrt(sigma) = 3:")
    for g2 in [1.0, 2.0, 4.0]:
        for Ns in [4, 8, 16, 32]:
            ratio = 0.5 * g2 * np.pi**2 / (c3 * Ns**3) / (m_s * 3.0 / Ns)
            lines.append(f"  g^2={g2} N_s={Ns}: Delta_PW/m0++ <= {ratio:.3e}")
    with open(OUT / "gribov_widths.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader(); wr.writerows(rows)
    (OUT / "gribov_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
