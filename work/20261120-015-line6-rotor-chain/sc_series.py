"""Strong-coupling (x = 1/g^2) Rayleigh-Schroedinger series for the O(3) rotor chain gap.
H/g = (1/2) sum L^2 + x V,  V = -sum n_i.n_{i+1}, periodic ring of L sites, momentum k = 0.
Ground state: sector (M = 0, k = 0), unperturbed |all l = 0>. Gap state: sector (M = 1, k = 0),
unperturbed (1/sqrt L) sum_i |l_i = 1, m_i = 1>; both nondegenerate in their sectors.
Coefficient of x^n equals the infinite-chain one for n <= L - 1 (a term winding the ring needs >= L bond
operators) and n <= 2 lmax - 1 (the gap state already has l = 1 at one site; reaching l = lmax + 1 there and
returning needs >= 2 lmax operators; checked: lmax = 1 and 2 differ at order 2).
Run: python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/sc_series.py
Writes sc_series.csv."""
import csv, os, sys
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rotor

HERE = os.path.dirname(os.path.abspath(__file__))


def k0_basis(L, lmax, M):
    d = rotor.site_ops(lmax)["d"]
    idx = rotor.m_sector(L, lmax, M)
    top = d ** (L - 1)
    shifts = [idx]
    s = idx.copy()
    for _ in range(L - 1):
        s = (s % top) * d + s // top
        shifts.append(s)
    S = np.stack(shifts)
    rep = S.min(axis=0)
    nsz = np.array([len(np.unique(S[:, c])) for c in range(S.shape[1])])
    reps, col = np.unique(rep, return_inverse=True)
    vals = 1 / np.sqrt(nsz)
    return sp.csr_matrix((vals, (np.arange(len(idx)), col)), shape=(len(idx), len(reps)))


def series(L, lmax, M, nmax):
    T, V, _ = rotor.chain_parts(L, "per", lmax, M)
    P = k0_basis(L, lmax, M)
    H0 = np.asarray((P.T @ (0.5 * T) @ P).diagonal()).ravel()
    off = P.T @ (0.5 * T) @ P - sp.diags(H0)
    assert abs(off).max() < 1e-12
    Vk = (-(P.T @ V @ P)).tocsr()
    e0 = H0.min()
    cand = np.nonzero(abs(H0 - e0) < 1e-12)[0]
    assert len(cand) == 1, cand
    i0 = cand[0]
    phi = np.zeros(len(H0)); phi[i0] = 1
    den = e0 - H0; den[i0] = np.inf
    R = 1 / den
    psi = [phi]; E = [e0]
    for n in range(1, nmax + 1):
        En = (Vk @ psi[n - 1])[i0]
        E.append(En)
        rhs = Vk @ psi[n - 1] - sum(E[k] * psi[n - k] for k in range(1, n + 1))
        p = R * rhs; p[i0] = 0
        psi.append(p)
    return np.array(E)


def main():
    rows = []
    nmax = 12
    for L, lmax in [(4, 2), (5, 2), (6, 2), (6, 1), (5, 3)]:
        E0 = series(L, lmax, 0, nmax)
        E1 = series(L, lmax, 1, nmax)
        gap = E1 - E0
        exact = min(L - 1, 2 * lmax - 1)  # gap state starts at l = 1, so l = lmax + 1 is reached at order lmax
        for n in range(nmax + 1):
            rows.append([L, lmax, n, f"{E0[n]:.12g}", f"{gap[n]:.12g}", int(n <= exact)])
        print(L, lmax, "gap coeffs", np.round(gap[:8], 8), "exact through", exact, flush=True)
        print(L, lmax, "E0 coeffs ", np.round(E0[:8], 8), flush=True)
    with open(os.path.join(HERE, "sc_series.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L_ring", "lmax", "order_n", "E0_coeff_per_g", "gap_coeff_per_g", "equals_infinite_chain"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
