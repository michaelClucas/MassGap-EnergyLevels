"""Spin-1 chains: AKLT (Knabe / Gosset-Mozgunov finite-size criteria) and Heisenberg (not frustration-free).
Setting: 1D lattice, spin 1, finite chains (open n sites, rings N sites), exact diagonalization in the Sz = 0
sector (every SU(2) multiplet has an Sz = 0 member, so gaps are full-space gaps).
Bond term family h(theta) = cos(theta) S.S + sin(theta) (S.S)^2; AKLT at tan(theta) = 1/3; Heisenberg at 0.
AKLT projector form: P2(i,i+1) = (1/2) S.S + (1/6)(S.S)^2 + 1/3.
Run: python tools/pin.py .venv/Scripts/python.exe -u work/20261120-015-line6-rotor-chain/spin1.py
Writes spin1_aklt_knabe.csv, spin1_path.csv."""
import csv, os, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
s2 = np.sqrt(2)
Sz = np.diag([1.0, 0.0, -1.0])
Sp = np.array([[0, s2, 0], [0, 0, s2], [0, 0, 0]])
Sm = Sp.T
SS = np.kron(Sz, Sz) + 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp))
SS2 = SS @ SS
P2 = 0.5 * SS + SS2 / 6 + np.eye(9) / 3
assert np.allclose(P2 @ P2, P2) and abs(np.trace(P2) - 5) < 1e-12


def chain(op2, L, ring):
    idx_cache = None
    mz = np.array([1, 0, -1])
    tot = np.zeros(1, dtype=int)
    for _ in range(L):
        tot = (tot[:, None] + mz[None, :]).ravel()
    idx = np.nonzero(tot == 0)[0]
    H = None
    bl = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if ring else [])
    for i, j in bl:
        if j == i + 1:
            t = sp.kron(sp.kron(sp.identity(3 ** i), sp.csr_matrix(op2)), sp.identity(3 ** (L - j - 1)), format="csr")
        else:
            U, s, Vt = np.linalg.svd(op2.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9))
            t = None
            for k in range(9):
                if s[k] < 1e-14:
                    continue
                Xa = (U[:, k] * s[k]).reshape(3, 3); Xb = Vt[k].reshape(3, 3)
                tt = sp.kron(sp.kron(sp.csr_matrix(Xb), sp.identity(3 ** (L - 2))), sp.csr_matrix(Xa), format="csr")
                t = tt if t is None else t + tt
        t = t[idx][:, idx]
        H = t if H is None else H + t
    return H.tocsr()


def low(H, k):
    if H.shape[0] <= 2000:
        return np.linalg.eigvalsh(H.toarray())[:k]
    return np.sort(spla.eigsh(H, k=k, which="SA", tol=1e-12)[0])


def main():
    t0 = time.time()
    # ---- AKLT: open-chain gaps eps_n (first nonzero eigenvalue), ring gaps, Knabe and GM bounds
    rows = []
    for n in range(3, 13):
        H = chain(P2, n, False)
        w = low(H, 8)
        eps = w[w > 1e-8][0]
        knabe = (n - 1) / (n - 2) * (eps - 1 / (n - 1))
        gm = 5 / 6 * (n * n + n) / (n * n - 4) * (eps - 6 / (n * (n + 1)))
        wr = low(chain(P2, n, True), 3)
        rows.append([n, f"{w[0]:.2e}", f"{eps:.8f}", f"{knabe:.6f}", f"{gm:.6f}", f"{wr[1] - wr[0]:.8f}"])
        print("AKLT", *rows[-1], f"{time.time() - t0:.0f}s", flush=True)
    with open(os.path.join(HERE, "spin1_aklt_knabe.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "open_E0", "open_gap_eps_n", "knabe_ring_lower", "gm_ring_lower(m>2n)", "ring_gap_n"])
        w.writerows(rows)
    # ---- path theta from Heisenberg (0) to AKLT (atan 1/3): frustration and ring gap
    th_aklt = np.arctan(1 / 3)
    rows = []
    for th in np.linspace(0, th_aklt, 6):
        hb = np.cos(th) * SS + np.sin(th) * SS2
        emin = np.linalg.eigvalsh(hb)[0]
        for N in (8, 10, 12):
            H = chain(hb, N, True)
            w = low(H, 3)
            frust = (w[0] - N * emin) / N  # >= 0; = 0 iff frustration-free
            rows.append([f"{th:.6f}", N, f"{emin:.6f}", f"{w[0] / N:.8f}", f"{frust:.6f}", f"{w[1] - w[0]:.8f}"])
            print("path", *rows[-1], flush=True)
    # distance of normalized bond terms from AKLT: ||h(theta)/c - P2|| with the best scale and shift
    with open(os.path.join(HERE, "spin1_path.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta", "N_ring", "bond_min_eig", "E0_per_bond", "frustration_per_bond", "ring_gap"])
        w.writerows(rows)
    # Heisenberg open-chain eps_n for the record (Knabe needs FF; shown to compare with 1/(n-1))
    for n in (6, 8, 10, 12):
        H = chain(SS - np.linalg.eigvalsh(SS)[0] * np.eye(9), n, False)
        w = low(H, 6)
        print("Heis open n", n, "E0 above sum of bond minima:", f"{w[0]:.6f}", "lowest gaps", np.round(w[1:4] - w[0], 6), flush=True)


if __name__ == "__main__":
    main()
