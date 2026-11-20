"""Exact diagonalization of the SU(2) Kogut-Susskind Hamiltonian on tiny 2+1D lattices.

Writes results.csv (lattice, sector, g2, j_max, dim, E0, E1, Delta), checks.txt and gap_vs_g2.png
in this folder. Usage: python work/20261120-003-su2-ks-ed/run_ed.py [--quick]
"""
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402

OUT = Path(__file__).resolve().parent
G2 = [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]


def log(msg):
    print(msg, flush=True)


# ------------------------------------------------------------ single plaquette, group manifold
def mathieu_levels(g2, n=6000, k=2):
    """Class functions psi(theta), U ~ exp(i theta n.sigma), theta in (0, pi); phi = sin(theta) psi.
    4 links each contribute the full Laplacian on U_p: (g^2/2)*4*C, C -> (-phi'' - phi)/4.
    H phi = (g^2/2)(-phi'' - phi) - (4/g^2) cos(theta) phi, phi(0) = phi(pi) = 0.
    Second-order finite differences, Richardson-extrapolated in the grid spacing."""
    def solve(n):
        h = np.pi / (n + 1)
        th = h * np.arange(1, n + 1)
        d = 0.5 * g2 * (2 / h**2 - 1) - (4 / g2) * np.cos(th)
        e = -0.5 * g2 / h**2 * np.ones(n - 1)
        return eigh_tridiagonal(d, e, select='i', select_range=(0, k - 1))[0]
    a, b = solve(n), solve(2 * n + 1)  # h -> h/2
    return (4 * b - a) / 3


def sc_pt_plaq1(g2):
    """Strong-coupling PT to O(1/g^6): Delta = 3 g^2/2 + 56/(15 g^6)."""
    return 1.5 * g2 + 56 / (15 * g2**3)


# ------------------------------------------------------------ reflection parity for the ladder
def ladder_parity(lat, basis, W):
    """Reflection x -> 2-x. All vertices are <= 3-valent, so intertwiners are unique and the
    reflection maps each basis state to +- another one. Signs are fixed by requiring
    P W P = W along a spanning tree of W's graph from the vacuum (+1), then checked on all
    entries. Returns sparse P."""
    lperm = {0: 1, 1: 0, 2: 3, 3: 2, 4: 6, 5: 5, 6: 4}
    idx = {st: i for i, st in enumerate(basis)}
    pi = np.empty(len(basis), int)
    for i, (tjs, ints) in enumerate(basis):
        nt = [0] * len(tjs)
        for l, lp in lperm.items():
            nt[lp] = tjs[l]
        pi[i] = idx[(tuple(nt), tuple([0] * lat.nv))]
    Wc = W.tocsr()
    s = np.zeros(len(basis))
    s[0] = 1
    from collections import deque
    q = deque([0])
    while q:
        i = q.popleft()
        row = Wc.getrow(i)
        for k, val in zip(row.indices, row.data):
            if s[k] == 0:
                s[k] = np.real(val / (Wc[pi[i], pi[k]] * s[i]))
                q.append(k)
    err_sign = np.max(np.abs(np.abs(s) - 1))
    s = np.sign(s)  # remove rounding so sector tests are exact
    P = sp.csr_matrix((s, (pi, np.arange(len(basis)))), shape=W.shape)
    err_comm = abs(P @ W - W @ P).max()
    err_sq = abs(P @ P - sp.identity(len(basis))).max()
    return P, max(err_comm, err_sq, err_sign)


def sector_basis(P, sign):
    """Sparse orthonormal basis of the P = sign subspace (P is a signed permutation)."""
    P = P.tocsc()
    n = P.shape[0]
    cols, rows, vals = [], [], []
    seen = np.zeros(n, bool)
    c = 0
    for i in range(n):
        if seen[i]:
            continue
        j = P.indices[P.indptr[i]]  # P e_i = s_i e_j
        si = P.data[P.indptr[i]]
        seen[i] = seen[j] = True
        if j == i:
            if si == sign:
                rows.append(i); cols.append(c); vals.append(1.0); c += 1
        else:
            sj = P.data[P.indptr[j]]
            rows += [i, j]; cols += [c, c]; vals += [1 / np.sqrt(2), sign * sj / np.sqrt(2)]
            c += 1
    assert len(set(rows)) == len(rows)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, c))


# ------------------------------------------------------------ main
def main(quick=False):
    global OUT
    if quick:  # smoke test: write to gitignored artifacts/, never over the committed results
        OUT = ROOT / "artifacts" / "20261120-003-su2-ks-ed-quick"
        OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []
    checks = []

    # ---- single plaquette
    lat = ks.single_plaquette()
    jlist = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80]
    basis, cas, W = ks.build(lat, max(jlist))
    order = np.argsort([st[0][0] for st in basis])
    cas, W = cas[order], W[order][:, order]
    for g2 in G2:
        for tjm in jlist:
            sl = slice(0, tjm + 1)
            H = ks.hamiltonian(cas[sl], W[sl, sl], g2)
            e = ks.lowest(H, 2)
            rows.append(("plaq1", "all", g2, tjm / 2, tjm + 1, e[0], e[1], e[1] - e[0]))
        m = mathieu_levels(g2)
        dev = max(abs(e[0] - m[0]), abs(e[1] - m[1]))
        checks.append(f"plaq1 g2={g2}: spin basis j_max=40 E0={e[0]:.10f} E1={e[1]:.10f}; "
                      f"group-manifold FD E0={m[0]:.10f} E1={m[1]:.10f}; max dev {dev:.1e}; "
                      f"SC-PT Delta={sc_pt_plaq1(g2):.6f} vs ED {e[1]-e[0]:.6f}")
    log(f"plaq1 done {time.time()-t0:.0f}s")

    # ---- two-plaquette ladder, parity resolved
    lat = ks.ladder2()
    tjl = [1, 2, 3, 4, 6, 8, 10, 12, 14, 16] if not quick else [1, 2, 4]
    for tjm in tjl:
        tb = time.time()
        basis, cas, W = ks.build(lat, tjm)
        P, perr = ladder_parity(lat, basis, W)
        checks.append(f"ladder2 j_max={tjm/2}: dim {len(basis)}, |W-W^T|={abs(W-W.T).max():.1e}, "
                      f"reflection consistency err {perr:.1e}")
        Bp, Bm = sector_basis(P, +1), sector_basis(P, -1)
        Hd = ks.hamiltonian(cas, W, 1.0).toarray().real
        un = np.sort(np.r_[np.linalg.eigvalsh(Bp.T @ Hd @ Bp), np.linalg.eigvalsh(Bm.T @ Hd @ Bm)])
        checks.append(f"  ladder2 j_max={tjm/2}: dims P+ {Bp.shape[1]}, P- {Bm.shape[1]}; "
                      f"spec(P+) u spec(P-) vs full spectrum at g2=1: max dev "
                      f"{np.max(np.abs(un - np.linalg.eigvalsh(Hd))):.1e}")
        for g2 in G2:
            H = ks.hamiltonian(cas, W, g2).tocsr()
            e = ks.lowest(H, 2)
            ev = ks.lowest((Bp.T @ H @ Bp).tocsr(), 2)
            od = ks.lowest((Bm.T @ H @ Bm).tocsr(), 2)
            rows.append(("ladder2", "all", g2, tjm / 2, len(basis), e[0], e[1], e[1] - e[0]))
            rows.append(("ladder2", "P+", g2, tjm / 2, len(basis), ev[0], ev[1], ev[1] - ev[0]))
            rows.append(("ladder2", "P-", g2, tjm / 2, len(basis), ev[0], od[0], od[0] - ev[0]))
        log(f"ladder2 j_max={tjm/2} dim={len(basis)} {time.time()-tb:.0f}s")

    # ---- periodic 2x2 (zero 't Hooft electric flux sector: BFS from the vacuum), all momenta
    lat = ks.periodic(2, 2)
    tjl = [1, 2, 3] if not quick else [1, 2]
    for tjm in tjl:
        tb = time.time()
        basis, cas, W = ks.build(lat, tjm)
        checks.append(f"periodic2x2 j_max={tjm/2}: dim {len(basis)}, |W-W^dag|={abs(W-W.conj().T).max():.1e}, "
                      f"max |Im W| {np.max(np.abs(W.data.imag)) if W.nnz else 0:.1e}")
        for g2 in G2:
            H = ks.hamiltonian(cas, W, g2)
            H = 0.5 * (H + H.conj().T)
            e = ks.lowest(H.real if np.max(np.abs(H.data.imag)) < 1e-12 else H, 5)
            rows.append(("periodic2x2", "all", g2, tjm / 2, len(basis), e[0], e[1], e[1] - e[0]))
            checks.append(f"  periodic2x2 j_max={tjm/2} g2={g2}: lowest 5 levels {np.round(e, 6).tolist()}")
        log(f"periodic2x2 j_max={tjm/2} dim={len(basis)} {time.time()-tb:.0f}s")

    with open(OUT / "results.csv", "w") as f:
        f.write("lattice,sector,g2,j_max,dim,E0,E1,Delta\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]:.12g},{r[6]:.12g},{r[7]:.12g}\n")
    (OUT / "checks.txt").write_text("\n".join(checks) + "\n")
    log("\n".join(checks))
    summary(rows)
    plot(rows)
    log(f"total {time.time()-t0:.0f}s")


def summary(rows):
    """Converged gap per (lattice, sector, g2): value at largest j_max and change from the previous."""
    from collections import defaultdict
    d = defaultdict(list)
    for r in rows:
        d[(r[0], r[1], r[2])].append((r[3], r[7]))
    log("\nlattice  sector  g2  j_max  Delta  |Delta(j_max)-Delta(prev)|")
    for key in sorted(d, key=lambda k: (k[0], k[1], k[2])):
        v = sorted(d[key])
        dd = abs(v[-1][1] - v[-2][1]) if len(v) > 1 else float('nan')
        log(f"{key[0]:12s} {key[1]:4s} {key[2]:5.2f} {v[-1][0]:5.1f} {v[-1][1]:.8f} {dd:.1e}")


def plot(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    sets = [("plaq1", "all", "1 plaquette"), ("ladder2", "P+", "ladder, P=+"),
            ("ladder2", "P-", "ladder, P=-"), ("periodic2x2", "all", "2x2 periodic")]
    for lat, sec, lab in sets:
        sel = [r for r in rows if r[0] == lat and r[1] == sec]
        if not sel:
            continue
        jm = max(r[3] for r in sel)
        jprev = max(r[3] for r in sel if r[3] < jm)
        prev = {r[2]: r[7] for r in sel if r[3] == jprev}
        pts = sorted((r[2], r[7], abs(r[7] - prev[r[2]])) for r in sel if r[3] == jm)
        x, y, dd = map(np.array, zip(*pts))
        line, = ax.plot(x, y, "-", lw=1, label=f"{lab} (j_max={jm:g})")
        ok = dd < 1e-3
        ax.plot(x[ok], y[ok], "o", color=line.get_color())
        ax.plot(x[~ok], y[~ok], "o", mfc="none", color=line.get_color())
    g = np.geomspace(0.2, 10, 100)
    ax.plot(g, 1.5 * g, "k--", lw=1, label="3g^2/2 (strong coupling)")
    for y, lab in [(4, "4"), (2 * np.sqrt(3), "2*sqrt(3)")]:
        ax.axhline(y, color="gray", lw=0.6, ls=":")
        ax.text(0.21, y * 1.03, lab, fontsize=7, color="gray")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("g^2"); ax.set_ylabel("Delta = E1 - E0 (lattice units)")
    ax.set_title("SU(2) Kogut-Susskind, 2+1D, finite volume; (open marker: not converged in j_max to 1e-3)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "gap_vs_g2.png", dpi=110)


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
