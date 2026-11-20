"""Skeptic checks of C13 (work/20261120-020). Setting: 2+1D SU(2) Kogut-Susskind, a = 1, finite
lattice (single plaquette, open 2-plaquette ladder), fixed g^2; plus 1D spin chains.

 A. HF: build the plaquette-product state f(U_p1) f(U_p2) pointwise, fit it in the ladder spin-network
    basis, and evaluate <H> with the ladder Hamiltonian (independent of the 020 shortcut E_HF := 2 E0(plaq)).
 B. Selected CI without exact eigenvectors: two-root CIPSI (Epstein-Nesbet importance) from the vacuum;
    and state-averaged selection with exact vectors (|v0|^2 + |v1|^2).
 C. MPS ladder: (i) reproducibility of 020's exact-targeted truncation (zero-weight ties);
    (ii) self-consistent DMRG-style targeting (environments from Ritz vectors, iterated to a fixed point).
 D. Non-degenerate D = 1 gapless chain: open ferromagnetic Heisenberg chain with a boundary field.
 E. EOM-CC sign pattern from 020's cc_plaquette.csv.
Writes *.csv and summary.txt here.
"""
import csv
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402

OUT = Path(__file__).resolve().parent
C13 = ROOT / "work" / "20261120-020-quantum-chemistry-methods"
G2 = [0.5, 1.0, 2.0, 3.0]
LINES = []
RNG = np.random.default_rng(20261120)


def log(s=""):
    print(s, flush=True)
    LINES.append(s)


def write_csv(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([f"{x:.10g}" if isinstance(x, float) else x for x in r])


def setup(latf, tj):
    lat = latf()
    basis, cas, W = ks.build(lat, tj)
    return lat, basis, cas, W


# ------------------------------------------------------------------ A. HF pointwise check
def su2(v):
    return ks.D_matrix(1, v)


def rotvec(U):
    """Rotation vector v with U = exp(-i v.sigma/2) (fundamental)."""
    c = np.clip(np.real(np.trace(U)) / 2, -1, 1)
    th = 2 * np.arccos(c)
    s = np.sin(th / 2)
    if abs(s) < 1e-14:
        return np.zeros(3)
    sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1])
    n = np.array([np.real(1j * np.trace(U @ m)) / (2 * s) for m in (sx, sy, sz)])
    return th * n


def hol(lat, p, Us):
    M = np.eye(2, dtype=complex)
    for l, d in lat.plaqs[p]:
        u = su2(Us[l])
        M = M @ (u if d > 0 else u.conj().T)
    return M


def hf_check():
    log("== A. HF: plaquette-product state evaluated with the ladder Hamiltonian ==")
    T = 6  # plaquette truncation tj <= 6 (converged to < 1e-10 at g^2 >= 2, see 020 ci_truncation jmax)
    plat, pb, pcas, pW = setup(ks.single_plaquette, T)
    llat, lb, lcas, lW = setup(ks.ladder2, 2 * T)
    keep = [i for i, st in enumerate(lb) if st[0][0] <= T and st[0][1] <= T]
    rows = []
    nsamp = 3 * len(keep)
    Us = RNG.normal(size=(nsamp, 7, 3)) * 2.0
    # ladder spin networks at the samples (column = basis state)
    A = np.array([[ks.wavefunction(llat, lb[i], Us[s]) for i in keep] for s in range(nsamp)])
    hv = [[rotvec(hol(llat, p, Us[s])) for p in (0, 1)] for s in range(nsamp)]
    pvals = {}  # plaquette spin networks at (W,1,1,1)
    for s in range(nsamp):
        for p in (0, 1):
            Uq = [hv[s][p], np.zeros(3), np.zeros(3), np.zeros(3)]
            pvals[(s, p)] = np.array([ks.wavefunction(plat, st, Uq) for st in pb])
    for g2 in G2:
        Hp = ks.hamiltonian(pcas, pW, g2).toarray().real
        w, V = np.linalg.eigh(Hp)
        f = V[:, 0]
        y = np.array([(pvals[(s, 0)] @ f) * (pvals[(s, 1)] @ f) for s in range(nsamp)])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = np.linalg.norm(A @ c - y) / np.linalg.norm(y)
        v = np.zeros(len(lb), complex); v[keep] = c
        Hl = ks.hamiltonian(lcas, lW, g2)
        ehf = np.real(v.conj() @ (Hl @ v)) / np.real(v.conj() @ v)
        el = ks.lowest(Hl, 2)
        rows.append((g2, nsamp, resid, float(np.linalg.norm(c)), ehf, 2 * w[0], ehf - 2 * w[0], el[0]))
        log(f"g2={g2:<4} fit resid {resid:.1e}, |c| = {np.linalg.norm(c):.12f}; <H_ladder> = {ehf:.10f}, "
            f"2 E0(plaq) = {2 * w[0]:.10f}, diff {ehf - 2 * w[0]:+.1e}; E0(ladder, tj<=12) = {el[0]:.6f}")
    write_csv("hf_check.csv", ["g2", "samples", "fit_residual", "norm", "E_HF_ladder", "2E0_plaq",
                               "diff", "E0_ladder"], rows)
    log("   (g^2 = 0.5, 1 need tj > 6 on the plaquette; the identity is algebraic, checked where converged)")
    log()


# ------------------------------------------------------------------ B. selected CI
def ritz2(Hd, idx):
    w, U = np.linalg.eigh(Hd[np.ix_(idx, idx)])
    return w, U


def cipsi(Hd, sizes):
    """Two-root CIPSI from the vacuum: add configurations with largest max_r |<a|H|psi_r>|^2/|E_r - H_aa|."""
    n = Hd.shape[0]
    sel = [0] + list(np.nonzero(Hd[0])[0][np.nonzero(Hd[0])[0] != 0])
    out = []
    for target in sizes:
        while len(sel) < target:
            w, U = ritz2(Hd, sel)
            ext = np.setdiff1d(np.arange(n), sel)
            imp = np.zeros(len(ext))
            for r in range(min(2, len(w))):
                psi = np.zeros(n); psi[sel] = U[:, r]
                coup = Hd[ext] @ psi
                imp = np.maximum(imp, coup ** 2 / np.maximum(np.abs(w[r] - Hd[ext, ext]), 1e-12))
            add = ext[np.argsort(-imp, kind="stable")][: max(1, min(target - len(sel), len(sel)))]
            sel = sorted(set(sel) | set(add.tolist()))
        w, _ = ritz2(Hd, sel[:target] if len(sel) > target else sel)
        out.append((len(sel), w[0], w[1]))
    return out


def ci_check():
    log("== B. Selected CI on the ladder (tj <= 12): honest selections ==")
    lat, basis, cas, W = setup(ks.ladder2, 12)
    rows = []
    sizes = [4, 6, 8, 12, 16, 24, 32, 48, 64]
    for g2 in G2:
        Hd = ks.hamiltonian(cas, W, g2).toarray().real
        w, V = np.linalg.eigh(Hd)
        D = w[1] - w[0]
        res = cipsi(Hd, sizes)
        for (m, e0, e1) in res:
            rows.append((g2, "cipsi2", m, e0, e1, e1 - e0, D, (e1 - e0) - D))
        wt = V[:, 0] ** 2 + V[:, 1] ** 2
        order = np.argsort(-wt, kind="stable")
        for m in sizes:
            ww, _ = ritz2(Hd, np.sort(order[:m]))
            rows.append((g2, "sel_avg", m, ww[0], ww[1], ww[1] - ww[0], D, (ww[1] - ww[0]) - D))
        # 020's sel_E1, and the separate state-specific variant: E0 from sel_E0 space, E1 from sel_E1 space
        for m in sizes:
            o0 = np.sort(np.argsort(-np.abs(V[:, 0]), kind="stable")[:m])
            o1 = np.sort(np.argsort(-np.abs(V[:, 1]), kind="stable")[:m])
            w0, _ = ritz2(Hd, o0); w1, _ = ritz2(Hd, o1)
            rows.append((g2, "sel_E1", m, w1[0], w1[1], w1[1] - w1[0], D, (w1[1] - w1[0]) - D))
            rows.append((g2, "state_specific", m, w0[0], w1[1], w1[1] - w0[0], D, (w1[1] - w0[0]) - D))
        for sch in ("cipsi2", "sel_avg", "sel_E1", "state_specific"):
            r = [x for x in rows if x[0] == g2 and x[1] == sch]
            log(f"g2={g2:<4} {sch:14s} gap err by dim {','.join(str(x[2]) for x in r)}: "
                + " ".join(f"{x[7]:+.2e}" for x in r))
    write_csv("ci_selected.csv", ["g2", "scheme", "dim", "E0_V", "E1_V", "gap_V", "gap_exact", "gap_error"], rows)
    log()


# ------------------------------------------------------------------ C. MPS ladder
def mps_setup(basis):
    a = np.array([st[0][0] for st in basis]); b = np.array([st[0][1] for st in basis])
    c = np.array([st[0][5] for st in basis])
    return a, b, c


def candidates(vecs, wts, a, b, c, nA, tol):
    cand = []
    for cc in np.unique(c):
        sel = np.where(c == cc)[0]
        rhoL = np.zeros((nA, nA)); rhoR = np.zeros((nA, nA))
        for wgt, v in zip(wts, vecs):
            if wgt == 0:
                continue
            M = np.zeros((nA, nA)); M[a[sel], b[sel]] = v[sel]
            rhoL += wgt * M @ M.T; rhoR += wgt * M.T @ M
        for side, rho in (("L", rhoL), ("R", rhoR)):
            w, U = np.linalg.eigh(rho)
            for i in range(nA):
                if w[i] > tol:
                    cand.append((w[i], side, int(cc), U[:, i]))
    return cand


def product_ritz(Hd, cand, Dk, a, b, c):
    keepL = sorted([x for x in cand if x[1] == "L"], key=lambda x: (-x[0], x[2]))[:Dk]
    keepR = sorted([x for x in cand if x[1] == "R"], key=lambda x: (-x[0], x[2]))[:Dk]
    vecs = []
    for (_, _, cl, uL) in keepL:
        for (_, _, cr, uR) in keepR:
            if cl != cr:
                continue
            sel = np.where(c == cl)[0]
            v = np.zeros(len(a)); v[sel] = uL[a[sel]] * uR[b[sel]]
            if np.linalg.norm(v) > 1e-12:
                vecs.append(v)
    if len(vecs) < 2:
        return None
    Q, _ = np.linalg.qr(np.array(vecs).T)
    w, U = np.linalg.eigh(Q.T @ Hd @ Q)
    return w, Q @ U


def mps_check():
    log("== C. MPS ladder (cut through the shared link, D multiplets per side) ==")
    lat, basis, cas, W = setup(ks.ladder2, 12)
    a, b, c = mps_setup(basis)
    nA = a.max() + 1
    rows = []
    Ds = (2, 3, 4, 6, 8, 12, 16, 24)
    targets = (((1.0, 0.0), "E0"), ((0.5, 0.5), "avg"), ((0.0, 1.0), "E1"))
    for g2 in G2:
        Hd = ks.hamiltonian(cas, W, g2).toarray().real
        w, V = np.linalg.eigh(Hd)
        D = w[1] - w[0]
        for wts, tag in targets:
            # (i) 020 recipe with exact vectors, zero-weight Schmidt vectors excluded (deterministic)
            cand = candidates([V[:, 0], V[:, 1]], wts, a, b, c, nA, 1e-13)
            nz = sum(1 for x in cand if x[1] == "L")
            for Dk in Ds:
                r = product_ritz(Hd, cand, Dk, a, b, c)
                if r is not None:
                    rows.append((g2, "exact_" + tag, Dk, nz, r[0][0], r[0][1], D, (r[0][1] - r[0][0]) - D, 0))
            # (ii) self-consistent: start from lowest-Casimir block states, environments from Ritz vectors
            for Dk in Ds:
                init = []
                for cc in np.unique(c):
                    for side, lab in (("L", a), ("R", b)):
                        for av in np.unique(lab[c == cc]):
                            u = np.zeros(nA); u[av] = 1.0
                            init.append((1.0 / (1 + av * (av + 2) / 4 + cc * (cc + 2) / 4), side, int(cc), u))
                cand2 = init
                prev = None
                for it in range(200):
                    r = product_ritz(Hd, cand2, Dk, a, b, c)
                    if r is None:
                        break
                    ev, X = r
                    vecs = [X[:, 0], X[:, 1] if X.shape[1] > 1 else X[:, 0]]
                    cur = (ev[0], ev[1] if len(ev) > 1 else np.inf)
                    if prev is not None and abs(cur[0] - prev[0]) + abs(cur[1] - prev[1]) < 1e-12:
                        break
                    prev = cur
                    cand2 = candidates(vecs, wts, a, b, c, nA, 1e-13)
                    # pad with fresh low-Casimir states so the block never shrinks below D (DMRG noise analogue)
                    have = {(x[1], x[2], int(np.argmax(np.abs(x[3])))) for x in cand2}
                    cand2 += [(1e-14 * x[0], *x[1:]) for x in init
                              if (x[1], x[2], int(np.argmax(x[3]))) not in have]
                if r is not None:
                    rows.append((g2, "scf_" + tag, Dk, it + 1, cur[0], cur[1], D, (cur[1] - cur[0]) - D,
                                 cur[0] - w[0]))
        for tag in ("exact_E0", "exact_avg", "exact_E1", "scf_E0", "scf_avg", "scf_E1"):
            r = [x for x in rows if x[0] == g2 and x[1] == tag]
            log(f"g2={g2:<4} {tag:9s} gap err D={','.join(str(x[2]) for x in r)}: "
                + " ".join(f"{x[7]:+.2e}" for x in r))
    write_csv("mps_check.csv", ["g2", "scheme", "D", "n_nonzero_or_iters", "E0_D", "E1_D", "gap_exact",
                                "gap_error", "E0_err"], rows)
    log()


def mps_zero_weight_diag():
    """How many 020 'kept' Schmidt vectors had zero weight (ties broken by floating-point noise)?"""
    lat, basis, cas, W = setup(ks.ladder2, 12)
    a, b, c = mps_setup(basis)
    nA = a.max() + 1
    out = []
    for g2 in G2:
        Hd = ks.hamiltonian(cas, W, g2).toarray().real
        w, V = np.linalg.eigh(Hd)
        for wts, tag in (((1.0, 0.0), "E0"), ((0.5, 0.5), "avg"), ((0.0, 1.0), "E1")):
            cand = candidates([V[:, 0], V[:, 1]], wts, a, b, c, nA, -1.0)
            L = sorted([x[0] for x in cand if x[1] == "L"], reverse=True)
            nz = sum(1 for x in L if x > 1e-13)
            out.append((g2, tag, nz, len(L)))
    log("   020 recipe: number of left Schmidt vectors with weight > 1e-13 (of all candidates): "
        + "; ".join(f"g2={g} {t}: {nz}/{n}" for g, t, nz, n in out))
    log("   -> any kept D above that count is chosen among zero-weight vectors by rounding noise")
    log()


# ------------------------------------------------------------------ D. non-degenerate gapless D = 1 example
def kron_op(ops, L):
    out = sp.identity(1, format="csr")
    for i in range(L):
        out = sp.kron(out, ops.get(i, sp.identity(2, format="csr")), format="csr")
    return out


def chain_check():
    log("== D. H = sum_i (1/4 - S_i.S_{i+1}) + (1/2 - S^z_1), open chain: unique product ground state ==")
    Sp = sp.csr_matrix([[0, 1], [0, 0]], dtype=float); Sm = Sp.T.tocsr()
    Sz = sp.csr_matrix([[0.5, 0], [0, -0.5]])
    rows = []
    for L in (6, 8, 10, 12, 14):
        n = 2 ** L
        H = sp.csr_matrix((n, n))
        for i in range(L - 1):
            H = H + 0.25 * sp.identity(n) - kron_op({i: Sz, i + 1: Sz}, L) \
                - 0.5 * (kron_op({i: Sp, i + 1: Sm}, L) + kron_op({i: Sm, i + 1: Sp}, L))
        H = H + 0.5 * sp.identity(n) - kron_op({0: Sz}, L)
        w = np.sort(eigsh(H.tocsc(), k=3, sigma=-0.1, which="LM", return_eigenvectors=False))
        up = np.zeros(n); up[0] = 1.0
        e_up = up @ (H @ up)
        rows.append((L, e_up, w[0], w[1], w[1] - w[0], (w[1] - w[0]) * (2 * L + 1) ** 2))
        log(f"L={L:2d} <up..up|H|up..up> = {e_up:.2e}, E0 = {w[0]:.2e}, E1 = {w[1]:.6f}, E2 = {w[2]:.6f}, "
            f"gap*(2L+1)^2 = {(w[1] - w[0]) * (2 * L + 1) ** 2:.4f}")
    write_csv("chain_nondegenerate.csv", ["L", "E_product", "E0", "E1", "gap", "gap_times_(2L+1)^2"], rows)
    log("   one-magnon: Dirichlet at site 1, Neumann at site L: eps = 1 - cos(pi/(2L+1)) ~ pi^2/(2(2L+1)^2)")
    log(f"   pi^2/2 = {np.pi ** 2 / 2:.4f}")
    log()


# ------------------------------------------------------------------ E. EOM-CC sign pattern
def eom_signs():
    log("== E. 020 EOM-CC gap-error sign pattern (n = 1..10, |err| > 1e-12) ==")
    rows = list(csv.DictReader(open(C13 / "cc_plaquette.csv")))
    for g2 in ("0.5", "1", "2", "3"):
        errs = [float(r["gap_error"]) for r in rows if r["g2"] == g2]
        ecc = [float(r["E_CC_minus_E0"]) for r in rows if r["g2"] == g2]
        s = "".join("+" if e > 1e-12 else "-" if e < -1e-12 else "0" for e in errs)
        s2 = "".join("+" if e > 1e-12 else "-" if e < -1e-12 else "0" for e in ecc)
        log(f"g2={g2:<4} EOM gap err signs {s}   E_CC - E0 signs {s2}")
    log()


def main():
    hf_check()
    ci_check()
    mps_zero_weight_diag()
    mps_check()
    chain_check()
    eom_signs()
    (OUT / "summary.txt").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
