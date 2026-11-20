"""Quantum-chemistry methods as gap estimators on tiny SU(2) Kogut-Susskind lattices.

Setting of every number: 2+1D Hamiltonian lattice, SU(2), a = 1, fixed g^2, finite lattice
(single plaquette, open 2-plaquette ladder), H as in src/ks_su2.py. Plus two 1D spin-chain
checks for the MPS/area-law direction. Nothing here is about the continuum.

Methods (all compared with the converged ED gap of the same finite lattice, which C5 certifies):
  CI   : Ritz values in the span of states <= k plaquette moves from the strong-coupling vacuum
         (the analogue of CIS, CISD, ... built on the HF/strong-coupling reference).
  JMAX : Ritz values with link spins <= j_max (the C1 truncation, also a CI hierarchy).
  HF   : best gauge-invariant product state of the ladder = product of plaquette class functions;
         energy 2 E0(plaq), Koopmans-type gap Delta(plaq).
  CC   : coupled-cluster method on the single plaquette, exp(S) |0>, S = sum_{tj<=n} s_tj chi_tj
         (the Schuette-Wichmann et al. ansatz restricted to one plaquette), projected equations;
         EOM-CC gap from the similarity-transformed H restricted to characters tj <= n.
  MPS  : ladder cut through the shared link; per electric sector c keep the D_c leading Schmidt
         vectors (left, right) of a target-weighted density matrix built from the exact E0/E1 states;
         Ritz in the product span (= best two-site DMRG centre with those environments).
  TEMPLE: Temple lower bound for E1 from the CI-k E1 vector with beta = CI-k E2 (the only
         number available without an a priori bound) -> checks whether the "bound" is false.
  CHAIN: TT-SVD MPS truncations of the critical and gapped transverse-field Ising chain and the
         ferromagnetic Heisenberg chain (exact D = 1 ground state, gapless).
Writes *.csv and summary.txt in this folder.
"""
import sys
from collections import deque
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import fsolve

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402

OUT = Path(__file__).resolve().parent
G2 = [0.5, 1.0, 2.0, 3.0]
TJ_EXACT = {"plaq1": 40, "ladder2": 12}  # converged at these g^2 (C1/C5)
LINES = []


def log(s=""):
    print(s, flush=True)
    LINES.append(s)


def write_csv(name, header, rows):
    with open(OUT / name, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(f"{x:.10g}" if isinstance(x, float) else str(x) for x in r) + "\n")


def ritz(Hd, idx, k=3):
    w = np.linalg.eigvalsh(Hd[np.ix_(idx, idx)])
    w = np.concatenate([w, [np.inf] * k])
    return w[:k]


def bfs_depth(W):
    A = (abs(W) > 0).tocsr()
    d = np.full(A.shape[0], -1)
    d[0] = 0
    q = deque([0])
    while q:
        i = q.popleft()
        for j in A.indices[A.indptr[i]:A.indptr[i + 1]]:
            if d[j] < 0:
                d[j] = d[i] + 1
                q.append(j)
    return d


# ------------------------------------------------------------------ CI hierarchies
def ci_runs():
    rows = []
    exact = {}
    for name, latf in (("plaq1", ks.single_plaquette), ("ladder2", ks.ladder2)):
        lat = latf()
        basis, cas, W = ks.build(lat, TJ_EXACT[name])
        depth = bfs_depth(W)
        tjmax_st = np.array([max(st[0]) for st in basis])
        for g2 in G2:
            Hc = ks.hamiltonian(cas, W, g2).toarray()
            assert np.max(np.abs(Hc.imag)) < 1e-12
            Hd = Hc.real
            w, V = np.linalg.eigh(Hd)
            E = w[:3]
            exact[(name, g2)] = (E, V[:, :3], Hd, basis)
            D = E[1] - E[0]
            for scheme, lab in (("depth", depth), ("jmax", tjmax_st)):
                for k in range(1, int(lab.max()) + 1):
                    idx = np.where(lab <= k)[0]
                    if scheme == "depth" and k > 12:
                        break
                    if scheme == "jmax" and k > 12:
                        break
                    e = ritz(Hd, idx)
                    rows.append((name, g2, scheme, k, len(idx), e[0], e[1], e[1] - e[0],
                                 E[0], E[1], D, (e[1] - e[0]) - D))
            # selected CI: keep the n basis states of largest weight in the exact E0 (resp. E1)
            # eigenvector, as a selected-CI / natural-configuration choice would
            for scheme, col in (("sel_E0", 0), ("sel_E1", 1)):
                order = np.argsort(-np.abs(V[:, col]))
                for n in (2, 3, 4, 6, 8, 12, 16, 24, 32):
                    if n >= len(basis):
                        break
                    e = ritz(Hd, np.sort(order[:n]))
                    rows.append((name, g2, scheme, n, n, e[0], e[1], e[1] - e[0],
                                 E[0], E[1], D, (e[1] - e[0]) - D))
    write_csv("ci_truncation.csv",
              ["lattice", "g2", "scheme", "k", "dim", "E0_ci", "E1_ci", "gap_ci", "E0", "E1",
               "gap_exact", "gap_error"], rows)
    return rows, exact


def summarize_ci(rows):
    log("== CI hierarchies (Ritz in truncated subspace; E_k(ci) >= E_k by HUM/min-max) ==")
    for name in ("plaq1", "ladder2"):
        for scheme in ("depth", "jmax", "sel_E0", "sel_E1"):
            for g2 in G2:
                r = [x for x in rows if x[0] == name and x[2] == scheme and x[1] == g2]
                errs = [x[11] for x in r]
                signs = "".join("+" if e > 1e-9 else ("-" if e < -1e-9 else "0") for e in errs)
                log(f"{name:8s} {scheme:6s} g2={g2:<4} exact gap {r[0][10]:.6f}; "
                    f"gap_ci - gap by k=1..: " + " ".join(f"{e:+.3g}" for e in errs[:9])
                    + f"   signs {signs[:12]}  (k = " + ",".join(str(x[3]) for x in r[:9]) + ")")
    log()


# ------------------------------------------------------------------ HF (product of plaquette states)
def hf_runs(exact):
    log("== HF: best gauge-invariant product state on the ladder = chi-product, E_HF = 2 E0(plaq) ==")
    rows = []
    for g2 in G2:
        Ep = exact[("plaq1", g2)][0]
        El = exact[("ladder2", g2)][0]
        ehf, dk = 2 * Ep[0], Ep[1] - Ep[0]
        rows.append((g2, ehf, El[0], ehf - El[0], dk, El[1] - El[0], dk - (El[1] - El[0])))
        log(f"g2={g2:<4} E_HF={ehf:.6f} E0={El[0]:.6f} (E_HF-E0={ehf - El[0]:.4g} >= 0); "
            f"Koopmans gap {dk:.6f} vs exact {El[1] - El[0]:.6f} (err {dk - (El[1] - El[0]):+.4g})")
    write_csv("hf_ladder.csv", ["g2", "E_HF", "E0", "E_HF_minus_E0", "gap_koopmans", "gap_exact",
                                "gap_error"], rows)
    log()


# ------------------------------------------------------------------ CC on the plaquette (character algebra)
class Chars:
    """Class functions on SU(2) as coefficient vectors over chi_{tj/2}, tj = 0..M (Haar-orthonormal)."""

    def __init__(self, M):
        self.M = M
        self.cas = np.array([t * (t + 2) / 4 for t in range(M + 1)])

    def mul(self, f, g):
        out = np.zeros(self.M + 1)
        for a in np.nonzero(f)[0]:
            for b in np.nonzero(g)[0]:
                for c in range(abs(a - b), a + b + 1, 2):
                    if c <= self.M:
                        out[c] += f[a] * g[b]
        return out

    def C(self, f):
        return self.cas * f

    def grad_dot(self, f, g):  # grad f . grad g = (f Cg + g Cf - C(fg))/2, C = -Laplacian
        return 0.5 * (self.mul(f, self.C(g)) + self.mul(g, self.C(f)) - self.C(self.mul(f, g)))


def hbar(ch, s, f, g2):
    """e^{-S} H e^{S} f for H = (g^2/2) 4 C - (2/g^2) chi_{1/2}, S multiplication by s."""
    chi = np.zeros(ch.M + 1); chi[1] = 1.0
    kin = ch.C(f) - 2 * ch.grad_dot(s, f) + ch.mul(ch.C(s), f) - ch.mul(ch.grad_dot(s, s), f)
    return 2 * g2 * kin - (2 / g2) * ch.mul(chi, f)


def cc_runs(exact):
    log("== Coupled cluster on the single plaquette: |psi> = exp(S)|0>, S = sum_{tj=1..n} s chi ==")
    rows = []
    for g2 in G2:
        Eex = exact[("plaq1", g2)][0]
        D = Eex[1] - Eex[0]
        guess = None
        for n in range(1, 11):
            ch = Chars(6 * n + 12)
            one = np.zeros(ch.M + 1); one[0] = 1.0

            def eqs(x):
                s = np.zeros(ch.M + 1); s[1:n + 1] = x
                return hbar(ch, s, one, g2)[1:n + 1]
            x0 = np.zeros(n)
            if guess is not None:
                x0[:len(guess)] = guess
            x, info, ier, _ = fsolve(eqs, x0, full_output=True, xtol=1e-13)
            s = np.zeros(ch.M + 1); s[1:n + 1] = x
            res = np.max(np.abs(eqs(x)))
            ecc = hbar(ch, s, one, g2)[0]
            guess = x
            # EOM-CC: matrix of Hbar on characters tj = 0..n (left projection on the same set)
            Hb = np.zeros((n + 1, n + 1))
            for k in range(n + 1):
                f = np.zeros(ch.M + 1); f[k] = 1.0
                Hb[:, k] = hbar(ch, s, f, g2)[:n + 1]
            ev = np.sort(np.linalg.eigvals(Hb).real)
            geom = ev[1] - ev[0]
            rows.append((g2, n, ecc, Eex[0], ecc - Eex[0], ev[0], ev[1], geom, D, geom - D,
                         res, int(ier == 1)))
        r = [x for x in rows if x[0] == g2]
        log(f"g2={g2:<4} E0={Eex[0]:.6f} gap={D:.6f}; E_CC-E0 n=1..10: "
            + " ".join(f"{x[4]:+.2e}" for x in r))
        log(f"          EOM gap-err n=1..10: " + " ".join(f"{x[9]:+.2e}" for x in r)
            + f"  (max residual {max(x[10] for x in r):.1e}, converged {sum(x[11] for x in r)}/10)")
    write_csv("cc_plaquette.csv", ["g2", "n", "E_CC", "E0", "E_CC_minus_E0", "eom_E0", "eom_E1",
                                   "eom_gap", "gap_exact", "gap_error", "residual", "converged"], rows)
    log()


# ------------------------------------------------------------------ MPS-like Schmidt truncation on the ladder
def mps_runs(exact):
    log("== Ladder 'DMRG centre' Ritz: keep D_total Schmidt multiplets across the shared link ==")
    rows = []
    for g2 in G2:
        E, V, Hd, basis = exact[("ladder2", g2)]
        D = E[1] - E[0]
        a = np.array([st[0][0] for st in basis]); b = np.array([st[0][1] for st in basis])
        c = np.array([st[0][5] for st in basis])
        nA = a.max() + 1
        for wts, tag in (((1.0, 0.0), "E0-target"), ((0.5, 0.5), "averaged"), ((0.0, 1.0), "E1-target")):
            # sector-resolved left/right density matrices
            cand = []  # (weight, side, c, vector)
            for cc in np.unique(c):
                sel = np.where(c == cc)[0]
                rhoL = np.zeros((nA, nA)); rhoR = np.zeros((nA, nA))
                for wgt, col in zip(wts, (0, 1)):
                    M = np.zeros((nA, nA)); M[a[sel], b[sel]] = V[sel, col].real
                    rhoL += wgt * M @ M.T; rhoR += wgt * M.T @ M
                for side, rho in (("L", rhoL), ("R", rhoR)):
                    w, U = np.linalg.eigh(rho)
                    for i in range(nA):
                        cand.append((w[i], side, cc, U[:, i]))
            for Dk in (1, 2, 3, 4, 6, 8, 12, 16, 24):
                keepL = sorted([x for x in cand if x[1] == "L"], key=lambda x: -x[0])[:Dk]
                keepR = sorted([x for x in cand if x[1] == "R"], key=lambda x: -x[0])[:Dk]
                vecs = []
                for (_, _, cl, uL) in keepL:
                    for (_, _, cr, uR) in keepR:
                        if cl != cr:
                            continue
                        sel = np.where(c == cl)[0]
                        v = np.zeros(len(basis)); v[sel] = uL[a[sel]] * uR[b[sel]]
                        if np.linalg.norm(v) > 1e-12:
                            vecs.append(v)
                if len(vecs) < 2:
                    continue
                Q, _ = np.linalg.qr(np.array(vecs).T)
                w = np.linalg.eigvalsh(Q.T @ Hd @ Q)
                rows.append((g2, tag, Dk, Q.shape[1], w[0], w[1], w[1] - w[0], E[0], E[1], D,
                             (w[1] - w[0]) - D))
        for tag in ("E0-target", "averaged", "E1-target"):
            r = [x for x in rows if x[0] == g2 and x[1] == tag]
            log(f"g2={g2:<4} {tag:9s} gap err by D={','.join(str(x[2]) for x in r)}: "
                + " ".join(f"{x[10]:+.2e}" for x in r))
    write_csv("mps_ladder.csv", ["g2", "target", "D_multiplets", "dim", "E0_D", "E1_D", "gap_D",
                                 "E0", "E1", "gap_exact", "gap_error"], rows)
    log()


# ------------------------------------------------------------------ Temple with an uncertified beta
def temple_runs(exact):
    log("== Temple lower bound for E1 with beta = CI-k E2 (upper bound on E2, wrong direction) ==")
    log("   Temple needs beta <= E2. Any certified beta is a lower bound on E2: the same kind of")
    log("   statement as the gap itself. 'false' = the computed 'bound' exceeds the true E1.")
    rows = []
    for name in ("plaq1", "ladder2"):
        lat = ks.single_plaquette() if name == "plaq1" else ks.ladder2()
        for g2 in G2:
            E, V, Hd, basis = exact[(name, g2)]
            depth = bfs_depth(sp.csr_matrix(np.abs(Hd - np.diag(np.diag(Hd))) > 0))
            for k in range(2, 13):
                idx = np.where(depth <= k)[0]
                if len(idx) < 3 or len(idx) == len(basis):
                    continue
                w, U = np.linalg.eigh(Hd[np.ix_(idx, idx)])
                phi = np.zeros(len(basis)); phi[idx] = U[:, 1]
                Hphi = Hd @ phi
                eta = phi @ Hphi
                var = Hphi @ Hphi - eta ** 2
                for beta, btag in ((w[2], "beta=E2_ci"), (E[2], "beta=E2_exact")):
                    lb = eta - var / (beta - eta) if beta > eta else -np.inf
                    rows.append((name, g2, k, btag, eta, var, beta, lb, E[1], lb - E[1],
                                 "false" if lb > E[1] + 1e-9 else "valid"))
    write_csv("temple.csv", ["lattice", "g2", "k", "beta_choice", "eta", "variance", "beta",
                             "E1_lower", "E1", "lower_minus_E1", "status"], rows)
    nfalse = sum(1 for r in rows if r[10] == "false" and r[3] == "beta=E2_ci")
    ntot = sum(1 for r in rows if r[3] == "beta=E2_ci")
    nfalse_ex = sum(1 for r in rows if r[10] == "false" and r[3] == "beta=E2_exact")
    log(f"beta = CI E2: {nfalse}/{ntot} 'bounds' false; beta = exact E2: {nfalse_ex} false (theorem)")
    for r in rows:
        if r[10] == "false" and r[3] == "beta=E2_ci":
            log(f"   false: {r[0]} g2={r[1]} k={r[2]} 'E1 >=' {r[7]:.6f} but E1 = {r[8]:.6f}")
    log()


# ------------------------------------------------------------------ 1D chains: MPS does not imply a gap
def kron_op(ops, L):
    out = sp.identity(1, format="csr")
    for i in range(L):
        out = sp.kron(out, ops.get(i, sp.identity(2, format="csr")), format="csr")
    return out


def tfim(L, h):
    X = sp.csr_matrix([[0, 1], [1, 0]], dtype=float); Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=float)
    H = sp.csr_matrix((2 ** L, 2 ** L))
    for i in range(L - 1):
        H = H - kron_op({i: Z, i + 1: Z}, L)
    for i in range(L):
        H = H - h * kron_op({i: X}, L)
    return H


def heis_ferro(L):  # H = -sum S_i.S_{i+1}, periodic
    Sp = sp.csr_matrix([[0, 1], [0, 0]], dtype=float); Sm = Sp.T.tocsr()
    Sz = sp.csr_matrix([[0.5, 0], [0, -0.5]])
    H = sp.csr_matrix((2 ** L, 2 ** L))
    for i in range(L):
        j = (i + 1) % L
        H = H - kron_op({i: Sz, j: Sz}, L) - 0.5 * (kron_op({i: Sp, j: Sm}, L) + kron_op({i: Sm, j: Sp}, L))
    return H


def tt_svd(psi, L, D):
    """Left-to-right SVD truncation to bond dimension D; returns the normalized truncated vector."""
    rest = psi.reshape(1, -1)
    cores = []
    r = 1
    for i in range(L - 1):
        rest = rest.reshape(r * 2, -1)
        U, s, Vt = np.linalg.svd(rest, full_matrices=False)
        k = min(D, len(s))
        cores.append(U[:, :k].reshape(r, 2, k))
        rest = s[:k, None] * Vt[:k]
        r = k
    cores.append(rest.reshape(r, 2, 1))
    v = cores[0]
    for c in cores[1:]:
        v = np.tensordot(v, c, axes=([-1], [0]))
    v = v.reshape(-1)
    return v / np.linalg.norm(v)


def chain_runs():
    from scipy.sparse.linalg import eigsh
    log("== 1D chains: small-D MPS approximations exist for gapless finite chains ==")
    rows = []
    for model, h in (("TFIM_critical_h1", 1.0), ("TFIM_gapped_h2", 2.0)):
        for L in (8, 12, 16):
            H = tfim(L, h)
            w, V = eigsh(H, k=2, which="SA", tol=1e-12)
            o = np.argsort(w); w = w[o]; psi = V[:, o[0]]
            for D in (1, 2, 4, 8):
                v = tt_svd(psi, L, D)
                e = v @ (H @ v)
                rows.append((model, L, D, w[0], e, (e - w[0]) / abs(w[0]), w[1] - w[0]))
            r = [x for x in rows if x[0] == model and x[1] == L]
            log(f"{model:17s} L={L:2d} gap={w[1] - w[0]:.4f}  rel E err D=1,2,4,8: "
                + " ".join(f"{x[5]:.1e}" for x in r))
    for L in (8, 12):
        H = heis_ferro(L)
        w = np.sort(eigsh(H, k=4, which="SA", tol=1e-12, return_eigenvectors=False))
        up = np.zeros(2 ** L); up[0] = 1.0
        e = up @ (H @ up)
        E0 = -L / 4
        # first excitation above the (2S+1)-fold ground multiplet: magnon at k = 2 pi/L
        n = 2 ** L
        # restrict to Sz = L/2 - 1 (one flipped spin) where the magnon lives
        idx = [1 << i for i in range(L)]
        H1 = H[idx][:, idx].toarray()
        m = np.sort(np.linalg.eigvalsh(H1))
        gap = m[m > E0 + 1e-9][0] - E0
        rows.append(("Heisenberg_ferro", L, 1, E0, e, (e - E0) / abs(E0), gap))
        log(f"Heisenberg_ferro  L={L:2d} product state |up...up> energy {e:.4f} = E0 {E0:.4f} (D=1 exact); "
            f"gap {gap:.5f} vs 1-cos(2pi/L) = {1 - np.cos(2 * np.pi / L):.5f}  (dim {n})")
    write_csv("mps_chains.csv", ["model", "L", "D", "E0", "E_D", "rel_err", "gap_finite_L"], rows)
    log()


def main():
    rows, exact = ci_runs()
    summarize_ci(rows)
    hf_runs(exact)
    cc_runs(exact)
    mps_runs(exact)
    temple_runs(exact)
    chain_runs()
    (OUT / "summary.txt").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
