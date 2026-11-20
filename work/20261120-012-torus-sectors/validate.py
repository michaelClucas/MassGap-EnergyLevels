"""Checks for src/torus_sym.py on the 2x2 torus (2+1D SU(2) KS, finite lattice):
1. direction-ordered legs give the same spectrum as ks.periodic(2,2) (zero flux, j_max = 1);
2. translations, sigma_x, sigma_d, rot90 are orthogonal and commute with W, per flux sector;
3. union of momentum-block spectra = full-sector spectrum (j_max = 1, all four flux sectors);
4. improved floor >= crude floor.
Usage (repo root): python work/20261120-012-torus-sectors/validate.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402
from src import certify as C  # noqa: E402
from src import torus_sym as T  # noqa: E402

FLUXES = [(0, 0), (1, 0), (0, 1), (1, 1)]


def main():
    tj = 2
    lat0 = ks.periodic(2, 2)
    lat = T.torus2x2()
    b0, c0, W0 = C.build_sector(lat0, tj, (0, 0))
    b1, c1, W1 = C.build_sector(lat, tj, (0, 0))
    g2 = 1.5
    e0 = np.linalg.eigvalsh(ks.hamiltonian(c0, W0, g2).toarray())
    e1 = np.linalg.eigvalsh(ks.hamiltonian(c1, W1, g2).toarray())
    print(f"[1] dim {len(b0)} vs {len(b1)}; max |spec diff| = {np.abs(e0 - e1).max():.2e}")
    ops = {"Tx": T.translation(1, 0), "Ty": T.translation(0, 1), "sigma_x": T.sigma_x(),
           "sigma_d": T.sigma_d(), "rot90": T.rot90()}
    worst = 0.0
    for fl in FLUXES:
        b, c, W = C.build_sector(lat, tj, fl)
        for name, op in ops.items():
            img_flux = fl[::-1] if name in ("sigma_d", "rot90") else fl
            bi = C.sector_states(lat, tj, img_flux)
            if img_flux != fl:
                continue
            R = T.op_matrix(lat, op, b)
            orth = abs(R.T @ R - np.eye(len(b))).max()
            comm = abs(R @ W - W @ R).max()
            worst = max(worst, orth, comm)
            print(f"[2] flux {fl} {name}: |R^T R - 1| = {orth:.1e}, |[R, W]| = {comm:.1e}")
    # translations are pure permutations
    b, c, W = C.build_sector(lat, tj, (0, 0))
    perm_ok = all(len(T.apply_op(lat, T.translation(1, 0), s)) == 1 and
                  abs(T.apply_op(lat, T.translation(1, 0), s)[0][1] - 1) < 1e-12 and
                  T.apply_op(lat, T.translation(1, 0), s)[0][0] == T.translate(s, (1, 0))
                  for s in b)
    print(f"[2] Tx acts as a pure permutation (apply_op == translate): {perm_ok}")
    for fl in FLUXES:
        b, c, W = C.build_sector(lat, tj, fl)
        full = np.sort(np.linalg.eigvalsh(ks.hamiltonian(c, W, g2).toarray()))
        blocks, n = T.build_blocks(fl, tj)
        sp_ = []
        for k, (reps, cas, Wk, S) in blocks.items():
            assert abs(Wk - Wk.T).max() < 1e-12
            if len(reps):
                sp_ += list(np.linalg.eigvalsh(ks.hamiltonian(cas, Wk, g2).toarray()))
        sp_ = np.sort(sp_)
        print(f"[3] flux {fl}: dim {n}, block dims {[len(v[0]) for v in blocks.values()]}, "
              f"max |spec diff| = {np.abs(full - sp_).max():.2e}")
    for g2 in (1.0, 1.5, 2.0, 3.0):
        f, t = T.split_floor(g2, 11.0)
        print(f"[4] g2={g2} c_Q=11: crude floor {0.5*g2*11 - 16/g2:.4f}, split floor {f:.4f} (t={t:.3f})")


if __name__ == "__main__":
    main()
