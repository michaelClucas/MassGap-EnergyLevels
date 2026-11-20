"""2x2 torus: what does the vacuum-seeded BFS leave out?

Seeds the same BFS (matrix elements from src/ks_su2.py, validated pointwise in 003) from
nontrivial 't Hooft electric-flux states: a fundamental Polyakov loop winding in x (flux (1,0)),
one in y (0,1), and one in each (1,1). Reports the lowest energy in each flux sector relative to
the zero-flux ground state E0. Also checks whether a j=1 winding loop (center-trivial) lies in the
vacuum-connected basis.
"""
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.ks_su2 import periodic, build, hamiltonian, lowest, plaquette_moves, intertwiners  # noqa: E402


def build_from(lat, tjmax, seed):
    index = {seed: 0}
    basis = [seed]
    rows, cols, vals = [], [], []
    q = deque([seed])
    while q:
        st = q.popleft()
        i = index[st]
        for p in range(len(lat.plaqs)):
            for nst, a in plaquette_moves(lat, st, p, tjmax):
                if nst not in index:
                    index[nst] = len(basis)
                    basis.append(nst)
                    q.append(nst)
                rows.append(index[nst]); cols.append(i); vals.append(a)
    n = len(basis)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    cas = np.array([sum(t * (t + 2) / 4 for t in s[0]) for s in basis])
    return basis, cas, W, index


def seed_state(lat, tj_on):
    tjs = [0] * len(lat.links)
    for l, t in tj_on.items():
        tjs[l] = t
    for v in range(lat.nv):
        assert intertwiners(lat.vertex_legs(v, tjs)).shape[0] >= 1, "seed violates Gauss law"
    return tuple(tjs), tuple([0] * lat.nv)


if __name__ == "__main__":
    lat = periodic(2, 2)
    # link order from periodic(): for y, for x: x-link then y-link -> x(0,0)=0, y(0,0)=1, x(1,0)=2, y(1,0)=3,
    # x(0,1)=4, y(0,1)=5, x(1,1)=6, y(1,1)=7
    seeds = {
        "flux(1,0)": {0: 1, 2: 1},            # fundamental loop along x at y=0
        "flux(0,1)": {1: 1, 5: 1},            # fundamental loop along y at x=0
        "flux(1,1)": {0: 1, 2: 1, 1: 1, 5: 1},
    }
    g2s = (1.0, 2.0, 3.0, 5.0, 10.0)
    for tjmax in (1, 2, 3):
        t0 = time.time()
        b0, c0, W0 = build(lat, tjmax)
        idx0 = {s: i for i, s in enumerate(b0)}
        E0 = {g2: lowest(hamiltonian(c0, W0, g2), k=2) for g2 in g2s}
        print(f"j_max={tjmax/2}: zero-flux dim {len(b0)}", flush=True)
        if tjmax >= 2:
            j1 = seed_state(lat, {0: 2, 2: 2})
            print(f"  j=1 winding loop along x in the vacuum-connected basis: {j1 in idx0}")
        for name, on in seeds.items():
            s = seed_state(lat, on)
            b, c, W, idx = build_from(lat, tjmax, s)
            overlap = sum(1 for st in b if st in idx0)
            line = []
            for g2 in g2s:
                e = lowest(hamiltonian(c, W, g2), k=1)[0]
                line.append(f"g2={g2}: {e - E0[g2][0]:.5f} (zero-flux gap {E0[g2][1]-E0[g2][0]:.5f})")
            print(f"  {name}: dim {len(b)}, states shared with zero-flux basis {overlap}")
            print("    E_flux - E0:  " + "; ".join(line), flush=True)
        print(f"  [{time.time()-t0:.0f}s]", flush=True)
    print("strong-coupling expectation: one fundamental winding loop of 2 links costs (g^2/2)*2*(3/4) = 3g^2/4,")
    print("half the zero-flux one-plaquette cost 3g^2/2.")
