"""Extract from the work-012 block pickles (artifacts/012-torus-sectors, gitignored) only what the
C9 certificates use: per (flux, k) block, Casimir diagonal on P, W on P x P and W on shell x P.
Zero flux: P = spins <= 3/2, shell = spin 2 (from blocks_f00_tj4). Flux (1,0), (1,1): P = spins <= 1,
shell = spin 3/2 (from blocks_f??_tj3). Writes artifacts/024-inertia/torus_blocks.npz.
Usage (repo root): python tools/pin.py .venv/Scripts/python.exe -u work/20261120-024-inertia-certification/extract_torus.py
"""
import gc
import pickle
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
ART = ROOT / "artifacts" / "012-torus-sectors"
OUT = ROOT / "artifacts" / "024-inertia"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    save = {}
    for fl, top in (((0, 0), 4), ((1, 0), 3), ((1, 1), 3)):
        with open(ART / f"blocks_f{fl[0]}{fl[1]}_tj{top}.pkl", "rb") as fh:
            data = pickle.load(fh)
        for k, (reps, cas, W, S) in data["blocks"].items():
            maxsp = np.array([max(r[0]) for r in reps])
            iP, ish = np.where(maxsp <= top - 1)[0], np.where(maxsp == top)[0]
            W = sp.csr_matrix(W)
            key = f"f{fl[0]}{fl[1]}_k{k[0]}{k[1]}"
            save[key + "_cas"] = np.asarray(cas)[iP]
            WP = W[iP][:, iP].tocoo(); WS = W[ish][:, iP].tocoo()
            for nm, M in (("WP", WP), ("WS", WS)):
                save[f"{key}_{nm}_r"], save[f"{key}_{nm}_c"], save[f"{key}_{nm}_v"] = M.row, M.col, M.data
                save[f"{key}_{nm}_shape"] = np.array(M.shape)
            print(f"{key}: P {len(iP)}, shell {len(ish)}", flush=True)
        del data
        gc.collect()
    np.savez_compressed(OUT / "torus_blocks.npz", **save)
    print("saved", flush=True)


if __name__ == "__main__":
    main()
