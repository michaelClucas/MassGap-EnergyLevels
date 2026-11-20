"""Build momentum-blocked KS matrices on the 2x2 torus for every Z2 flux sector, complete basis
up to spin tjmax/2, and save them to artifacts/012-torus-sectors/ (gitignored). Slow (pure Python
plaquette moves): ~10-20 min for tjmax = 4 with 2 processes.
Usage (repo root): python work/20261120-012-torus-sectors/build.py <tjmax> <fx,fy> [<fx,fy> ...]
"""
import pickle
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import torus_sym as T  # noqa: E402

ART = ROOT / "artifacts" / "012-torus-sectors"


def main():
    ART.mkdir(parents=True, exist_ok=True)
    tj = int(sys.argv[1])
    for f in sys.argv[2:]:
        fl = tuple(int(x) for x in f.split(","))
        t = time.time()
        log = lambda m: print(m, flush=True)  # noqa: E731
        log(f"flux {fl} tjmax {tj}: building")
        blocks, n = T.build_blocks(fl, tj, nproc=2, log=log)
        with open(ART / f"blocks_f{fl[0]}{fl[1]}_tj{tj}.pkl", "wb") as fh:
            pickle.dump(dict(flux=fl, tjmax=tj, dim=n, blocks=blocks), fh)
        log(f"flux {fl} tjmax {tj}: dim {n}, blocks {[len(b[0]) for b in blocks.values()]}, "
            f"{time.time() - t:.0f}s")


if __name__ == "__main__":
    main()
