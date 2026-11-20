"""Certified two-sided bounds on E0, E1 and Delta = E1 - E0 of the untruncated SU(2) KS Hamiltonian
on the single plaquette, open 2-plaquette ladder and periodic 2x2 lattice (2+1D, finite lattice,
fixed g^2), from truncated spin-network matrices. Method: README.md and src/certify.py.

Usage (repo root): python work/20261120-006-certified-gap-bounds/certify_gaps.py [--quick]
                   [--lattices plaq1,ladder2,periodic2x2]
Writes results_<lattice>.csv per lattice, merges all present into results.csv, appends to checks.txt.
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src import ks_su2 as ks  # noqa: E402
from src import certify as C  # noqa: E402

OUT = Path(__file__).resolve().parent
QUICK = "--quick" in sys.argv
G2 = [1.0, 3.0] if QUICK else [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
# lattice -> (constructor, top 2j built, list of 2j_max to certify)
# periodic2x2 at j_max = 3/2 (dim 9418): upper bounds only. A certified lower bound needs one
# inertia count (LDL^T) per level at n = 9418, which took > 3 min each here; see README.
PLAN = {
    "plaq1": (ks.single_plaquette, 9 if QUICK else 33, [2, 4, 8] if QUICK else [2, 4, 8, 16, 32]),
    "ladder2": (ks.ladder2, 5 if QUICK else 13, [2, 4] if QUICK else [2, 4, 8, 12]),
    "periodic2x2": (ks.periodic, 3, [1, 2, 3]),
}
SECTOR = {"plaq1": "gauge-invariant (all states)",
          "ladder2": "gauge-invariant (all states; both reflection parities)",
          "periodic2x2": "gauge-invariant zero Z2-flux; all momenta"}
FIELDS = ["lattice", "sector", "g2", "j_max", "dim", "E0_lo", "E0_hi", "E1_lo", "E1_hi",
          "Delta_lo", "Delta_hi", "d_floor", "d_source", "c_Q", "method"]


def log(msg, fh=None):
    print(msg, flush=True)
    if fh:
        fh.write(msg + "\n")


def run_lattice(name, chk):
    ctor, top, tjs = PLAN[name]
    lat = ctor(2, 2) if name == "periodic2x2" else ctor()
    t0 = time.time()
    basis, cas, W = C.build_sector(lat, top)
    maxsp = np.array([max(st[0]) for st in basis])
    log(f"[{name}] sector basis up to j={top/2}: dim {len(basis)} ({time.time()-t0:.1f}s)", chk)
    c_top, src_top = C.casimir_floor(lat, top)
    rows = []
    for tj in tjs:
        iP = np.where(maxsp <= tj)[0]
        vac = C.vacuum_component_size(basis, W, maxsp <= tj)
        c_Q, src = C.casimir_floor(lat, tj)
        log(f"[{name}] j_max={tj/2}: dim P {len(iP)}, vacuum-connected {vac}, "
            f"c_Q={c_Q:g} ({src})", chk)
        kw = {}
        if tj < top:
            kw["shell"] = np.where(maxsp == tj + 1)[0]
            iQ1 = np.where(maxsp > tj)[0]
            kw["iQ1"] = iQ1
            kw["at_top_Q1"] = C.plaq_count_at_max(lat, [basis[i] for i in iQ1],
                                                  [True] * len(iQ1), top)
            kw["c_Q_top"] = c_top
        else:
            kw["at_max_P"] = C.plaq_count_at_max(lat, basis, maxsp <= tj, tj)
        for g2 in G2:
            t1 = time.time()
            r = C.certify(cas, W, g2, len(lat.plaqs), iP, c_Q,
                          lower=not (name == "periodic2x2" and tj == top), **kw)
            lo, hi = r["lo"], r["hi"]
            if not np.isfinite(lo).all() and tj == top and name == "periodic2x2":
                r["method"] = "rayleigh-ritz-upper-only"
            row = dict(lattice=name, sector=SECTOR[name], g2=g2, j_max=tj / 2, dim=r["dim"],
                       E0_lo=lo[0], E0_hi=hi[0], E1_lo=lo[1], E1_hi=hi[1],
                       Delta_lo=lo[1] - hi[0], Delta_hi=hi[1] - lo[0], d_floor=r["d"],
                       d_source=r["dsrc"], c_Q=c_Q, method=r["method"])
            rows.append(row)
            log(f"  g2={g2:5}: E0 in [{lo[0]:.10g}, {hi[0]:.10g}]  E1 in [{lo[1]:.10g}, "
                f"{hi[1]:.10g}]  Delta in [{row['Delta_lo']:.8g}, {row['Delta_hi']:.8g}]  "
                f"d={r['d']:.4g} ({r['dsrc']})  {time.time()-t1:.1f}s", chk)
    with open(OUT / f"results_{name}{'_quick' if QUICK else ''}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(rows)
    return rows


def validate_plaq(chk):
    """Single plaquette: certified intervals must contain the iteration-1 j_max=40 ED values
    (which agree with the independent group-manifold solver to 1.2e-8)."""
    ref = {}
    with open(ROOT / "work/20261120-003-su2-ks-ed/results.csv") as f:
        for r in csv.DictReader(f):
            if r["lattice"] == "plaq1" and float(r["j_max"]) == 40.0:
                ref[float(r["g2"])] = (float(r["E0"]), float(r["E1"]))
    f = OUT / f"results_plaq1{'_quick' if QUICK else ''}.csv"
    bad = 0
    for r in csv.DictReader(open(f)):
        g2 = float(r["g2"])
        if g2 not in ref:
            continue
        e0, e1 = ref[g2]
        for e, a, b in ((e0, "E0_lo", "E0_hi"), (e1, "E1_lo", "E1_hi")):
            lo, hi = float(r[a]), float(r[b])
            # the reference itself is accurate to ~1e-12 (converged in j_max), allow that
            if not (lo - 1e-9 <= e <= hi + 1e-9):
                bad += 1
                log(f"  VIOLATION plaq1 g2={g2} j_max={r['j_max']}: {e} not in [{lo}, {hi}]", chk)
    log(f"[plaq1 validation] reference (j_max=40 ED) inside every certified interval: "
        f"{'yes' if bad == 0 else f'NO ({bad} violations)'}", chk)


def merge():
    rows = []
    for name in list(PLAN) + ([] if QUICK else ["periodic2x2_j32"]):
        f = OUT / f"results_{name}{'_quick' if QUICK else ''}.csv"
        if f.exists():
            rows += list(csv.DictReader(open(f)))
    with open(OUT / f"results{'_quick' if QUICK else ''}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(rows)


def main():
    names = list(PLAN)
    for a in sys.argv[1:]:
        if a.startswith("--lattices"):
            names = a.split("=", 1)[1].split(",")
    if QUICK:
        names = [n for n in names if n != "periodic2x2"]
    with open(OUT / f"checks{'_quick' if QUICK else ''}.txt", "a") as chk:
        for n in names:
            run_lattice(n, chk)
            if n == "plaq1":
                validate_plaq(chk)
    merge()


if __name__ == "__main__":
    main()
