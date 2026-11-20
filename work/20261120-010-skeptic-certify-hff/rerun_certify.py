"""Skeptic rerun of selected src/certify.py points (C5), written to this folder only.

1. Rerun plaq1 / ladder2 / periodic2x2 (j_max <= 1 only) points with the 006 code path and compare
   with 006 results.csv (reproducibility).
2. Soundness cross-check without new compute: for each lattice and g^2, every certified interval
   (all j_max) must intersect the tightest one; a lower bound above another row's upper bound would
   refute the method.
3. Rounding slack: on the plaquette, distance between 006 float bounds and the exact-rational bounds
   of plaq_exact.csv (how much of the 1e-11 width is floating-point margin).
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

HERE = Path(__file__).resolve().parent
REF = ROOT / "work/20261120-006-certified-gap-bounds/results.csv"


def rerun(name, ctor, top, tjs, g2s):
    lat = ctor()
    basis, cas, W = C.build_sector(lat, top)
    maxsp = np.array([max(st[0]) for st in basis])
    c_top, _ = C.casimir_floor(lat, top)
    out = []
    for tj in tjs:
        iP = np.where(maxsp <= tj)[0]
        c_Q, src = C.casimir_floor(lat, tj)
        kw = {}
        if tj < top:
            iQ1 = np.where(maxsp > tj)[0]
            kw = dict(shell=np.where(maxsp == tj + 1)[0], iQ1=iQ1, c_Q_top=c_top,
                      at_top_Q1=C.plaq_count_at_max(lat, [basis[i] for i in iQ1],
                                                    [True] * len(iQ1), top))
        else:
            kw = dict(at_max_P=C.plaq_count_at_max(lat, basis, maxsp <= tj, tj))
        for g2 in g2s:
            r = C.certify(cas, W, g2, len(lat.plaqs), iP, c_Q, **kw)
            out.append((name, g2, tj / 2, len(iP), r["lo"][0], r["hi"][0], r["lo"][1],
                        r["hi"][1], r["lo"][1] - r["hi"][0], r["hi"][1] - r["lo"][0], r["d"],
                        r["dsrc"], c_Q, r["method"]))
    return out


def load_ref():
    rows = list(csv.DictReader(open(REF)))
    for r in rows:
        for k in ("g2", "j_max", "E0_lo", "E0_hi", "E1_lo", "E1_hi", "Delta_lo", "Delta_hi"):
            r[k] = float(r[k])
    return rows


def main():
    t0 = time.time()
    rows = []
    rows += rerun("plaq1", ks.single_plaquette, 17, [4, 8, 16], [0.5, 1.0, 3.0, 10.0])
    rows += rerun("ladder2", ks.ladder2, 5, [2, 4], [1.0, 3.0, 10.0])
    rows += rerun("ladder2", ks.ladder2, 9, [8], [1.0])
    # torus: basis only up to j = 1 here (skeptic memory/time limit); j_max = 1 uses the CRUDE
    # shell, unlike 006 (exact shell from a j = 3/2 basis), so only j_max = 1/2 is a like-for-like
    # comparison; j_max = 1 must merely be a valid (looser) interval containing 006's.
    rows += rerun("periodic2x2", lambda: ks.periodic(2, 2), 2, [1, 2], [3.0, 5.0, 10.0])
    ref = load_ref()
    fields = ["lattice", "g2", "j_max", "dim", "E0_lo", "E0_hi", "E1_lo", "E1_hi", "Delta_lo",
              "Delta_hi", "d", "d_src", "c_Q", "method"]
    with open(HERE / "rerun_certify.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        w.writerows(rows)
    print("== 1. rerun vs 006 results.csv")
    for row in rows:
        name, g2, jm = row[0], row[1], row[2]
        m = [r for r in ref if r["lattice"] == name and r["g2"] == g2 and r["j_max"] == jm]
        tag = ""
        if m:
            r = m[0]
            diffs = [abs(a - b) for a, b in zip(row[4:10], [r[k] for k in fields[4:10]])
                     if np.isfinite(a) and np.isfinite(b)]
            tag = f"max |diff| vs 006 = {max(diffs):.2e}" if diffs else "no finite entries"
        print(f"  {name:12} g2={g2:5} j_max={jm:4}: Delta in [{row[8]:.10g}, {row[9]:.10g}] "
              f"d={row[10]:.4g} ({row[11]}) {tag}")
    print("== 2. cross-j_max consistency of all 006 rows (intervals must pairwise intersect)")
    bad = 0
    for key in sorted({(r["lattice"], r["g2"]) for r in ref}):
        rs = [r for r in ref if (r["lattice"], r["g2"]) == key]
        for q in ("E0", "E1", "Delta"):
            lo = max(r[q + "_lo"] for r in rs)
            hi = min(r[q + "_hi"] for r in rs)
            if lo > hi:
                bad += 1
                print(f"  INCONSISTENT {key} {q}: max lower {lo} > min upper {hi}")
    print(f"  inconsistencies: {bad}")
    print("== 3. plaquette float bounds vs exact rational bounds (positive slack = sound)")
    ex = {float(r["g2"]): r for r in csv.DictReader(open(HERE / "plaq_exact.csv"))}
    slack = []
    for r in ref:
        if r["lattice"] != "plaq1" or r["g2"] not in ex:
            continue
        e = ex[r["g2"]]
        for q in ("E0", "E1"):
            if np.isfinite(r[q + "_lo"]):
                slack.append((float(e[q + "_lo"]) - r[q + "_lo"], r["g2"], r["j_max"], q, "lo"))
            slack.append((r[q + "_hi"] - float(e[q + "_hi"]), r["g2"], r["j_max"], q, "hi"))
    slack.sort()
    for s in slack[:5]:
        print(f"  smallest slack {s[0]:.2e} at g2={s[1]} j_max={s[2]} {s[3]}_{s[4]}")
    print(f"done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
