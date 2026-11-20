"""Summarize rotor_certified.csv: per (L, bc, lmax) the smallest g with a certified gap lower bound > 0,
and the gap interval there; compares with the finite-ring rotor-tower scale g/L and exp(-2 pi/g)/g."""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for fn in ("rotor_certified.csv", "rotor_certified_large.csv"):
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        rows += [dict(r, src=fn) for r in csv.DictReader(open(p))]
groups = {}
for r in rows:
    groups.setdefault((int(r["L"]), r["bc"], int(r["lmax"]), r["src"][16:-4] or "main"), []).append(r)
print("L bc lmax dim | smallest certified g: gap in [lo, hi] | next g (fails): floor vs E1_hi | g/L, exp(-2pi/g)/g at smallest g")
for k, rs in sorted(groups.items()):
    ok = [r for r in rs if float(r["gap_lo"]) > 0]
    bad = [r for r in rs if float(r["gap_lo"]) <= 0]
    if ok:
        r = min(ok, key=lambda r: float(r["g"]))
        g = float(r["g"])
        s = f"g={g}: [{float(r['gap_lo']):.5g}, {float(r['gap_hi']):.5g}]"
        nb = [b for b in bad if float(b["g"]) < g]
        f = max(nb, key=lambda b: float(b["g"])) if nb else None
        fs = f"g={f['g']}: floor {float(f['floor']):.4g} vs E1_hi {float(f['E1_hi']):.4g}" if f else "-"
        print(*k, r["dim"], "|", s, "|", fs, "|", f"{g / k[0]:.3g}, {np.exp(-2 * np.pi / g) / g:.2e}")
    else:
        print(*k, rs[0]["dim"], "| none")
