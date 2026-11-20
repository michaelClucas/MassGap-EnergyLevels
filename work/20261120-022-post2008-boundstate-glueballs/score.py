"""022: score post-2008 bound-state glueball models (models.csv) against continuum lattice data in data/.

Conventions (same as skeptic 019):
- data loading: `load` imported from work/20261120-019-skeptic-c12/model.py (m/sqrt_sigma rows only; lower bounds excluded).
- ratio test: model m/m(0++) vs lattice m/m(0++), lattice ratio error from m/sqrt(sigma) errors in quadrature
  (correlations ignored), pull_lat = (r - lv)/le, pull_5% = (r - lv)/hypot(le, 0.05 lv).
  Copied from 019 compare.py (loop over STATES); summary metrics n, #|pull|<=1, chi2/n copied from 019 rank.py.
- absolute test (models that give m/sqrt(sigma)): same pull formulas on m/sqrt(sigma) directly.
No lattice numbers are typed here. hbar*c is a physical constant, used only to convert FBSE MeV -> r0 units.
"""
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "20261120-019-skeptic-c12"))
import model as M019  # noqa: E402  (provides load(); DATA path = repo data/)

MODEL_ERR = 0.05
HBARC_MEV_FM = 197.3269804

LAT = {
    ("3", "4"): M019.load("glueball_su3_4d_athenodorou_teper_2020.csv", 3, 4, "AT2020"),
    ("inf", "4"): M019.load("glueball_sun_4d_athenodorou_teper_2021.csv", "inf", 4, "AT2021"),
    ("3", "3"): None,
    ("inf", "3"): M019.load("glueball_sun_3d_athenodorou_teper_2016.csv", "inf", 3, "AT2016"),
}
LAT_N3_AT2021 = M019.load("glueball_sun_4d_athenodorou_teper_2021.csv", 3, 4, "AT2021")
LATNAME = {("3", "4"): "AT2020 SU(3) 4D", ("inf", "4"): "AT2021 SU(inf) 4D", ("inf", "3"): "AT2016 SU(inf) 2+1D"}


def scalar(quantity, unit_filter=None):
    with open(M019.DATA / "string_tension_scales_su3_4d.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["source_id"] == "AT2020" and r["quantity"] == quantity and (unit_filter is None or r["units"] == unit_filter):
                return float(r["value"]), float(r["err_stat"] or 0)
    raise KeyError(quantity)


R0_FM, _ = scalar("r0", "fm")
R0SQS, R0SQS_E = scalar("r0*sqrt_sigma")

models = {}
with open(HERE / "models.csv", newline="") as f:
    for r in csv.DictReader(f):
        models.setdefault(r["model"], dict(N=r["N"], dim=r["dim"], units=r["units"], d={}))
        models[r["model"]]["d"][(r["state"], int(r["level"]))] = float(r["value"])

out = []


def summarize(tag, name, pulls):
    n = len(pulls)
    if n == 0:
        print(f"  {tag:10s} n=0")
        return
    pl = [p[0] for p in pulls]
    pb = [p[1] for p in pulls]
    print(f"  {tag:10s} n={n:2d} within1sig={sum(abs(x) <= 1 for x in pl):2d} "
          f"chi2/n lat={sum(x * x for x in pl) / n:9.1f} chi2/n 5%={sum(x * x for x in pb) / n:7.2f}")


def compare(name, test, mvals, lat, label):
    pulls = []
    for s in sorted(mvals):
        if s not in lat:
            continue
        r = mvals[s]
        lv, le = lat[s]
        pl = (r - lv) / le
        pb = (r - lv) / math.hypot(le, MODEL_ERR * lv)
        pulls.append((pl, pb))
        out.append(dict(model=name, test=test, lattice=label, JPC=s[0], level=s[1], model_value=round(r, 4),
                        lattice_value=round(lv, 4), lattice_err=round(le, 4), pull_lattice=round(pl, 2),
                        pull_budget5pct=round(pb, 2)))
        print(f"    {s[0]}{'*' * s[1]:3s} model {r:7.3f}  lattice {lv:7.3f}({le:.3f})  pull {pl:+7.1f}  pull5% {pb:+6.2f}")
    summarize(test, name, pulls)


def ratio_lat(lat):
    m0, e0 = lat[("0++", 0)]
    return {s: (v / m0, (v / m0) * math.hypot(e / v, e0 / m0)) for s, (v, e) in lat.items() if s != ("0++", 0)}


for name, md in models.items():
    key = (md["N"], md["dim"])
    lat = LAT.get(key)
    d = md["d"]
    print(f"\n== {name}  vs {LATNAME.get(key)}")
    if ("0++", 0) in d:
        base = d[("0++", 0)]
        rat = {s: v / base for s, v in d.items() if s != ("0++", 0)}
        compare(name, "ratio", rat, ratio_lat(lat), LATNAME[key])
    if md["units"] == "m/sqrt_sigma":
        compare(name, "absolute", d, lat, LATNAME[key])
    if md["units"] == "MeV" and name.startswith("FBSE") and ("0++", 0) in d:
        # FBSE masses quoted with r0 = 0.472 fm (AT2020 convention): m/sqrt(sigma) = m r0 / (r0 sqrt(sigma))
        v = d[("0++", 0)] * R0_FM / HBARC_MEV_FM / R0SQS
        print(f"  m(0++)/sqrt(sigma) from r0: {v:.3f}")
        compare(name, "m0/sqrtsig", {("0++", 0): v}, {("0++", 0): lat[("0++", 0)]}, LATNAME[key])

# N-dependence: BCMS2025 Casimir scaling, M(N)/sqrt(sigma) = eta(adj,N)/sqrt2 * M_inf, eta(3) = 3/2.
print("\n== Constituent-BCMS2025 extended to SU(3) via Casimir scaling (parameter-free given the SU(inf) fit)")
fac = 1.5 / math.sqrt(2)
b = models["Constituent-BCMS2025"]["d"]
compare("Constituent-BCMS2025->SU(3)", "absolute", {s: v * fac for s, v in b.items()}, LAT[("3", "4")], LATNAME[("3", "4")])
base = b[("0++", 0)]
compare("Constituent-BCMS2025->SU(3)", "ratio", {s: v / base for s, v in b.items() if s != ("0++", 0)},
        ratio_lat(LAT[("3", "4")]), LATNAME[("3", "4")])

# Pure N-dependence test of the Casimir law: lattice M(3)/M(inf) (both AT2021) vs 3/sqrt(8)
print("\n== Casimir-scaling law M(3)/M(inf) = 3/sqrt(8) vs AT2021 (ground states with both values)")
pred = 3 / math.sqrt(8)
lat_nd = {}
for s, (v3, e3) in LAT_N3_AT2021.items():
    if s in LAT[("inf", "4")] and s[1] == 0:
        vi, ei = LAT[("inf", "4")][s]
        lat_nd[s] = (v3 / vi, (v3 / vi) * math.hypot(e3 / v3, ei / vi))
compare("CasimirLaw(BCMS2025/HLLPV2017)", "N-ratio", {s: pred for s in lat_nd}, lat_nd, "AT2021 M(3)/M(inf)")
print("  (ratio-of-ratios: the law predicts m(2++)/m(0++) independent of N)")
for N in ("3", "inf"):
    L = LAT_N3_AT2021 if N == "3" else LAT[("inf", "4")]
    r = ratio_lat(L)[("2++", 0)]
    print(f"    AT2021 N={N}: m(2++)/m(0++) = {r[0]:.4f}({r[1]:.4f})")

with open(HERE / "comparison.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)

# FBSE with its own quoted (extrapolation-only) errors added in quadrature to the lattice ratio errors.
print("\n== FBSE-HFS2025 ratios with model errors (extrapolation only, uncorrelated quadrature) + lattice errors")
errs = {}
with open(HERE / "models.csv", newline="") as f:
    for r in csv.DictReader(f):
        if r["model"] == "FBSE-HFS2025":
            errs[(r["state"], int(r["level"]))] = (float(r["value"]), float(r["err"]))
m0, e0 = errs[("0++", 0)]
lr = ratio_lat(LAT[("3", "4")])
pc = []
for s, (v, e) in sorted(errs.items()):
    if s == ("0++", 0) or s not in lr:
        continue
    r = v / m0
    me = r * math.hypot(e / v, e0 / m0)
    p = (r - lr[s][0]) / math.hypot(lr[s][1], me)
    pc.append(p)
    print(f"    {s[0]}{'*' * s[1]:3s} {r:.3f}({me:.3f}) vs {lr[s][0]:.3f}  pull {p:+.2f}")
print(f"  n={len(pc)} within1sig={sum(abs(x) <= 1 for x in pc)} within2sig={sum(abs(x) <= 2 for x in pc)} "
      f"chi2/n={sum(x * x for x in pc) / len(pc):.2f}")
