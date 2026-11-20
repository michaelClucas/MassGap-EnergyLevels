"""Fit the constituent-gluon model on a few states, predict the rest. Writes fits.csv and fit.log text to stdout.

Run from repo root:  python tools/pin.py .venv/Scripts/python.exe work/20261120-013-resultB-constituent-gluon/fit.py
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as M  # noqa: E402

HERE = Path(__file__).resolve().parent
MODEL_ERR = 0.05  # stated model-error budget: 5% of the lattice mass, added in quadrature

# ------------------------------------------------------------- lattice data
SU3_4D = M.load("glueball_su3_4d_athenodorou_teper_2020.csv", 3, 4, "AT2020")
INF_4D = M.load("glueball_sun_4d_athenodorou_teper_2021.csv", "inf", 4, "AT2021")
SU3_3D = M.load("glueball_sun_3d_athenodorou_teper_2016.csv", 3, 3, "AT2016")
INF_3D = M.load("glueball_sun_3d_athenodorou_teper_2016.csv", "inf", 3, "AT2016")
g2N_3 = 1.0 / M.load_scalar("glueball_sun_3d_athenodorou_teper_2016.csv", 3, "sqrt_sigma/g2N")[0]
g2N_inf = 1.0 / M.load_scalar("glueball_sun_3d_athenodorou_teper_2016.csv", "inf", "sqrt_sigma/g2N")[0]

rows = []


def model_level(ch, jpc, lev):
    v = ch.get(jpc, [])
    return v[lev] if lev < len(v) else (float("nan"), "no state in model")


def compare(setting, variant, ch, lat, fitset, k):
    chi_l, chi_b, npred, nwithin = 0.0, 0.0, 0, 0
    print(f"\n== {setting} / {variant}  (k = {k}; fit set {sorted(fitset)})")
    print(f"{'JPC':5s} {'lev':>3s} {'lattice':>12s} {'model':>8s} {'pull_lat':>9s} {'pull_5%':>8s}  role  model state")
    for (jpc, lev), (val, err) in sorted(lat.items(), key=lambda kv: kv[1][0]):
        E, lab = model_level(ch, jpc, lev)
        role = "fit" if (jpc, lev) in fitset else "pred"
        if math.isnan(E):
            pl = pb = float("nan")
        else:
            pl = (E - val) / err
            pb = (E - val) / math.hypot(err, MODEL_ERR * val)
        if role == "pred":
            npred += 1
            if not math.isnan(E):
                chi_l += pl ** 2
                chi_b += pb ** 2
                nwithin += abs(pl) <= 1.0
        print(f"{jpc:5s} {lev:3d} {val:7.3f}({err:5.3f}) {E:8.3f} {pl:9.1f} {pb:8.2f}  {role:4s}  {lab}")
        rows.append(dict(setting=setting, variant=variant, k=k, JPC=jpc, level=lev, lattice=val, err=err,
                         model=round(E, 4) if not math.isnan(E) else "", role=role,
                         pull_lattice=round(pl, 2) if not math.isnan(pl) else "",
                         pull_budget5pct=round(pb, 3) if not math.isnan(pb) else "", model_state=lab))
    npm = sum(1 for (j, l) in lat if (j, l) not in fitset and not math.isnan(model_level(ch, j, l)[0]))
    print(f"predictions: {npred} lattice states ({npm} with a model counterpart); within 1 sigma(lat): {nwithin}; "
          f"chi2/pred lattice-only = {chi_l / max(npm, 1):.1f}, with 5% budget = {chi_b / max(npm, 1):.2f}")
    return dict(chi_l=chi_l / max(npm, 1), chi_b=chi_b / max(npm, 1), nwithin=nwithin, npm=npm, npred=npred)


def fit4d(variant, kin, free, fitset, start):
    """free: names among m, kappa, C0; others fixed at C0=0."""
    def unpack(x):
        p = dict(m=None, kappa=None, C0=0.0)
        p.update(zip(free, x))
        return (p["m"], p["kappa"], p["C0"])

    def res(x):
        ch = M.by_channel(M.spectrum4d(unpack(x), N=3, kin=kin))
        return [(model_level(ch, j, l)[0] - SU3_4D[(j, l)][0]) / SU3_4D[(j, l)][1] for (j, l) in fitset]

    lo = dict(m=0.3, kappa=0.0, C0=-15.0)
    hi = dict(m=5.0, kappa=4.0, C0=15.0)
    sol = least_squares(res, start, bounds=([lo[f] for f in free], [hi[f] for f in free]), x_scale=1.0,
                        diff_step=1e-3)
    p = unpack(sol.x)
    print(f"\n#### 4D fit {variant}: m/sqrt(sigma) = {p[0]:.4f}, kappa = N alpha_s = {p[1]:.4f} "
          f"(alpha_s = {p[1] / 3:.4f}), C0/sqrt(sigma) = {p[2]:.4f}; fit residual chi2 = {2 * sol.cost:.3g}")
    ch3 = M.by_channel(M.spectrum4d(p, N=3, kin=kin))
    s3 = compare("4D SU(3) AT2020", variant, ch3, SU3_4D, fitset, len(free))
    # large N at fixed 't Hooft coupling: kappa = lambda/(4 pi) fixed, m and C0 fixed, sigma_A/sigma -> 2
    chinf = M.by_channel(M.spectrum4d(p, N=10 ** 6, kin=kin))
    sinf = compare("4D SU(inf) AT2021", variant + " (same params, lambda fixed)", chinf, INF_4D, set(), 0)
    return p, s3, sinf


def fit3d(variant, kin, fitset, start):
    def res(x):
        ch = M.by_channel(M.spectrum3d(tuple(x), 3, g2N_3, kin=kin))
        return [(model_level(ch, j, l)[0] - SU3_3D[(j, l)][0]) / SU3_3D[(j, l)][1] for (j, l) in fitset]

    sol = least_squares(res, start, bounds=([0.2, -20], [6, 20]), diff_step=1e-3)
    p = tuple(sol.x)
    print(f"\n#### 2+1D fit {variant}: m/sqrt(sigma) = {p[0]:.4f}, C0/sqrt(sigma) = {p[1]:.4f}; "
          f"g^2N/sqrt(sigma) = {g2N_3:.4f} (lattice input); fit chi2 = {2 * sol.cost:.3g}")
    ch3 = M.by_channel(M.spectrum3d(p, 3, g2N_3, kin=kin))
    s3 = compare("3D SU(3) AT2016", variant, ch3, SU3_3D, fitset, 2)
    # SU(inf): m and the log-scale constant fixed in units of sqrt(sigma); g^2N from lattice at N=inf
    chinf = M.by_channel(M.spectrum3d(p, 10 ** 6, g2N_inf, kin=kin))
    sinf = compare("3D SU(inf) AT2016", variant + " (same m, C0)", chinf, INF_3D, set(), 0)
    return p, s3, sinf


if __name__ == "__main__":
    F2 = {("0++", 0), ("2++", 0)}
    F3 = {("0++", 0), ("2++", 0), ("0-+", 0)}
    summary = []
    for kin in ("rel", "nr"):
        p, a, b = fit4d(f"{kin} k=2 (m,kappa; C0=0)", kin, ["m", "kappa"], F2, [1.0, 0.6])
        summary.append(("4D", kin, 2, p, a, b))
        p, a, b = fit4d(f"{kin} k=3 (m,kappa,C0)", kin, ["m", "kappa", "C0"], F3, [1.0, 0.6, 0.0])
        summary.append(("4D", kin, 3, p, a, b))
        p, a, b = fit3d(f"{kin} k=2 (m,C0)", kin, F2, [1.0, 0.0])
        summary.append(("3D", kin, 2, p, a, b))
    with open(HERE / "fits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n#### summary (per variant: SU(3) then SU(inf); chi2 per predicted state with a model counterpart)")
    for dim, kin, k, p, a, b in summary:
        print(f"{dim} {kin} k={k} params={np.round(p, 4).tolist()}  SU3: {a['nwithin']}/{a['npm']} within 1sig, "
              f"chi2/pred lat={a['chi_l']:.1f} 5%={a['chi_b']:.2f} | SUinf: {b['nwithin']}/{b['npm']}, "
              f"chi2/pred lat={b['chi_l']:.1f} 5%={b['chi_b']:.2f}")
