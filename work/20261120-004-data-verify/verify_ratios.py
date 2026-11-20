"""Independent recomputation of the headline glueball ratios from data/ (stdlib only).

Written from scratch (does not import or copy work/20261120-002-glueball-data/ratios.py).
Errors: first and second quoted errors in quadrature per input, inputs uncorrelated,
linear propagation for ratios. Run from the repo root.
"""
import csv
import math

D = "data/"


def load(fname):
    with open(D + fname, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def pick(rows, N, jpc, quantity, table_prefix, level="0"):
    hits = [r for r in rows if r["N"] == N and r["JPC"] == jpc and r["quantity"] == quantity
            and (level is None or r["level"] == level) and r["table"].startswith(table_prefix)]
    assert len(hits) == 1, (N, jpc, quantity, table_prefix, len(hits))
    r = hits[0]
    e = math.hypot(float(r["err_stat"] or 0), float(r["err_sys"] or 0))
    return float(r["value"]), e


def ratio(a, b):
    v = a[0] / b[0]
    return v, v * math.hypot(a[1] / a[0], b[1] / b[0])


def fmt(x):
    return f"{x[0]:.4f} +- {x[1]:.4f}"


scales = load("string_tension_scales_su3_4d.csv")
conv = pick(scales, "3", "", "r0*sqrt_sigma", "Table 17", level=None)  # AT2020 r0*sqrt(sigma)

out = []


def row(label, m0, m2, mps, msig):
    out.append((label, ratio(m2, m0) if m2 else None, ratio(mps, m0) if mps else None, msig))


def add(label, fname, N, q, tab, tab2=None, conv_r0=False):
    rows = load(fname)
    m0 = pick(rows, N, "0++", q, tab)
    m2 = pick(rows, N, "2++", q, tab2 or tab)
    try:
        mps = pick(rows, N, "0-+", q, tab)
    except AssertionError:
        mps = None
    msig = ratio(m0, conv) if conv_r0 else m0
    row(label, m0, m2, mps, msig)


add("4D SU(3) AT2020 T17", "glueball_su3_4d_athenodorou_teper_2020.csv", "3", "m/sqrt_sigma", "Table 17")
add("4D SU(3) AT2021 T34", "glueball_sun_4d_athenodorou_teper_2021.csv", "3", "m/sqrt_sigma", "Table 34")
add("4D SU(3) MP1999 T7", "glueball_su3_4d_morningstar_peardon_1999.csv", "3", "r0*m", "Table VII", conv_r0=True)
add("4D SU(3) Chen2006 T21", "glueball_su3_4d_chen_2006.csv", "3", "r0*m", "Table XXI", conv_r0=True)
add("4D SU(3) LTW2004 T14", "glueball_sun_4d_lucini_teper_wenger_2004.csv", "3", "m/sqrt_sigma", "Table 14")
add("4D SU(8) LTW2004 T14", "glueball_sun_4d_lucini_teper_wenger_2004.csv", "8", "m/sqrt_sigma", "Table 14")
add("4D SU(8) AT2021 T35", "glueball_sun_4d_athenodorou_teper_2021.csv", "8", "m/sqrt_sigma", "Table 35")
add("4D SU(inf) AT2021 T38", "glueball_sun_4d_athenodorou_teper_2021.csv", "inf", "m/sqrt_sigma", "Table 38")
add("4D SU(inf) LTW2004 T14", "glueball_sun_4d_lucini_teper_wenger_2004.csv", "inf", "m/sqrt_sigma", "Table 14")
add("3D SU(3) AT2016 B4", "glueball_sun_3d_athenodorou_teper_2016.csv", "3", "m/sqrt_sigma", "Table B4")
add("3D SU(3) Teper1998 T28", "glueball_sun_3d_teper_1998.csv", "3", "m/sqrt_sigma", "Table 28")
add("3D SU(inf) AT2016 B10/B11", "glueball_sun_3d_athenodorou_teper_2016.csv", "inf", "m/sqrt_sigma", "Table B10", "Table B11")
add("3D SU(inf) Teper1998 T34", "glueball_sun_3d_teper_1998.csv", "inf", "m/sqrt_sigma", "Table 34")

print(f"{'dataset':28s} {'m2++/m0++':>18s} {'m0-+/m0++':>18s} {'m0++/sqrt(sigma)':>18s}")
res = {}
for label, r2, rps, ms in out:
    res[label] = (r2, rps, ms)
    print(f"{label:28s} {fmt(r2) if r2 else 'n/a':>18s} {fmt(rps) if rps else 'n/a':>18s} {fmt(ms):>18s}")


def tension(a, b, name):
    print(f"{name}: diff {a[0]-b[0]:+.3f}, {abs(a[0]-b[0])/math.hypot(a[1], b[1]):.2f} sigma")


print("\nDiscrepancies")
A = res["4D SU(3) AT2020 T17"]
tension(res["4D SU(3) MP1999 T7"][2], A[2], "m0++/sqrt(sigma) MP1999 vs AT2020")
tension(res["4D SU(3) Chen2006 T21"][2], A[2], "m0++/sqrt(sigma) Chen2006 vs AT2020")
tension(res["4D SU(3) LTW2004 T14"][2], A[2], "m0++/sqrt(sigma) LTW2004 vs AT2020")
tension(res["4D SU(3) LTW2004 T14"][0], A[0], "m2++/m0++ LTW2004 vs AT2020")
tension(res["4D SU(3) MP1999 T7"][0], A[0], "m2++/m0++ MP1999 vs AT2020")
tension(res["4D SU(8) LTW2004 T14"][2], res["4D SU(8) AT2021 T35"][2], "SU(8) m0++/sqrt(sigma) LTW2004 vs AT2021")
tension(res["4D SU(inf) LTW2004 T14"][2], res["4D SU(inf) AT2021 T38"][2], "SU(inf) m0++/sqrt(sigma) LTW2004 vs AT2021")

r0 = pick(scales, "3", "", "r0", "eqn", level=None)
hbarc = 197.3269804  # MeV fm, CODATA
print(f"\n1/r0 = {hbarc / r0[0]:.1f} +- {hbarc * r0[1] / r0[0]**2:.1f} MeV (paper: 418(5))")
s = conv[0] / r0[0] * hbarc
se = s * math.hypot(conv[1] / conv[0], r0[1] / r0[0])
print(f"sqrt(sigma) = 1.160/r0 = {s:.1f} +- {se:.1f} MeV (paper: 485(6))")
