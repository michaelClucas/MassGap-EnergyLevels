"""Headline glueball ratios from data/*.csv with first-order error propagation. Stdlib only.

Run from repo root:  python work/20261120-002-glueball-data/ratios.py
Errors: err_stat and err_sys combined in quadrature; inputs treated as uncorrelated
(ratios of two m/sqrt_sigma values from one paper therefore slightly overstate the error,
since the common sqrt_sigma uncertainty partly cancels).
"""
import csv
import math
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")


def load(fn):
    with open(os.path.join(DATA, fn), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def get(rows, N, JPC, quantity, level=0, table=None):
    hits = [r for r in rows if r["N"] == str(N) and r["JPC"] == JPC and r["quantity"] == quantity
            and r["level"] == str(level) and (table is None or table in r["table"])]
    if len(hits) != 1:
        raise LookupError((N, JPC, quantity, level, table, len(hits)))
    r = hits[0]
    e = math.hypot(float(r["err_stat"] or 0), float(r["err_sys"] or 0))
    return float(r["value"]), e


def scalar(rows, quantity):
    hits = [r for r in rows if r["quantity"] == quantity and not r["JPC"]]
    r = hits[0]
    return float(r["value"]), math.hypot(float(r["err_stat"] or 0), float(r["err_sys"] or 0))


def div(a, b):
    v = a[0] / b[0]
    return v, v * math.hypot(a[1] / a[0], b[1] / b[0])


def fmt(x):
    return "%.3f(%.3f)" % x if x else "n/a"


def line(label, m0, m2, mps, msig):
    r2 = div(m2, m0) if m2 else None
    rps = div(mps, m0) if mps else None
    print("| %-46s | %-14s | %-14s | %-14s |" % (label, fmt(r2), fmt(rps), fmt(msig)))


def main():
    at20 = load("glueball_su3_4d_athenodorou_teper_2020.csv")
    mp = load("glueball_su3_4d_morningstar_peardon_1999.csv")
    ch = load("glueball_su3_4d_chen_2006.csv")
    ltw = load("glueball_sun_4d_lucini_teper_wenger_2004.csv")
    at21 = load("glueball_sun_4d_athenodorou_teper_2021.csv")
    t98 = load("glueball_sun_3d_teper_1998.csv")
    at16 = load("glueball_sun_3d_athenodorou_teper_2016.csv")
    r0s = scalar(at20, "r0*sqrt_sigma")  # AT2020 continuum r0*sqrt(sigma), used to convert r0 units

    print("| %-46s | %-14s | %-14s | %-14s |" % ("setting / source", "m(2++)/m(0++)", "m(0-+)/m(0++)", "m(0++)/sqrt(s)"))
    print("|" + "-" * 48 + "|" + "-" * 16 + "|" + "-" * 16 + "|" + "-" * 16 + "|")

    q = "m/sqrt_sigma"
    T17 = "Table 17"
    m0, m2, mps = (get(at20, 3, j, q, table=T17) for j in ("0++", "2++", "0-+"))
    line("4D SU(3) AT2020 (Tab 17)", m0, m2, mps, m0)
    for name, rows, tab in (("MP1999 (Tab VII)", mp, "Table VII"), ("Chen2006 (Tab XXI)", ch, "Table XXI")):
        m0, m2, mps = (get(rows, 3, j, "r0*m", table=tab) for j in ("0++", "2++", "0-+"))
        line("4D SU(3) " + name + "; m/sqrt(s) via AT2020 r0*sqrt(s)", m0, m2, mps, div(m0, r0s))
    m0, m2, mps = (get(at21, 3, j, q, table="Table 34") for j in ("0++", "2++", "0-+"))
    line("4D SU(3) AT2021 (Tab 34)", m0, m2, mps, m0)
    m0, m2 = (get(ltw, 3, j, q) for j in ("0++", "2++"))
    line("4D SU(3) LTW2004 (Tab 14)", m0, m2, None, m0)

    m0, m2, mps = (get(at21, "inf", j, q, table="Table 38") for j in ("0++", "2++", "0-+"))
    line("4D SU(inf) AT2021 (Tab 38)", m0, m2, mps, m0)
    m0, m2 = (get(ltw, "inf", j, q) for j in ("0++", "2++"))
    line("4D SU(inf) LTW2004 (Tab 14)", m0, m2, None, m0)

    m0, m2, mps = (get(at16, 3, j, q, table="Table B4") for j in ("0++", "2++", "0-+"))
    line("3D SU(3) AT2016 (Tab B4)", m0, m2, mps, m0)
    m0, m2, mps = (get(t98, 3, j, q, table="Table 28") for j in ("0++", "2++", "0-+"))
    line("3D SU(3) Teper1998 (Tab 28)", m0, m2, mps, m0)
    m0 = get(at16, "inf", "0++", q, table="Table B10")
    m2 = get(at16, "inf", "2++", q, table="Table B11")
    mps = get(at16, "inf", "0-+", q, table="Table B10")
    line("3D SU(inf) AT2016 (Tabs B10-B11)", m0, m2, mps, m0)
    m0, m2, mps = (get(t98, "inf", j, q, table="Table 34") for j in ("0++", "2++", "0-+"))
    line("3D SU(inf) Teper1998 (Tab 34)", m0, m2, mps, m0)


if __name__ == "__main__":
    main()
