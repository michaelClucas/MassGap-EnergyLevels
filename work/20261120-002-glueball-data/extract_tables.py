"""Extract glueball / string-tension tables from arXiv LaTeX sources into data/*.csv.

Stdlib only. Numbers are parsed from the authors' LaTeX table source, never typed by hand.

Usage (from repo root):
    python work/20261120-002-glueball-data/extract_tables.py <SRC_DIR> [<OUT_DIR>]

<SRC_DIR> holds one subdirectory per paper, named by arXiv id with '/' -> '_'
(e.g. hep-lat_9901004), each containing the unpacked e-print from
https://arxiv.org/e-print/<arXiv id>. <OUT_DIR> defaults to data/.
"""
import csv
import os
import re
import sys
from decimal import Decimal

BS = chr(92)  # backslash
COLS = ["source_id", "N", "dim", "JPC", "quantity", "value", "err_stat", "err_sys",
        "units", "table", "level", "row_label", "note"]


# ---------------------------------------------------------------- LaTeX helpers
def read(path):
    with open(path, encoding="latin-1") as f:
        lines = [ln.replace(BS + "%", "").split("%")[0] for ln in f]
    return "".join(lines)


def table_block(tex, label):
    """Text of the tabular inside the table float that carries \\label{label}."""
    i = tex.index(BS + "label{" + label + "}")
    start = tex.rfind(BS + "begin{table", 0, i)
    end = tex.find(BS + "end{table", i)
    blk = tex[start:end]
    a = blk.index(BS + "begin{tabular}")
    b = blk.index(BS + "end{tabular}")
    return blk[a:b]


def rows(block):
    body = block.split("}", 1)[1]
    # drop column spec of tabular: first {...} group after begin{tabular}
    body = block[len(BS + "begin{tabular}"):]
    depth = 0
    for k, ch in enumerate(body):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                body = body[k + 1:]
                break
    out = []
    for r in body.split(BS + BS):
        r = r.replace(BS + "hline", "").strip()
        if not r:
            continue
        out.append([c.strip() for c in r.split("&")])
    return out


def clean(s):
    for t in (BS + "bb", BS + "hspace{1.25ex}", "$", "{", "}", " ", BS + ","):
        s = s.replace(t, "")
    return "".join(s.split())


VE = re.compile(r"^(\[?)(\d+\.?\d*)((?:\([\d.]+\)\)?)+)(\]?)(.*)$")


def parse_ve(cell):
    """'4.21(11)(4)' -> ('4.21', ['0.11','0.04'], rest). Returns None if no number."""
    c = clean(cell)
    lb = c.startswith(BS + "ge")
    if lb:
        c = c[len(BS + "ge"):]
    m = VE.match(c)
    if not m:
        return None
    val = m.group(2)
    errs = re.findall(r"\(([\d.]+)\)", m.group(3))
    nd = len(val.split(".")[1]) if "." in val else 0
    out = []
    for e in errs:
        if "." in e:
            out.append(e)
        else:
            out.append(str(Decimal(e).scaleb(-nd)))
    rest = (m.group(1) + m.group(4) + m.group(5)).replace(BS + "dagger", "dagger")
    rest = rest.replace("^", "")
    return val, out, rest, lb


def jpc_from(cell):
    """'$0^{++\\ast}$' / '$0^{*++}$' / '$0^{++}$ {gs}' -> ('0++', level, extra)."""
    c = clean(cell).replace(BS + "ast", "*")
    m = re.match(r"^(\d+)\^?(\**)([+-]{2})(\**)(.*)$", c)
    if not m:
        raise ValueError("bad JPC cell: " + cell)
    lev = len(m.group(2)) + len(m.group(4))
    extra = m.group(5)
    if extra == "gs":
        pass
    elif extra.startswith("ex"):
        lev = int(extra[2:] or 1)
    return m.group(1) + m.group(3), lev, extra


def row(**kw):
    r = {k: "" for k in COLS}
    r.update({k: str(v) for k, v in kw.items()})
    return r


def flags(rest):
    notes = []
    if "**" in rest:
        notes.append("author mark ** (significant uncertainty in J assignment)")
    elif "*" in rest:
        notes.append("author mark * (see source caption/text)")
    if "+" in rest:
        notes.append("author mark + (see source text)")
    if "dagger" in rest:
        notes.append("dagger: interpretation requires further study (per source)")
    if rest.startswith("[") or "]" in rest:
        notes.append("square brackets: average of N=12,16, not an extrapolation (per source)")
    return "; ".join(notes)


# ---------------------------------------------------------------- per-paper
def mp1999(src):
    tex = read(os.path.join(src, "hep-lat_9901004", "spectrum.tex"))
    out = []
    tab = "Table VII (label tab:finalmasses)"
    for r in rows(table_block(tex, "tab:finalmasses"))[1:]:
        jpc, lev, _ = jpc_from(r[0])
        other = clean(r[1]).replace(BS + "dots", "...")
        v1, e1, rest1, _ = parse_ve(r[2])
        v2, e2, rest2, _ = parse_ve(r[3])
        note = "other J possible: " + other if other else ""
        f = flags(rest1)
        note = "; ".join(x for x in (note, f) if x)
        out.append(row(source_id="MP1999", N=3, dim=4, JPC=jpc, quantity="r0*m", value=v1,
                       err_stat=e1[0], err_sys=e1[1], units="dimensionless", table=tab,
                       level=lev, row_label=clean(r[0]),
                       note="err_stat=continuum-extrapolation stat; err_sys=anisotropy" + ("; " + note if note else "")))
        out.append(row(source_id="MP1999", N=3, dim=4, JPC=jpc, quantity="m", value=v2,
                       err_stat=e2[0], err_sys=e2[1], units="MeV", table=tab, level=lev,
                       row_label=clean(r[0]),
                       note="err_stat=combined r0*m error; err_sys=from 1/r0=410(20) MeV" + ("; " + note if note else "")))
    tab = "Table VIII (label tab:massratios)"
    for r in rows(table_block(tex, "tab:massratios")):
        lab = clean(r[0]).replace(BS + "ast", "*").replace("[2mm]", "")
        m = re.match(r"^m\((\d)\^?(\*?)([+-]{2})\)/m\((\d[+-]{2})\)$", lab.replace("^", ""))
        if not m:
            m = re.match(r"^m\((\d)(\*?)([+-]{2})\)/m\((\d[+-]{2})\)$", lab.replace("^", ""))
        v, e, rest, _ = parse_ve(r[1].replace("[2mm]", ""))
        out.append(row(source_id="MP1999", N=3, dim=4, JPC=m.group(1) + m.group(3),
                       quantity="m/m(%s)" % m.group(4), value=v, err_stat=e[0],
                       units="dimensionless", table=tab, level=len(m.group(2)), row_label=lab))
    # scale convention quoted in caption of Table VII
    out.append(row(source_id="MP1999", N=3, dim=4, quantity="1/r0", value="410", err_stat="20",
                   units="MeV", table=tab.replace("VIII (label tab:massratios)", "VII caption (label tab:finalmasses)"),
                   note="scale convention used by the authors to convert r0*m to MeV"))
    return out


def chen2006(src):
    tex = read(os.path.join(src, "hep-lat_0510074", "paper.tex"))
    out = []
    tab = "Table XXI (label summary1; 21st table float in source)"
    for r in rows(table_block(tex, "summary1"))[1:]:
        jpc, lev, _ = jpc_from(r[0])
        v1, e1, _, _ = parse_ve(r[1])
        v2, e2, _, _ = parse_ve(r[2])
        n = "source typo '10.25(4))(10)' read as 10.25(4)(10)" if "))" in r[1] else ""
        out.append(row(source_id="Chen2006", N=3, dim=4, JPC=jpc, quantity="r0*m", value=v1,
                       err_stat=e1[0], err_sys=e1[1], units="dimensionless", table=tab, level=lev,
                       row_label=clean(r[0]),
                       note="err_stat=continuum-extrapolation stat; err_sys=1% anisotropy" + ("; " + n if n else "")))
        out.append(row(source_id="Chen2006", N=3, dim=4, JPC=jpc, quantity="m", value=v2,
                       err_stat=e2[0], err_sys=e2[1], units="MeV", table=tab, level=lev,
                       row_label=clean(r[0]),
                       note="err_stat=combined r0*m error; err_sys=from 1/r0=410(20) MeV"))
    tab = "Table VIII (label comparison; column 'this work')"
    for r in rows(table_block(tex, "comparison"))[1:]:
        jpc, lev, _ = jpc_from(r[1])
        v, e, _, _ = parse_ve(r[2])
        out.append(row(source_id="Chen2006", N=3, dim=4, JPC=jpc, quantity="r0*m", value=v,
                       err_stat=e[0], units="dimensionless", table=tab, level=lev,
                       row_label=clean(r[0]).replace("^", ""),
                       note="per cubic irrep R^PC (row_label), before combining into J"))
    out.append(row(source_id="Chen2006", N=3, dim=4, quantity="1/r0", value="410", err_stat="20",
                   units="MeV", table="Table XXI caption (label summary1)",
                   note="scale convention used by the authors"))
    return out


PC4 = ["++", "-+", "+-", "--"]


def at_jgrid(rows_, sid, N, tab, quantity, units, note_extra=""):
    out = []
    for r in rows_:
        lab = clean(r[0])
        m = re.match(r"^(\d)(gs|ex\d*)$", lab)
        if not m:
            continue
        J, lv = m.group(1), m.group(2)
        lev = 0 if lv == "gs" else int(lv[2:] or 1)
        for k, pc in enumerate(PC4):
            if k + 1 >= len(r):
                continue
            p = parse_ve(r[k + 1])
            if not p:
                continue
            v, e, rest, lb = p
            q = quantity + ("_lower_bound" if lb else "")
            note = "; ".join(x for x in (flags(rest), note_extra) if x)
            out.append(row(source_id=sid, N=N, dim=4, JPC=J + pc, quantity=q, value=v,
                           err_stat=e[0], units=units, table=tab, level=lev,
                           row_label=lab + " / P,C=" + pc, note=note))
    return out


def at2020(src):
    tex = read(os.path.join(src, "2007.06422", "paper.tex"))
    out = []
    rs = rows(table_block(tex, "table_MK_J"))
    out += at_jgrid(rs, "AT2020", 3, "Table 17 (label table_MK_J)", "m/sqrt_sigma", "dimensionless")
    for r in rs:
        if "r_0" in r[0]:
            v, e, _, _ = parse_ve(r[0].split("=")[1])
            out.append(row(source_id="AT2020", N=3, dim=4, quantity="r0*sqrt_sigma", value=v,
                           err_stat=e[0], units="dimensionless", table="Table 17 (label table_MK_J) last row",
                           note="continuum; r0 from interpolating formulae of refs HWalpha-1998, SNRS-2001 (per text)"))
    out += at_jgrid(rows(table_block(tex, "table_MGeV_J")), "AT2020", 3,
                    "Table 18 (label table_MGeV_J)", "m", "GeV",
                    "uses r0=0.472(5) fm; errors statistical only (per caption)")
    out.append(row(source_id="AT2020", N=3, dim=4, quantity="r0", value="0.472", err_stat="0.005",
                   units="fm", table="eqn label eqn_r0MeV (Sec. on physical units)",
                   note="lattice QCD (>=3 light flavours) value quoted from ref RS-2014; equals 1/(418(5) MeV)"))
    out.append(row(source_id="AT2020", N=3, dim=4, quantity="sqrt_sigma", value="485", err_stat="6",
                   units="MeV", table="eqn label eqn_KMeV", note="= 1.160(6)/r0 with r0=0.472(5) fm"))
    out.append(row(source_id="AT2020", N=3, dim=4, quantity="sqrt_sigma", value="440", units="MeV",
                   table="text before eqn_r0MeV",
                   note="approximate traditional value from Regge slope + rotating string model ('something like 440 MeV'); no error quoted"))
    return out


def ltw2004(src):
    tex = read(os.path.join(src, "hep-lat_0404008", "gksSUNgtcwio.tex"))
    out = []
    tab = "Table 14 (label table_continuumG)"
    cols = [("0++", 0), ("0++", 1), ("2++", 0)]
    for r in rows(table_block(tex, "table_continuumG"))[2:]:
        g = clean(r[0])
        m = re.match(r"^SU\((\d+|" + re.escape(BS) + r"infty)\)(_an)?$", g)
        if not m:
            continue
        N = "inf" if "infty" in m.group(1) else m.group(1)
        for (jpc, lev), cell in zip(cols, r[1:]):
            v, e, rest, _ = parse_ve(cell)
            notes = []
            if m.group(2):
                notes.append("anisotropic-lattice SU(4) calculation (separate row in source)")
            if "star" in rest or "*" in rest:
                notes.append("star: poor fit (per caption)")
            if N == "inf":
                notes.append("fit m/sqrt_sigma(N)=m(inf)+c/N^2 to continuum values (eqn eqn_Nmk)")
            out.append(row(source_id="LTW2004", N=N, dim=4, JPC=jpc, quantity="m/sqrt_sigma",
                           value=v, err_stat=e[0], units="dimensionless", table=tab, level=lev,
                           row_label=g + " / " + ["0++", "0++*", "2++"][cols.index((jpc, lev))],
                           note="; ".join(notes)))
    return out


def at2021(src):
    tex = read(os.path.join(src, "2106.00364", "paper_v4b.tex"))
    out = []
    for lab, num in (("table_MKJ_N2-5", 34), ("table_MKJ_N6-12", 35)):
        rs = rows(table_block(tex, lab))
        hdr = [clean(c) for c in rs[1]]
        Ns = [re.search(r"SU\((\d+)\)", h).group(1) for h in hdr[1:]]
        for r in rs[2:]:
            jpc, lev, _ = jpc_from(r[0])
            for N, cell in zip(Ns, r[1:]):
                p = parse_ve(cell)
                if not p:
                    continue
                v, e, rest, _ = p
                out.append(row(source_id="AT2021", N=N, dim=4, JPC=jpc, quantity="m/sqrt_sigma",
                               value=v, err_stat=e[0], units="dimensionless",
                               table="Table %d (label %s)" % (num, lab), level=lev,
                               row_label=clean(r[0]) + " / SU(%s)" % N, note=flags(rest)))
    out += at_jgrid(rows(table_block(tex, "table_MK_J_SUN")), "AT2021", "inf",
                    "Table 38 (label table_MK_J_SUN)", "m/sqrt_sigma", "dimensionless",
                    "large-N extrapolation with O(1/N^2) correction; 2++,1-+,1-- gs exclude smallest N (per text)")
    return out


def teper1998(src):
    tex = read(os.path.join(src, "hep-lat_9804008", "papern.tex"))
    out = []
    rs = rows(table_block(tex, "table_mcont"))
    for r in rs[2:]:
        jpc, lev, _ = jpc_from(r[0])
        for N, cell in zip((2, 3, 4, 5), r[1:]):
            p = parse_ve(cell)
            if not p:
                continue
            v, e, _, _ = p
            out.append(row(source_id="Teper1998", N=N, dim=3, JPC=jpc, quantity="m/sqrt_sigma",
                           value=v, err_stat=e[0], units="dimensionless",
                           table="Table 28 (label table_mcont)", level=lev,
                           row_label=clean(r[0]).replace(BS + "ast", "*") + " / SU(%d)" % N))
    for lab, num, q in (("table_masskN", 34, "m/sqrt_sigma"), ("table_massgN", 33, "m/g2N")):
        for r in rows(table_block(tex, lab))[1:]:
            jpc, lev, _ = jpc_from(r[0])
            v, e, _, _ = parse_ve(r[1])
            out.append(row(source_id="Teper1998", N="inf", dim=3, JPC=jpc, quantity=q, value=v,
                           err_stat=e[0], units="dimensionless", table="Table %d (label %s)" % (num, lab),
                           level=lev, row_label=clean(r[0]).replace(BS + "ast", "*"),
                           note="linear fit in 1/N^2 using N>=%s; slope %s; CL %s%%" % (
                               clean(r[3]), clean(r[2]), clean(r[4]))))
    return out


def at2016(src):
    tex = read(os.path.join(src, "1609.03873", "paper.tex"))
    out = []
    for r in rows(table_block(tex, "table_sigf_cont"))[1:]:
        g = clean(r[0])
        m = re.match(r"^SU\((\d+|" + re.escape(BS) + r"infty)\)\^?([ab]?)$", g)
        N = "inf" if "infty" in m.group(1) else m.group(1)
        vol = {"a": "medium volume", "b": "large volume"}.get(m.group(2), "")
        for q, cell in (("sqrt_sigma/g2N", r[1]), ("m/g2N", r[3])):
            v, e, _, _ = parse_ve(cell)
            out.append(row(source_id="AT2016", N=N, dim=3, JPC="0++" if q == "m/g2N" else "",
                           quantity=q, value=v, err_stat=e[0], units="dimensionless",
                           table="Table 18 (label table_sigf_cont)", level=0 if q == "m/g2N" else "",
                           row_label=g, note=vol))
    groups = [("table_mksu2_cont", 20, 2), ("table_mksu3_cont", 21, 3), ("table_mksu4_cont", 22, 4),
              ("table_mksu6_cont", 23, 6), ("table_mksu8_cont", 24, 8), ("table_mksu12_cont", 25, 12),
              ("table_mksu16_cont", 26, 16)]
    for lab, num, N in groups:
        rs = rows(table_block(tex, lab))[2:]
        w = 4 if N == 8 else 3
        cur = [None, None]
        cnt = [0, 0]
        for r in rs:
            for side in (0, 1):
                cells = r[side * w:(side + 1) * w]
                if len(cells) < 2:
                    continue
                if clean(cells[0]):
                    cur[side] = jpc_from(cells[0])[0]
                    cnt[side] = 0
                p = parse_ve(cells[1])
                if not p:
                    continue
                v, e, _, _ = p
                note = "chi2/ndf " + clean(cells[-1])
                prev = [x for x in out if x["N"] == str(N) and x["JPC"] == cur[side]]
                if prev and prev[-1]["value"] == v and prev[-1]["err_stat"] == e[0]:
                    note += "; identical to previous level in source (likely duplicated row)"
                out.append(row(source_id="AT2016", N=N, dim=3, JPC=cur[side], quantity="m/sqrt_sigma",
                               value=v, err_stat=e[0], units="dimensionless",
                               table="Table %d (label %s)" % (num, lab), level=cnt[side],
                               row_label="%s level %d (side %s)" % (cur[side], cnt[side], "LR"[side]),
                               note=note))
                cnt[side] += 1
    for lab, num in (("table_mksuN_J0", 27), ("table_mksuN_J2", 28), ("table_mksuN_J1", 29)):
        cur = None
        cnt = 0
        for r in rows(table_block(tex, lab))[2:]:
            if clean(r[0]):
                cur = jpc_from(r[0])[0]
                cnt = 0
            v, e, rest, _ = parse_ve(r[1])
            nge = clean(r[-2])
            extra = ("fit in 1/N^2 (+1/N^4 if c2 given) using N>=" + nge) if nge else ""
            out.append(row(source_id="AT2016", N="inf", dim=3, JPC=cur, quantity="m/sqrt_sigma",
                           value=v, err_stat=e[0], units="dimensionless",
                           table="Table %d (label %s)" % (num, lab), level=cnt,
                           row_label="%s level %d" % (cur, cnt),
                           note="; ".join(x for x in (flags(rest), extra) if x)))
            cnt += 1
    return out


FILES = {
    "glueball_su3_4d_morningstar_peardon_1999.csv": mp1999,
    "glueball_su3_4d_chen_2006.csv": chen2006,
    "glueball_su3_4d_athenodorou_teper_2020.csv": at2020,
    "glueball_sun_4d_lucini_teper_wenger_2004.csv": ltw2004,
    "glueball_sun_4d_athenodorou_teper_2021.csv": at2021,
    "glueball_sun_3d_teper_1998.csv": teper1998,
    "glueball_sun_3d_athenodorou_teper_2016.csv": at2016,
}


def main():
    src = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else "data"
    os.makedirs(outdir, exist_ok=True)
    scales = []
    for fn, fun in FILES.items():
        rs = fun(src)
        scales += [r for r in rs if r["dim"] == "4" and not r["JPC"]]
        with open(os.path.join(outdir, fn), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=COLS, lineterminator="\n")
            w.writeheader()
            w.writerows(rs)
        print(fn, len(rs))
    # copy of the scale-setting rows (no J^PC) from the 4D files, for convenience
    with open(os.path.join(outdir, "string_tension_scales_su3_4d.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(scales)
    print("string_tension_scales_su3_4d.csv", len(scales))


if __name__ == "__main__":
    main()
