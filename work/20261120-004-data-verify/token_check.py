"""Machine cross-check: every CSV row's printed token v(e) must occur in the text of the PDF page
that holds its table (stdlib only).

Usage (repo root): python token_check.py <csv> <pages.txt> "<table prefix>=<page>" ...
<pages.txt> is `pdftotext -layout paper.pdf pages.txt` output (pages split on form feed).
Token search ignores whitespace. This only shows each value+error string is printed on the right
page; row/column assignment was checked by eye against the rendered page (see README).
"""
import csv
import re
import sys


def tokens(value, err):
    """Printed forms: last-digit notation 9.47(116), or literal 1.7(1.5) when err >= 1."""
    d = len(value.split(".")[1]) if "." in value else 0
    e = float(err)
    out = [f"{value}({round(e * 10**d)})"]
    if e >= 1 and d > 0:
        out.append(f"{value}({err})")
    return out


csv_path, txt_path, *maps = sys.argv[1:]
pages = open(txt_path, encoding="utf-8", errors="replace").read().split("\f")
pages = [re.sub(r"\s+", "", p) for p in pages]
tmap = {k: int(v) for k, v in (m.split("=") for m in maps)}
n_ok = n_bad = 0
with open(csv_path, newline="", encoding="utf-8") as f:
    for i, r in enumerate(csv.DictReader(f), start=2):
        tab = r["table"].split(" (")[0]
        if tab not in tmap or not r["err_stat"]:
            continue
        ts = tokens(r["value"], r["err_stat"])
        t = " or ".join(ts)
        if any(x in pages[tmap[tab] - 1] for x in ts):
            n_ok += 1
        else:
            n_bad += 1
            print(f"NOT FOUND line {i}: {tab} N={r['N']} {r['JPC']} lvl{r['level']} token {t}")
print(f"{csv_path}: {n_ok} found, {n_bad} not found")
