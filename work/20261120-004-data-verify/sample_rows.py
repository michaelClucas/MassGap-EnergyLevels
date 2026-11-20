"""Draw the random audit sample of data/ rows (stdlib only).

Pools all rows of the seven glueball CSVs (string_tension_scales_su3_4d.csv is a copy of rows
already in them), excludes the 0++/2++/0-+ ground-state and scale rows checked exhaustively,
then draws 30 rows with random.seed(0). Run from the repo root.
"""
import csv
import glob
import random

HEADLINE = {"0++", "2++", "0-+"}

rows = []
for path in sorted(glob.glob("data/glueball_*.csv")):
    with open(path, newline="", encoding="utf-8") as f:
        for i, r in enumerate(csv.DictReader(f), start=2):
            r["_file"], r["_line"] = path.split("/")[-1].split("\\")[-1], i
            rows.append(r)

def is_headline(r):
    return (r["JPC"] in HEADLINE and r["level"] == "0") or r["JPC"] == ""

pool = [r for r in rows if not is_headline(r)]
random.seed(0)
sample = random.sample(pool, 30)
print(f"total rows {len(rows)}, non-headline pool {len(pool)}")
for r in sorted(sample, key=lambda r: (r["_file"], r["_line"])):
    print(f'{r["_file"]}:{r["_line"]} | N={r["N"]} {r["JPC"]} lvl{r["level"]} {r["quantity"]} '
          f'= {r["value"]}({r["err_stat"]}){"(" + r["err_sys"] + ")" if r["err_sys"] else ""} | {r["table"].split(" (")[0]}')
