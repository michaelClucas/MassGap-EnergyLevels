"""Shared CSV writer for the inertia checks of work 024."""
import csv

import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent
FIELDS = ["claim", "work", "system", "setting", "g", "j_or_l_max", "level_k", "quantity",
          "claimed_lower", "dim", "lanczos_used", "matrix", "n_below_lo", "n_below_hi", "allowed",
          "delta", "residual", "method", "verdict", "note"]


def write(name, rows):
    with open(OUT / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in r.items()})


def res_fields(res):
    return dict(n_below_lo=res["n_lo"], n_below_hi=res["n_hi"], delta=res["delta"],
                residual=float(np.nanmax([res["r_lo"], res["r_hi"]])), method=res["method"], verdict=res["verdict"])
