"""Self-test of src/inertia.py against dense eigenvalues: random symmetric (dense and sparse paths),
near-degenerate shifts, and a sparse 2D Laplacian. Usage (repo root):
python tools/pin.py .venv/Scripts/python.exe work/20261120-024-inertia-certification/selftest.py"""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src import inertia as I  # noqa: E402

rng = np.random.default_rng(0)
bad = 0
for n in (50, 400, 1500):
    B = rng.standard_normal((n, n)); A = (B + B.T) / 2
    w = np.linalg.eigvalsh(A)
    for s in np.quantile(w, [0.01, 0.3, 0.5, 0.9]) + 1e-3:
        true = int(np.sum(w < s))
        for dm in ((10**9, 0) if n <= 400 else (10**9,)):  # dense path; sparse path only n <= 400 (dense random fill)
            r = I.inertia_below(sp.csr_matrix(A), s, dense_max=dm)
            ok = r["n_lo"] == true == r["n_hi"]
            bad += not ok
            print(f"n={n} s={s:+.4f} path={'dense' if dm else 'sparse'} true={true} "
                  f"bracket=[{r['n_lo']},{r['n_hi']}] r={max(r['r_lo'], r['r_hi']):.1e} {'ok' if ok else 'MISMATCH'}")
# shift within 1e-15 of an eigenvalue: must be inconclusive or correct, never wrong
B = rng.standard_normal((300, 300)); A = (B + B.T) / 2; w = np.linalg.eigvalsh(A)
r = I.inertia_below(A, w[10]); print("shift at eigenvalue 10:", r["n_lo"], r["n_hi"], "(true count below is 10 or 11)")
bad += not (r["n_lo"] in (10, 11, None) and r["n_hi"] in (10, 11, None) and (r["n_lo"] or 0) <= (r["n_hi"] or 99))
# sparse 2D Laplacian 60x60 (n = 3600), known spectrum
m = 60
T = sp.diags([-np.ones(m - 1), 2 * np.ones(m), -np.ones(m - 1)], [-1, 0, 1])
Lp = (sp.kron(T, sp.identity(m)) + sp.kron(sp.identity(m), T)).tocsc()
ev1 = 2 - 2 * np.cos(np.pi * np.arange(1, m + 1) / (m + 1))
w = np.sort((ev1[:, None] + ev1[None, :]).ravel())
for s in (w[0] + 1e-4, (w[5] + w[6]) / 2, 1.0):
    r = I.inertia_below(Lp, s, dense_max=0)
    true = int(np.sum(w < s)); ok = r["n_lo"] == true == r["n_hi"]; bad += not ok
    print(f"laplacian s={s:.5f} true={true} bracket=[{r['n_lo']},{r['n_hi']}] r={r['r_hi']:.1e} {r['method']} {'ok' if ok else 'MISMATCH'}")
print("SELFTEST", "PASS" if bad == 0 else f"FAIL ({bad})")
