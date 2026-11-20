"""Torus j_max=2 rerun (zero-flux sector, existing code) to test 'converged for g^2 >= 3'."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.ks_su2 import periodic, build, hamiltonian, lowest

t0 = time.time()
basis, cas, W = build(periodic(2, 2), 4)
print(f"dim {len(basis)} built in {time.time()-t0:.0f}s", flush=True)
for g2 in (2.0, 3.0, 5.0, 10.0):
    w = lowest(hamiltonian(cas, W, g2), k=3)
    print(f"g2={g2} E0={w[0]:.10f} E1={w[1]:.10f} Delta={w[1]-w[0]:.10f}", flush=True)
