# 20261120-012: 2x2 torus, all electric-flux sectors, momentum and point-group blocking, certified bounds

## Question

On the SU(2) Kogut-Susskind 2x2 periodic torus (2+1D), what is the honest lowest gap once all Z2 electric ('t Hooft) flux sectors are included, what is the gap inside the zero-flux sector, how do the low levels split by momentum and C4v irrep, and how far toward weak coupling can certified two-sided bounds be pushed?

Setting of every number below: 2+1D Hamiltonian lattice, gauge group SU(2), a = 1, periodic 2x2 lattice (8 links, 4 plaquettes), finite volume, fixed g^2, H = (g^2/2) sum_l j_l(j_l+1) - (2/g^2) sum_p Tr U_p (same H as `src/ks_su2.py`). Nothing here is a continuum or infinite-volume statement. Finite-truncation ED values are observations; certified intervals are proofs about this finite lattice up to floating-point rounding (not interval arithmetic).

## Method

**Basis and sectors.** `src/torus_sym.py` (new module; `src/ks_su2.py`, `src/certify.py` unchanged). The flux sector (fx, fy) is the parity of sum 2j over the x-links at x = 0 (resp. y-links at y = 0), conserved by every Tr U_p. Each sector basis is enumerated completely up to spin j_max with `certify.sector_states` (all spin configurations admitting an invariant at each vertex, all intertwiner labels), so P is the full truncated sector, as the Feshbach Lemma requires. The skeptic's winding-loop BFS (005) is not needed: enumeration reaches every state, and `build_sector` closure was already checked in 006.

**Translations (momentum).** `torus2x2()` reorders the four legs at every vertex by direction (E-out, N-out, W-in, S-in). Since `ks.intertwiners` is a deterministic function of the leg tuple, every vertex then uses the same intertwiner basis for the same spins, and translations act on spin networks as pure permutations, no phases. Momentum states are the usual orbit sums over Z2 x Z2 with characters (-1)^(kx gx + ky gy), k in {0, pi}^2; H_k is built by applying the plaquette moves to orbit representatives only (4x fewer than the full BFS), in 2 processes.

**Point group.** `apply_op` implements any lattice automorphism, including link reversal via D^j_{mn}(U^-1) = (-1)^{n-m} D^j_{-n,-m}(U), and overlaps the transformed vertex tensor with the image vertex's intertwiner basis. sigma_x (x -> -x) and sigma_d (x <-> y) are used to label levels at k = (0,0) and (pi,pi) in the zero-flux sector by C4v irrep (A1: ++, A2: --, B1: sigma_x +, sigma_d -, B2: sigma_x -, sigma_d +; E: 2-dim, traces 0). The point group is used for labelling only, not for further blocking.

**Certified bounds.** The Schur-complement Lemma of `work/20261120-006-*/README.md`, applied per (flux, k) block (P, Q and the shell are translation invariant, so the block decomposition is exact). Zero flux: P = spins <= 3/2, exact shell X from the complete spin-2 states. Flux sectors: P = spins <= 1, exact shell from spin 3/2 (the j_max = 2 flux build was not run, see Compute). Upper bounds: min-max, using the j_max = 2 (zero flux) or 3/2 (flux) block Ritz values.

New: an improved floor on Q H Q (proved here, `split_floor`). On the 2x2 torus every link lies in exactly two plaquettes, so for t in [0, 1]

    H = (1-t)(g^2/2) sum_l C_l + sum_p h_p(t),   h_p(t) = t (g^2/4) sum_{l in p} C_l - (2/g^2) Tr U_p.

On Q (gauge-invariant states with some spin > j_max) the first term is >= (1-t)(g^2/2) c_Q. Each h_p(t) acts on L^2(SU(2)^4) of its four links; its heat semigroup is positivity improving, so its ground state is unique and positive, hence invariant under gauge transformations at the four corners, hence a class function of U_p, on which sum_{l in p} C_l acts as 4 j(j+1). So inf spec h_p(t) = lowest eigenvalue of the tridiagonal character-basis matrix (diagonal t g^2 j(j+1), off-diagonal -2/g^2), bounded from below by the same Lemma with a truncation at 2j = 400. Then d = max_t [(1-t)(g^2/2) c_Q + 4 e_p(t)] >= crude floor (t = 0). Gain: +0.9 at g^2 = 3, +1.3 at g^2 = 2, +1.8 at g^2 = 1.5 for c_Q = 11 (`validate.log`). (proved, using the Perron-Frobenius theorem for positivity-improving semigroups, cited from memory: Reed-Simon IV, Sec. XIII.12.)

Choice of E: when U < d, the standard L = f(U) of 006; when U >= d, bisection for the largest E < min(d, U) with lambda_k(A - X/(d-E)) >= E (same Lemma, any such E is a valid lower bound). `refine.py` applies the bisection also to the zero-flux k = 0 block at g^2 = 1.5, where the standard choice is poor because U is close to d.

## Validation (`validate.log`, `check_01.log`)

- Direction-ordered torus vs `ks.periodic(2,2)`, zero flux, j_max = 1: spectra equal to 4e-14.
- Tx, Ty, sigma_x, sigma_d, rot90: orthogonal and commute with W to 4e-15 in every flux sector they preserve; Tx is a pure permutation.
- Union of momentum-block spectra = full-sector spectrum to 6e-14, all four flux sectors (j_max = 1).
- (0,1) spectrum = (1,0) spectrum (blocks k <-> (ky, kx)), j_max = 3/2, to 1.4e-13. So (0,1) is not rebuilt at j_max = 2.
- Zero-flux k = 0 gaps at j_max = 2 reproduce the skeptic's unblocked j_max = 2 run (005) to all printed digits: 2.625455, 4.128011, 7.235141, 14.860551 at g^2 = 2, 3, 5, 10.

## Results

### Spectrum by sector and momentum (ED, observation)

E0 = global ground state (zero flux, k = 0, A1, at every g^2). Entries are E - E0 for the lowest level of each block (second level for the ground block). Zero flux at j_max = 2; flux sectors at j_max = 3/2 (so at g^2 <= 1.5 the flux - zero-flux differences mix truncations; the consistent j_max = 3/2 values are in brackets).

| g^2 | (0,0) k=0 E1 [irrep] | next (0,0) k=0 A1 | (0,0) k=(pi,0) | (0,0) k=(pi,pi) | (1,0) k=0 torelon | (1,0) k=(0,pi) | (1,1) k=0 |
|---|---|---|---|---|---|---|---|
| 1 | 0.883 [B1] (j 3/2: 1.129) | 1.091 | 2.662 | 3.426 [B2] | 0.436 (0.296) | 2.421 | 0.796 (0.655) |
| 1.5 | 1.658 [B1] | 1.841 | 2.992 | 3.242 [B2] | 0.560 (0.555) | 1.945 | 1.011 |
| 2 | 2.6255 [A1] | 2.6255 | 3.443 | 3.483 | 1.0332 | 1.928 | 1.729 |
| 3 | 4.1280 [A1] | 4.1280 | 4.640 | 4.650 | 2.0342 | 2.388 | 3.471 |
| 5 | 7.2351 [A1] | 7.2351 | 7.531 | 7.533 | 3.6964 | 3.780 | 6.832 |
| 10 | 14.8606 [A1] | 14.8606 | 15.004 | 15.004 | 7.4932 | 7.504 | 14.657 |

- **Overall gap** Delta_all (lowest level above E0 in any sector): always the (1,0)/(0,1) torelon at k = 0, about half the zero-flux gap at g^2 >= 2 (strong-coupling value 3g^2/4 vs 3g^2/2), falling to ~0.3-0.45 at g^2 = 1 (not converged there).
- **Zero-flux gap** Delta_00: at k = 0 for every g^2; nonzero-momentum zero-flux states lie higher at all couplings sampled. The lowest zero-flux excitation is A1 for g^2 >= 2 and switches to B1 at g^2 <= 1.5 (A1 then second, 0.18-0.21 higher). Full level lists with irreps: `results.csv` (`ed-level*` rows).
- Convergence: zero flux j_max 3/2 -> 2 changes Delta_00 by 6e-5 at g^2 = 3, 5e-3 at g^2 = 2, 0.05 at 1.5, 0.25 at 1 (not converged). Flux sectors only to j_max = 3/2; 005 shows j_max 1 -> 3/2 changes of 4e-3 (g^2 = 3), 0.04 (2), 0.31 (1).

### Certified bounds (proved for the untruncated H on this lattice, up to rounding)

| g^2 | Delta_all (all four flux sectors, all k) | Delta_00 (zero flux, all k) | Delta_00 at k = 0 | E0 |
|---|---|---|---|---|
| 1 | none | none | none | upper only, -6.999978 |
| 1.5 | none | [0.7694, 1.7016] | [0.7694, 1.7016] | [-2.864552, -2.820871] |
| 2 | none (lower -0.80) | [2.60552, 2.62583] | same | [-1.292267, -1.291894] |
| 3 | [2.02927, 2.03423] | [4.12793, 4.12801] | same | [-0.3931398, -0.3931392] |

Previously (006): zero-flux torus certified only at g^2 >= 3 ([4.1257, 4.1281] at g^2 = 3, from an unblocked j_max = 3/2 LDL run). Now: g^2 = 2 and 1.5 are certified in the zero-flux sector, the g^2 = 3 interval is 30x narrower, and the overall gap including torelons is certified at g^2 = 3. The g^2 = 1.5 row uses `refine.csv` (E1 >= -2.051487, E0 >= -2.864552) with the k != 0 rows of `results.csv`.

Where and why no certified lower bound is obtained (all "floor d below E1"-type failures):
- g^2 = 1, every sector: d = -7.86 (zero flux, c_Q = 11) and -8.70 (flux, c_Q = 7.5) lie below the ground-state upper bound -7.000. Any E < d is then below E0, so no gap information is possible; rows marked `none(d<=E0_hi)`.
- g^2 = 1.5, flux sectors: d = -2.68 (P at j_max = 1, c_Q = 7.5) is above E0 but below the torelon level (-2.26), and the bisection only gives E0(1,0) >= -4.20, below E0. So Delta_all has no lower bound, although Delta_00 does.
- g^2 = 2, flux sectors: d = 1.27 clears the torelon level, but the standard choice with d - U = 0.64 to 1.5 and the large shell coupling at j_max = 1 gives E0(1,0, k=(0,pi)) >= -2.09 and E0(1,1, k=0) >= -1.64, below E0 = -1.29. A j_max = 3/2 P in the flux sectors (needs the j_max = 2 flux build) or bisection there would very likely fix this; not done (time).

### Comparison with the lattice glueball picture (qualitative only)

In the zero-flux sector at k = 0 the lowest excitation is A1 (the lattice 0+ channel) for g^2 >= 2, as for the lightest glueball in 2+1D SU(2) lattice studies (0++ lightest, 2++ heavier). For g^2 <= 1.5 a B1 level (one of the lattice 2+ channels, square-lattice B1 being a component of J = 2) drops below the A1. On a 2x2 torus a natural candidate is a two-torelon state (two winding fundamental loops whose fluxes cancel), combined antisymmetrically between x and y: its energy 0.883 at g^2 = 1 is close to twice the j_max-2-zero-flux / j_max-3/2-flux torelon difference (2 x 0.436 = 0.87). This is a conjecture from energies alone; no overlap with an explicit torelon-pair state was computed. Finite-volume torelon and ditorelon states contaminating the 0++ and 2++ channels at small L is the known small-volume behaviour in Euclidean lattice glueball work. None of this says anything about continuum glueball masses: the lattice is 2x2 (L ~ one correlation length or less), and the B1 < A1 ordering is presumably a small-volume effect.

## Verdict

- Proved (modulo floating-point rounding, finite 2x2 lattice, 2+1D SU(2) KS, fixed g^2): the zero-flux gap is certified at g^2 = 1.5, 2, 3 and the full gap including all four flux sectors at g^2 = 3 (table above). The overall gap is the (1,0)/(0,1) torelon at k = 0, ~0.49x the zero-flux gap at g^2 = 3.
- Observation (finite truncation): the global first excited state is a torelon at every g^2 sampled (1 to 10); within zero flux it is k = 0, A1 for g^2 >= 2 and B1 for g^2 <= 1.5.
- No certified bound at g^2 = 1 in any sector: the electric floor with c_Q = 11 lies below E0. Getting there needs a larger c_Q (P at j_max = 2: c_Q = 15.5 gives split floor -6.46, still below E1 ~ -6.1) or a nested/iterated floor; with dense blocks only it is out of reach.
- New tooling: `src/torus_sym.py` (direction-ordered torus basis, permutation translations, general automorphisms with link reversal, momentum blocks, split floor).

## Open issues

1. Flux sectors at j_max = 2 not built (each ~13 min and ~5 GB peak with the current tuple-heavy builder); this is what blocks Delta_all at g^2 = 2 and convergence of the torelon at g^2 <= 1.5.
2. Builder memory: the zero-flux j_max = 2 build peaked at ~5.4 GB (main process holding all move tuples), briefly above the 5 GB budget. Streaming the moves into arrays per chunk would fix it.
3. Point group used for labels only. Blocking by C4v at k = 0 would cut the j_max = 2 block from 23089 to ~3000 and make dense certification at P = j_max 2 possible (needs projectors with the non-permutation vertex maps).
4. Ditorelon interpretation of the B1 level is conjectural.
5. Still floating-point, not interval arithmetic (as in 006).
6. `src/reproduce.py` not extended: the j_max = 2 build takes ~14 min; `validate.py` (~1 min) and `check_01.py` are quick and could be added.

## Files and reproduce

- `validate.py` (~1 min): checks above -> `validate.log`.
- `build.py 3 0,0 1,0 0,1 1,1` (~1.5 min) and `build.py 4 0,0` (~14 min, ~5.4 GB peak): blocks to `artifacts/012-torus-sectors/` (gitignored). `run.log` is the j_max = 2 build log (the (1,0) j_max = 2 build was started and stopped).
- `check_01.py` -> `check_01.log`.
- `analyze.py 4` (~15 min): ED by block with irreps and certified bounds -> `results.csv`, `analyze.log`. With only the j_max = 3/2 builds, `analyze.py 3` gives the j_max = 3/2 ED numbers (`analyze_tj3.log` holds the ED part of such a run from an earlier code version; its partial certified lines are superseded).
- `refine.py 1.5` -> `refine.csv`, `refine.log`.

`results.csv` columns: kind (ed, ed-levelN, gap-ed, certified, gap-certified), sector (flux), momentum, irrep (zero flux, k = (0,0)/(pi,pi); "deg2(sx~0, sd~0)" is the 2-dim E irrep with traces at eigensolver precision), g2, j_max (for certified rows: j_max of P), dim, E0, E1, E0_lo/E0_hi/E1_lo/E1_hi, d_floor, method. Certified gap rows for g^2 = 1.5 in `results.csv` predate `refine.py`; the table above uses the refined values.
