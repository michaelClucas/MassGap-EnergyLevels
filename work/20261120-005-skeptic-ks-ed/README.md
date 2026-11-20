# 20261120-005: skeptic review of the SU(2) KS exact diagonalization (claims C1, C2, C3)

## Question

Do claims C1-C3 in `CLAIMS.md` (from `work/20261120-003-su2-ks-ed/`, `src/ks_su2.py`) survive an independent attack on normalization, perturbation theory, weak-coupling limits, j_max truncation, and sector selection?

Setting throughout: 2+1D Hamiltonian lattice, SU(2), a = 1, finite volume (open 1-plaquette, open 2-plaquette ladder, periodic 2x2). Every number here is a finite-volume, finite-spacing observation unless marked proved.

## Method

- `plaq_independent.py` (no import from `src/`): single plaquette as a tridiagonal matrix in the character basis, H_jj = 2 g^2 j(j+1), H_{j,j+-1/2} = -2/g^2. Exact rational Rayleigh-Schrodinger PT to 4th order. Separate finite-difference solve of the theta equation. Anharmonic correction to the weak-coupling limit derived by hand. Output `plaq_independent.log`.
- `ladder_independent.py` (no import from `src/`): ladder reduced to the theta graph (outer 3-link paths have one spin each since corner vertices are 2-valent). Basis f_abc(P1, P2) with my own Racah-formula CG coefficients; Tr P matrix elements from the CG series plus Haar orthogonality; checked pointwise against explicit D matrices at random group elements. Compared to every ladder row of `results.csv`, then pushed to j_max = 16. Output `ladder_independent.log`.
- `torus_sectors.py` (uses `src/ks_su2.py` matrix elements, which 003 validated pointwise): reruns the BFS seeded from a fundamental winding (Polyakov) loop instead of the vacuum, giving the lowest state in each nonzero 't Hooft electric-flux sector. Output `torus_sectors.log`.
- `torus_jmax2.py`: zero-flux torus at j_max = 2 with the existing code, to test "converged for g^2 >= 3". Output `run_torus_jmax2.log`.

All runs pinned, 1 thread.

## Results per claim

### Normalization (affects all three): survives

Wilson action S = (2/g^2) sum (N - Re Tr U_p) with beta = 2N/g^2 gives (1/4) F^2 in the naive continuum limit (checked: Tr U ~ N - g^2 a^4 F^a F^a / 4). Its Hamiltonian limit is (g^2/2) sum E^2 - (1/g^2) sum Tr(U_p + U_p^dag) + const, i.e. Hamer's x = 2/g^4 convention, with E^2 the Casimir (j(j+1)) for generators sigma/2. That is exactly what 003 uses. Strong-coupling gap 4 * (3/4) * g^2/2 = 3g^2/2: correct (proved, elementary).

### C2 (strong-coupling series, single plaquette): survives, now independently derived

Exact PT (proved: nondegenerate RS PT on the tridiagonal matrix; odd orders vanish because the perturbation has zero diagonal and the chain is bipartite):

    E0 = -8/(3 g^6) + ...,   E(j=1/2) = 3g^2/2 + 16/(15 g^6) + ...
    Delta = 3g^2/2 + 56/(15 g^6) - 20144/(3375 g^14) + O(g^-22)

Hand derivation of the 2nd-order term agrees: (2/g^2)^2 [2/(3g^2) - 2/(5g^2) + 2/(3g^2)] = 56/(15 g^6). Adding the new g^-14 term moves PT onto ED: g^2 = 3: 4.635542 vs ED 4.635626 (PT2 4.638272); g^2 = 5: 7.5297903 vs 7.5297906. Basis-free finite-difference solve matches the spin basis to <= 2e-9.

Wording nit: "matches ED to the printed digits at g^2 = 10" is weak evidence (next term is 6e-7 there). The g^2 = 3, 5 agreement with the 4th-order term is the stronger check. The "series derivation not independently checked" caveat can be dropped.

### C3 (weak-coupling harmonic limits): survives; plaquette limit upgraded, ladder confirmed numerically

- Single plaquette (proved at the level of standard semiclassics: 1D Schrodinger operator on (0, pi) with Dirichlet BC, harmonic well at theta = 0 with mass 1/g^2 and omega = 2; odd oscillator levels n = 1, 3 give Delta -> 2 omega = 4). First correction, from the quartic term -theta^4/(6 g^2) with <theta^4>_n = (3/4) l^4 (2n^2 + 2n + 1), l^2 = g^2/2:

      Delta = 4 - (5/8) g^2 + O(g^4)

  ED: (4 - Delta)/g^2 = 0.6258, 0.6266, 0.6291 at g^2 = 0.01, 0.02, 0.05, extrapolating to 0.625. So the "approached with O(g^2) corrections" statement now has its coefficient.
- Ladder: harmonic analysis re-derived independently. Link mass 4/g^2 per color component, stiffness (4/g^2) C^T C, so omega^2 = eig(C C^T) = eig([[4, +-1], [+-1, 4]]) = 3, 5. Each normal mode is a color triplet; Gauss law at the root of a gauge-fixing tree leaves global SU(2) conjugation, so physical states are singlets: lowest excitations 2 omega_1 = 2 sqrt3 (P+, either mode squared is even) and omega_1 + omega_2 = sqrt3 + sqrt5 (P-, the two modes have opposite reflection parity). Same as 003. Numerical extrapolation: see ladder section below.
- Is the open single plaquette a sensible stepping stone? As a code test yes: tree gauge fixing leaves one SU(2) loop variable, and the 4 link Casimirs all act as the full Laplacian on it, which is why the kinetic coefficient is 4 g^2/2. As physics no: it has no notion of volume or of noncontractible loops (no torons), and its weak-coupling gap is the lattice-cutoff oscillator 2 omega = 4 in lattice units. 003 already says this; the claim's last clause is accurate.

### Ladder, independent ED (C1, C3)

Independent theta-graph construction (own CG, own D matrices): pointwise identity (Tr P1 + Tr P2) f = sum W f holds to 1.1e-14; W symmetric to 6e-16; [P, W] = 0 exactly. Against all 200 ladder rows of `results.csv` (P+ and P-, every g^2, j_max 1/2..8): basis dimensions identical, max |E - E_ref| = 5e-11. So the 003 ladder numbers, including the recoupling signs, are reproduced by code that shares nothing with `src/ks_su2.py`.

Truncation beyond 003 (gaps from the global E0):

| j_max | g^2=0.2 P+ | g^2=0.2 P- | g^2=0.3 P- | g^2=0.5 P- |
|---|---|---|---|---|
| 8 | 3.40053519 | 3.84704484 | 3.78369296 | 3.65030945 |
| 9 | 3.40053417 | 3.84701684 | 3.78369295 | 3.65030945 |
| 10 | 3.40053412 | 3.84701614 | 3.78369295 | 3.65030945 |

003's g^2 = 0.2 P- value 3.847 is right to ~3e-5 (their own step estimate 7e-4 was conservative).

Weak coupling (j_max = 10): (2 sqrt3 - Delta_P+)/g^2 = 0.312, 0.318 and (sqrt3 + sqrt5 - Delta_P-)/g^2 = 0.580, 0.606 at g^2 = 0.1, 0.2: smooth O(g^2) approach to the harmonic values, consistent with C3 (observation; g^2 = 0.1 not checked at j_max > 10). At g^2 <= 0.05, j_max = 10 is not converged (gaps above the limit), so no smaller couplings are claimed.

### C1: weakened (wording), numbers survive

1. **Sector exclusion on the torus is material.** The torus basis is the zero 't Hooft electric-flux sector. Seeding the same BFS from a fundamental winding loop gives the lowest state in flux (1,0) / (0,1) at

   | g^2 | zero-flux gap (003, j_max 3/2) | E_flux(1,0) - E0, j_max 1 | j_max 3/2 |
   |---|---|---|---|
   | 10 | 14.8606 | 7.4932 | 7.4932 |
   | 5 | 7.2351 | 3.6965 | 3.6964 |
   | 3 | 4.1281 | 2.0378 | 2.0342 |
   | 2 | 2.6309 | 1.0729 | 1.0331 |
   | 1 | 1.1286 | 0.6105 | 0.2961 (not converged) |

   At strong coupling this is a single 2-link winding string, cost 3g^2/4, half the zero-flux 3g^2/2. So the lowest gauge-invariant excitation on the 2x2 torus is a torelon (flux-tube) state, not the zero-flux E1. Excluding flux sectors is standard (they are superselected by the center symmetry, and their energy grows like sigma L, so they are not glueballs), but then the claim must say "zero-flux sector". The flux-(1,1) sector sits just below the zero-flux gap (14.657 vs 14.861 at g^2 = 10).
   A j = 1 winding loop (center-trivial) is in the vacuum-connected basis, consistent with 003's identification of that basis as the full zero-flux sector (not proved).
2. **Momentum.** The torus E1 is the lowest zero-flux state over all momenta. At strong coupling it is a member of the one-plaquette multiplet, not necessarily a p = 0 glueball. Already in 003 open issue 2; CLAIMS.md should carry it.
3. **"Gap > 0" is not a numerical finding.** On any finite lattice the KS Hamiltonian is a Schrodinger operator on a compact manifold (discrete spectrum) with a positivity-improving heat semigroup, so the ground state is nondegenerate and Delta > 0 at every g (cited in 003: Reed-Simon IV XIII.12). This holds per flux sector too. The content of C1 is the values, not the sign. Status "numerical" is fine for the values.
4. **Torus convergence.** Survives, and is better than stated. Zero-flux torus at j_max = 2 (dim 90899, existing code, `run_torus_jmax2.log`): Delta = 2.625455 (g^2 = 2), 4.128011 (3), 7.235141 (5), 14.860551 (10), against 003's j_max = 3/2 values 2.630887, 4.128073, 7.235141, 14.860551. So the change is 6e-5 at g^2 = 3 and 5e-3 at g^2 = 2: "converged for g^2 >= 3" holds to ~1e-4, and g^2 = 2 is good to ~1e-2.
5. Ladder P- at g^2 = 0.2: 003's own last-step change is 7e-4, so "3.847" is good to about 1e-3. See ladder section for the j_max = 16 value.

## Verdict

- Normalization: survives.
- C1: numbers survive; wording weakened. Torus value is a zero-flux-sector, all-momentum gap; the true lowest gauge-invariant excitation on the 2x2 torus is a torelon at about half that energy for g^2 >= 2.
- C2: survives, independently derived, next term -20144/(3375 g^14) added.
- C3: survives. Single-plaquette limit 4 with correction -(5/8) g^2 derived; ladder limits re-derived and checked numerically.

Proposed CLAIMS.md wording:

- C1: "Gap Delta(g^2) > 0 (guaranteed at finite volume; the content is the values) ... 2x2 torus: zero 't Hooft electric-flux sector only, all momenta mixed, j_max = 3/2 converged to ~1e-4 for g^2 >= 3 (checked at j_max = 2). Nonzero-flux (torelon) states lie lower, at about half the zero-flux gap at strong coupling (`work/20261120-005-skeptic-ks-ed/`)."
- C2: "Delta = 3g^2/2 + 56/(15 g^6) - 20144/(3375 g^14) + O(g^-22); derivation independently checked (005); status proved (PT series) plus numerical agreement."
- C3: "single plaquette: Delta = 4 - (5/8) g^2 + O(g^4) (leading term and first correction derived, standard semiclassics); ladder harmonic limits re-derived independently, numerically confirmed (005)."

## Open issues

1. Rigor of the harmonic limit: the plaquette case is a textbook semiclassical statement (Simon 1983, semiclassical analysis of low-lying eigenvalues, cited from memory, unchecked). The ladder case relies on the singlet projection argument at leading order; no error bound.
2. Whether the vacuum-connected BFS basis equals the full zero-flux sector is checked only for the j = 1 winding loop, not proved.
3. Torus j_max = 2 took 21 min and ~5 GB (pure Python BFS). Anything larger needs momentum blocking.
4. Torus flux sectors were run only to j_max = 3/2; values for g^2 <= 2 there are not converged.
