# State

Iteration: 1 (2026-11-20)

## Status
Iteration 1 done: literature survey, lattice glueball data, Python env, first SU(2) Kogut-Susskind exact diagonalization. No proofs yet. Lean toolchain not installed (needs user OK to install `elan`).

Compute rule changed during iteration 1: this program owns logical CPUs 9-11 only; all compute runs via `python tools/pin.py <cmd>`, and `tools/repin.ps1` runs at the start of each iteration and after subagents return. Repin found nothing to fix in iteration 1.

## Settled
- Data: `data/` holds continuum-extrapolated glueball spectra from 7 papers (4D SU(3), 4D SU(N) N=2..12 and N->inf, 2+1D SU(N)), extracted by script from arXiv LaTeX sources (`work/20261120-002-glueball-data/`). Not yet independently verified.
  - Headline (Athenodorou-Teper): 4D SU(3) m(2++)/m(0++) = 1.437(11), m(0-+)/m(0++) = 1.549(16), m(0++)/sqrt(sigma) = 3.405(21); 4D SU(inf) 1.497(8), 1.534(11), 3.072(14); 2+1D SU(3) 1.658(5), 2.114(10), 4.368(7).
  - Known discrepancies between sources (SU(3) m(0++)/sqrt(sigma) AT2020 vs MP1999/Chen ~2 sigma; LTW2004 vs AT2020 on 2++/0++ and SU(8)). See work README.
- Prior art (`literature/`, `work/20261120-001-literature-survey/`), so we don't claim novelty for:
  - KS glueballs by coupled cluster, variational, ED of SU(2)/SU(3) plaquette chains, qubit plaquette chain with continuum limit (2512.11068).
  - Martingale gap proof for frustration-free deformed 3+1D Z2 gauge theory (2409.10612). No non-Abelian version found.
  - Bakry-Emery lattice gap at strong coupling (Shen-Zhu-Zhu 2204.12737); fails for U(N) flat directions (2510.22788). Continuum orbit-space version only formal.
  - theta as Bloch momentum (Jackiw-Rebbi, Callan-Dashen-Gross 1976); Koller-van Baal finite-volume spectra.
  - Gluon BCS vacuum models; Cornwall dynamical mass; <A^2> gap equation.
  - Rigorous lattice gap known only at strong coupling (Osterwalder-Seiler, Shen-Zhu-Zhu, Cao-Nissim-Sheffield 2505.16585, 2509.04688).
- SU(2) KS ED (`work/20261120-003-su2-ks-ed/`, `src/ks_su2.py`): gauge-invariant spin-network basis validated to 1e-12..1e-15 (pointwise plaquette action, gauge invariance, group-manifold cross-check 1.2e-8). Single plaquette matches strong-coupling PT Delta = 3g^2/2 + 56/(15 g^6). Weak-coupling gap on tiny open lattices tends to a g-independent harmonic value (4 for one plaquette), set by geometry, not continuum physics.

## Live lines of attack (ranked)
1. Line 2, spectral-gap methods: frustration-free non-Abelian gauge model (SU(2)_k string-net or q-deformed KS), gap via Knabe / Lemm-Mozgunov finite-size criterion with certified ED, then stability toward KS. No prior non-Abelian work found (web only; needs INSPIRE citation check of 2409.10612).
2. Line 1, certified two-sided gap bounds (Temple/Weinstein + truncation error) for small SU(2) KS lattices. Feeds line 2.
3. Line 3, Payne-Weinberger lower bound on the finite-lattice Gribov region (convex; log det FP concave). Expected to vanish in infinite volume: candidate finite-volume result or Result-C obstruction.
4. Line 6, O(3) quantum rotor chain as test bed for every method.
5. Line 5, interaction-driven: low priority, mostly published.

## Dead / Result-C candidates
- Line 4 (theta/Bloch): likely obstruction. theta-band structure describes vacuum energy across theta sectors, not the glueball gap within one sector. Needs an explicit write-up to count.

## Open questions
- Does Yarotsky (math-ph/0412040) already give the strong-coupling KS gap including Gauss law?
- 2x2 torus ED needs j_max >= 2 (~1e5 states) and momentum blocking to reach g^2 ~ 1. Toron behaviour (Delta ~ g^2 at weak coupling on the torus) untested.
- Lüscher 1983 / Koller-van Baal references in the ED README were cited from memory; verify.
- Unvetted claimed full proofs (2506.00284, 2611.11534, a Zenodo record) listed in literature; not relied on.

## Next planned tasks (iteration 2)
- Data verifier: independent second agent re-extracts a sample of `data/` rows from the papers and checks SOURCES.md.
- Mathematical physicist + literature: design a frustration-free non-Abelian lattice gauge Hamiltonian and state the Knabe/Lemm-Mozgunov finite-size condition it needs; INSPIRE citation check of 2409.10612, 2204.12737, Dell'Antonio-Zwanziger 1991; check Yarotsky vs Gauss law.
- Computational physicist: certified Temple/Weinstein two-sided gap bounds with truncation-error control for plaquette, ladder, 2x2 torus; momentum blocking to push the torus to j_max = 2.
- Mathematical physicist: write the line-4 obstruction precisely; scope the Payne-Weinberger/Gribov bound.
- Lean: blocked on user approval to install `elan`.
