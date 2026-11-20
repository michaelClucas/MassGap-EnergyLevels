# 20261120-030: revision of the Result-C draft after skeptics 028 and 029

## Question
Apply every required change from the technical skeptic (`028`) and the literature skeptic (`029`) to `work/20261120-025-resultC-draft/claim.md`, keep every number traceable in `number_map.csv`, and stop `src/reproduce.py` runs from overlapping (028 M5).

## Method
Text edits only. No computation, no reproduce run (the full run 027 was executing). Replacement text from the skeptic reports was used verbatim where given. Files changed: `work/20261120-025-resultC-draft/claim.md`, `work/20261120-025-resultC-draft/number_map.csv` (N12, N13, N14 value and note columns only; match strings unchanged), `src/reproduce.py` (lock). New: this README and `new_references.md`.

## Result: issue by issue

### Technical skeptic (028)

| issue | resolution | what changed |
|---|---|---|
| H1 headline overstates blocks | fixed | New "block type" column (proved / shown for our version only / argued / open) with the skeptic's per-row reading; Answer reworded to "each failure is shown ... for the route as the program instantiated it; rows marked argued or open name what would have to change", headline now "obstruction, partial", names row 3 as open and the SDP/bootstrap family as not examined. |
| H2 U(1) misstated | fixed | Row 2 now says the U(1) analogue of the rate c min(1, 1/beta) is false on free-boundary boxes [N19]; volume-uniform positivity neither shown nor refuted. |
| H3 Feynman / KKN missing | fixed | Added to "did NOT examine" with refs (Nucl. Phys. B 188, 479; hep-th/9705087) and m = e^2 N/(2 pi). |
| H3 Moncrief/Mondal not in "open" | fixed | Added to "What remains open"; row 3 marked open-conditional; row 3(a) notes a strongly log-concave J psi_0^2 gives a diameter-free bound. |
| M1 Omega vs Lambda | fixed | Row 3(d) restricted to continuum small T^3 (diam(Lambda) >= 2 pi sqrt(3L), Delta_PW <= g^2/(24L) [N29]); N26/N30 labelled as Omega quantities, not upper bounds on the Lambda constant. |
| M2 j_max ~ g^-2 | fixed | Replaced with j_max >~ sqrt(N_p)/g from the norm-bound floor, marked heuristic, dim ~ j_max^{#links}; also listed under Known weaknesses. |
| M3 inward rounding | fixed | [2.02927, 2.03424], [4.12793, 4.12802], >= 0.7693 in claim.md and number_map.csv (same source rows; notes record the outward rounding and the old value). |
| M4 1/N | fixed | Row 5 adds 1/N as a free parameter and why the matrix-model / planar argument covers it (014 sec. 6, 018). |
| M5 concurrent reproduce | fixed in code, damage not repaired | `src/reproduce.py` now takes `artifacts/.reproduce.lock` (O_CREAT/O_EXCL, PID and args inside), removes it in a `finally`, exits 2 with a message naming the PID if the lock exists; `--check-only` also locks, `--list` does not. Not run or tested here. The running 027 process loaded the old module and holds no lock, so the lock protects only runs started after it finishes. The 027 sandbox damage from 028's quick run still needs the repair 028 gives (`--only` 002 004 005 007 008, then `--check-only`) after 027 ends. |
| L1 Froehlich-Spencer threshold | fixed | "FS prove m = 0 for all sufficiently large beta (non-explicit threshold); beta_XY = 2 lies above the numerically known KT point." |
| L2 series radius | fixed | "Finite radius" removed; row 7 says the radius was not computed (5 coefficients) and a truncated series cannot reproduce e^{-2 pi/g}. |
| L3 Hessian constant | fixed | "attained in d = 2 [N17] (d = 3 open)". |
| L4 beta mismatch N20 vs N21 | partially | Row 2 now states the two numbers are at different beta and that a common-beta pair was not computed. A same-beta number needs a run (`016/check_c8.py`), not done under the no-compute brief. |
| L5 z = 2 argued | fixed | Row 2(d) marked argued, with Chevyrev-Qu-Shen 2608.27828 as consistent evidence. |
| L6 CLAIMS.md C5 stale on Lean | not fixed (out of scope) | Overseer owns CLAIMS.md: `schur_count` / `schur_count_sorted` exist in `lean/MassGap/MinMax.lean`. |
| L7 KS as tight-binding; Gorantla-Huang, Lemm-Mozgunov | fixed | One sentence under setting abbreviations; Gorantla-Huang 2409.10612 and Lemm-Mozgunov 1801.08915 in the literature check. |
| L8 GM prefactor unchecked | not fixed (needs the paper) | Listed under Known weaknesses. |
| 011 not cited in row 2 | fixed | `011` added to row 2. |

### Literature skeptic (029)

| issue | resolution | what changed |
|---|---|---|
| 1 Row 6 contradicted (RDM/SDP/bootstrap) | fixed | Row 6 uses 029's replacement wording, extended: v2RDM gives E_0 lower bounds; the lattice-YM positivity bootstrap bounds observables (Wilson loops), not the gap; gap bootstraps exist for the TFIM and qubit Hamiltonians; none bounds a YM gap. Family added to "did NOT examine" and to "What remains open" as the leading open direction. Row 2(a) narrowed: non-frustration-freeness blocks Knabe/GM/martingale, not every gap-certification method. Row 6 block type: proved (wavefunction methods) + open. |
| 2 Novelty of certified intervals | fixed | "no prior certified intervals" replaced by 029's text naming arXiv:2611.11534 and arXiv:2503.15539 as unrefereed overlapping preprints, overlap in numbers unconfirmed. Row 1(a) cites the 2611.11534 hard-cutoff obstruction. Comparison with 2611.11534 added to "What remains open". |
| 3 Unvetted-proofs list | fixed in claim.md | Jacobsen 2506.00284 and Faizal-Shabir 2606.19362 (Fortschr. Phys. 74 (2026) e70097) listed as claimed full proofs, not relied on; 2611.11534 and 2503.15539 listed as partial. STATE.md needs the same change (overseer). |
| 4 AKLT published bound stronger | fixed | Row 2 says our GM numbers recompute the curve in 2411.03680 Fig. 2 and quotes Rai et al. >= 0.34976 (Knabe 0.24806) as stronger; Known weaknesses notes eps_12 is superseded. Rozman-Molnar-Schuch in the literature check. |
| 5 RK attribution | fixed | Henley cond-mat/0311345 and Castelnovo-Chamon-Mudry-Pujol 2005 credited; Fendley-Fradkin kept as a remark; 2608.27828 added to row 2(d). |
| 6 HKS 1979 O(3) series | partially | Cited as prior art in row 7 and the literature check. Whether our five coefficients agree with it was not checked (paper not read); stated as such. |
| 6 FS wording | fixed | As L1 above. |
| 6 Chatterjee 2602.00436 | fixed | Row 5 and literature check. `literature/02` update left to overseer. |
| 6 "from memory" citations | fixed | Known weaknesses narrowed to: bibliographic data resolved, content not re-read; Maehly unresolved. |
| 6 `literature/02` journal refs | not fixed (out of scope) | Listed in `new_references.md` for the overseer. |
| 6 FBSE journal ref | fixed | Row 8: 2503.03821 = Eur. Phys. J. C 85 (2025). |
| Tedesco 2607.13228 adjacent to row 3 | fixed | Cited in row 3 as adjacent, no contradiction. |

## Verdict
All High and Medium issues from both skeptics are applied in `claim.md`. Not fixed: L4 (needs a same-beta run), L6 and STATE.md/`literature/` edits (overseer-owned), L8 and the HKS coefficient comparison (need the papers). The revised claim is weaker in its headline (candidate C, partial: row 3 open, SDP/bootstrap family not examined) and should go to a fresh skeptic before the overseer accepts it.

## Open issues
- `src/reproduce.py` lock is untested; first use should be the 027 repair runs after 027 finishes. Check that `artifacts/.reproduce.lock` is gitignored via `artifacts/` (it is inside it).
- 027's number check will report N1-N5, N9, N31 missing until the 028 M5 repair is run.
- Same-beta pair for N20/N21; GM prefactor; HKS coefficients; 2611.11534 body comparison.
- Overseer: update CLAIMS.md C5 (Lean), C9 printed intervals if they differ, STATE.md unvetted-proofs list, `literature/` from `new_references.md`.
