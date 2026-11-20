# Sources for data/

All numbers were parsed from the LaTeX source of the arXiv e-print (`https://arxiv.org/e-print/<id>`,
fetched 2026-11-20) by `work/20261120-002-glueball-data/extract_tables.py`. Table numbers below
are the order of table floats in that source (equal to the printed numbering of the arXiv
version); the LaTeX `\label` is given too, since journal versions may renumber.
"Continuum" means the authors' a -> 0 extrapolation; "N=inf" is the authors' large-N fit.
No unit conversions were applied to stored values; the only derived numbers are in
`work/20261120-002-glueball-data/ratios.py` output, not in `data/`.

Error notation: `v(e)` means e in units of the last printed digit of v; `v(x.y)` is literal.

---

## MP1999 -> `glueball_su3_4d_morningstar_peardon_1999.csv`

- arXiv: hep-lat/9901004. C. J. Morningstar, M. Peardon, "The glueball spectrum from an anisotropic lattice study", 1999 (arXiv version; journal reference not verified here).
- Setting: 3+1D SU(3) pure gauge, improved anisotropic action, continuum limit, infinite-volume (authors' finite-volume checks).
- Table VII (`tab:finalmasses`): r0*m_G and m_G in MeV for 15 states. r0*m errors: (continuum-extrapolation stat)(anisotropy, ~1%). MeV errors: (combined r0*m)(from r0^-1 = 410(20) MeV). Column "Other J" stored in `note`; dagger (0*++, 2*-+) = interpretation tentative.
- Table VIII (`tab:massratios`): mass ratios (quantity `m/m(0++)` or `m/m(2++)`); errors assume channels uncorrelated (per text).
- Scale: r0^-1 = 410(20) MeV (Table VII caption). Stored as `1/r0`.

## Chen2006 -> `glueball_su3_4d_chen_2006.csv`

- arXiv: hep-lat/0510074. Y. Chen, A. Alexandru, S. J. Dong, T. Draper, I. Horvath, F. X. Lee, K. F. Liu, N. Mathur, C. Morningstar, M. Peardon, S. Tamhankar, B. L. Young, J. B. Zhang, "Glueball spectrum and matrix elements on anisotropic lattices", arXiv Oct 2005, journal Phys. Rev. D 73 (2006) 014516 (per independent verification; source id kept as Chen2006; journal reference not verified here).
- Setting: 3+1D SU(3), improved anisotropic action, 5 lattice spacings, continuum.
- Table XXI (`summary1`, final spectrum): r0*M_G (stat)(1% anisotropy) and M_G in MeV (combined r0*M)(r0^-1 = 410(20) MeV), 12 J^PC. Source typo `10.25(4))(10)` for 3-- read as 10.25(4)(10).
- Table VIII (`comparison`), column "this work" only: r0*M_G per cubic irrep R^PC (row_label) before combining into J. The column quoting MP1999 was not extracted (see discrepancies).
- Scale: r0^-1 = 410(20) MeV.

## AT2020 -> `glueball_su3_4d_athenodorou_teper_2020.csv`

- arXiv: 2007.06422. A. Athenodorou, M. Teper, "The glueball spectrum of SU(3) gauge theory in 3+1 dimensions", 2020 (arXiv version; journal reference not verified here).
- Setting: 3+1D SU(3), Wilson plaquette action, beta 5.6924-6.50, continuum, sqrt(sigma) = fundamental string tension.
- Table 17 (`table_MK_J`): M/sqrt(sigma) per J^PC, gs/ex1/ex2. `*` = J assignment "likely", `**` = "significant uncertainty" (Table 16 convention, `table_M_J_R`). Last row: r0*sqrt(sigma) = 1.160(6).
- Table 18 (`table_MGeV_J`): masses in GeV using r0 = 0.472(5) fm; errors statistical only.
- Scales: r0 = 0.472(5) fm = 1/418(5) MeV (eqn `eqn_r0MeV`, quoted from a lattice QCD review the paper cites as RS-2014); sqrt(sigma) = 1.160(6)/r0 = 485(6) MeV (eqn `eqn_KMeV`); "sqrt(sigma) ~ 440 MeV" from Regge slope + rotating string (text just before `eqn_r0MeV`, no error, explicitly rough). All three copied into `string_tension_scales_su3_4d.csv`.

## LTW2004 -> `glueball_sun_4d_lucini_teper_wenger_2004.csv`

- arXiv: hep-lat/0404008. B. Lucini, M. Teper, U. Wenger, "Glueballs and k-strings in SU(N) gauge theories: calculations with improved operators", 2004 (arXiv version; journal reference not verified here).
- Setting: 3+1D SU(N), N=2,3,4,6,8, Wilson action; SU(4) also on an anisotropic lattice (row `SU(4)_an`, stored as N=4 with note).
- Table 14 (`table_continuumG`): continuum m/sqrt(sigma) for 0++, 0++* (level 1), 2++, and N=inf from m(N) = m(inf) + c/N^2 (eqn `eqn_Nmk`). Star = poor fit.
- Not extracted: Table 15 (`table_l16G`, fixed a*sqrt(sigma) ~ 0.2, not continuum).

## AT2021 -> `glueball_sun_4d_athenodorou_teper_2021.csv`

- arXiv: 2106.00364. A. Athenodorou, M. Teper, "SU(N) gauge theories in 3+1 dimensions: glueball spectrum, string tensions and topology", 2021 (arXiv version; journal reference not verified here).
- Setting: 3+1D SU(N), N=2,3,4,5,6,8,10,12, Wilson action, continuum.
- Table 34 (`table_MKJ_N2-5`) and Table 35 (`table_MKJ_N6-12`): M(J^PC)/sqrt(sigma). Marks: per text (around the cubic-irrep tables) `*` = poor fit 2.5 < chi2/ndf < 3.5; the LaTeX source comments under these tables add J-assignment remarks (e.g. "su4 4++ * = could be 3++", "su2 2++ex1 + = only T2 for J=2"). Kept as "author mark" in `note`.
- Table 38 (`table_MK_J_SUN`): N=inf, fit M/sqrt(sigma) = M(inf) + c/N^2 (eqn `eqn_MKN`); 2++, 1-+, 1-- ground states exclude the smallest N (per text). Entries printed as `>= x(e)` are stored with quantity `m/sqrt_sigma_lower_bound`.
- r0*sqrt(sigma) = 1.160(6) and r0 = 0.472(5) fm are quoted from AT2020 (not duplicated here).
- Not extracted: Table 36 (`table_MKJ_lowbound`, N-independent bounds), Table 37 (`table_MJ_N3`, SU(3) GeV, duplicates AT2020 Table 18 plus bounds), cubic-irrep Tables 23-31.

## Teper1998 -> `glueball_sun_3d_teper_1998.csv`

- arXiv: hep-lat/9804008. M. J. Teper, "SU(N) gauge theories in 2+1 dimensions", 1998 (arXiv version; journal reference not verified here).
- Setting: 2+1D SU(N), N=2,3,4,5, Wilson action, continuum. In 2+1D g^2 has mass dimension 1.
- Table 28 (`table_mcont`): continuum m/sqrt(sigma); `*`, `**` = first, second excitation (level 1, 2).
- Table 34 (`table_masskN`): N=inf of m/sqrt(sigma), linear fit in 1/N^2; slope, N>= and CL in `note`.
- Table 33 (`table_massgN`): N=inf of m/(g^2 N), same fit form.

## AT2016 -> `glueball_sun_3d_athenodorou_teper_2016.csv`

- arXiv: 1609.03873. A. Athenodorou, M. Teper, "SU(N) gauge theories in 2+1 dimensions: glueball spectra and k-string tensions", 2016 (arXiv version; journal reference not verified here).
- Setting: 2+1D SU(N), N=2,3,4,6,8,12,16, continuum. sigma = sigma_f, the fundamental string tension. Parity doubling M(J^{+C}) = M(J^{-C}) for J != 0 holds in the continuum (stated in the paper, eqn following "parity doubling for states with J != 0"); both partners are listed separately.
- Table B1 (`table_sigf_cont`): sqrt(sigma_f)/(g^2 N) and the mass gap m_G/(g^2 N) (stored as JPC 0++, level 0), N=inf included; superscripts a,b = medium/large volume.
- Tables B3-B9 (`table_mksu{2,3,4,6,8,12,16}_cont`): M/sqrt(sigma) for all levels; level = row order within the J^PC block, chi2/ndf in `note`; for SU(8) the c1 coefficient column is not stored.
- Tables B10-B12 (`table_mksuN_J0/J2/J1`): N=inf, fits with c1/N^2 (+ c2/N^4 where given). Values in square brackets are averages of N=12,16, not extrapolations (flagged in `note`).
- Not extracted: Table B2 (k-string tension ratios).

- Table numbers are the printed arXiv appendix labels (B1-B12). Corrected 2026-11-20 after independent verification (work/20261120-004-data-verify); values unchanged.

---

## Gaps and transcription remarks

- No table could not be fetched; all seven papers were read from source.
- Chen2006 Table VIII quotes MP1999 values in its last column; for 1-- it prints 9.50(4), while MP1999 Table VII gives 1-- = 9.40(6) (9.50(4) is MP1999's 2*-+). Not extracted; noted only.
- AT2016 Table B4 (SU(3)) prints the 2-- level-3 row `12.09(12)` twice; both rows stored, the second flagged.
- AT2020 Table 17 vs AT2021 Table 34 for SU(3): identical except 2+- gs, 8.74(12)* vs 8.91(15); AT2021 also omits 1-+ ex1, ex2 and 4+- gs.
- Physical units: MP1999/Chen2006 use r0^-1 = 410(20) MeV; AT2020/2021 use r0 = 0.472(5) fm = 1/418(5) MeV. Do not mix MeV columns across sources without rescaling.
