# IMPORTANT 
This is the only authoritative source for the P+ mass corrections published in my
Union3.1 paper. The derivative repository hosted by Choudhury26 and its associated Zenodo
are *not* a valid nor approved source of these files.

That author is distributing *without permission* an LLM-derived reverse reconstruction 
of my own published table. Do not cite that paper if using its associated data products 
as they have never measured a host galaxy mass or SN luminosity and cannot provide 
adequate support or context for the measurements that were in fact made by this author.

# panplus_CorrectedLowZ

A correction to the Pantheon+SH0ES bias corrections (and therefore the
corrected B-band magnitudes, `m_b_corr`) for low-redshift supernovae whose
host-galaxy mass sits close to the M = 10 (log solar mass) threshold used in
the BS21 bias-correction simulations.

## Background

The Pantheon+SH0ES bias correction (`biasCor_m_b`) is computed separately for
low-mass (`HOST_LOGMASS <= 10`) and high-mass (`HOST_LOGMASS > 10`) host
galaxies, as a cubic function of SALT2 color `c`. At low redshift
(`0.01 < zCMB < 0.15`), host masses are systematically underestimated, which
moves some supernovae onto the wrong side of the M = 10 threshold and gives
them the wrong bias correction.

This repository re-derives the cubic bias-correction-vs-color relations for
the low- and high-mass populations directly from the published Pantheon+SH0ES
catalog, shifts each low-z host mass by `delta_mass` (see "Mass offset"
below), and recomputes `biasCor_m_b` (and therefore `m_b_corr`) for any
supernova that crosses the M = 10 threshold once the shift is applied.
Supernovae that don't cross the threshold are left untouched.

## Contents

- `extract_lowz_bias_correction.py` - the correction pipeline. Reads
  `Pantheon+SH0ES.dat` and writes the two output files below.
- `Pantheon+SH0ES.dat` - the published Pantheon+SH0ES catalog (input).
- `Pantheon+SH0ES_HoytFixLowZbiasCorr.dat` - the full Pantheon+SH0ES catalog
  with `m_b_corr` updated for the affected low-z supernovae.
- `Pantheon+SH0ES_HoytLowZbiasCorr_corrections.dat` - just the affected
  supernovae, with their old/new host mass, bias correction, and `m_b_corr`.
- `biascorrVScolor_panplus_correctedMasses.png` - diagnostic plot of the
  bias-correction-vs-color fits and which supernovae were shifted.
- `verify_covariance_unchanged.py` - diagnostic showing the published
  Pantheon+SH0ES covariance matrix does not need to be recomputed alongside
  the bias-correction shift (see "Covariance matrix" below). Writes
  `verify_covariance_unchanged.png`.
- `compute_mass_offset.py` - derives `delta_mass` and generates the low-z/
  high-z mass comparison figures below (see "Mass offset" below). Not
  runnable as-is - see that section for why.
- `mass_offset.json` - the `delta_mass` value `extract_lowz_bias_correction.py`
  reads, with its full provenance (see "Mass offset" below).
- `mass_compare_lowZ.png` / `mass_compare_highZ.png` - Pantheon+SH0ES vs.
  Prospector host mass comparison, for the low-z and high-z samples used to
  derive `delta_mass`.
- `panPlus_union3_IDmatch.csv` - Pantheon+ CID <-> Union3 SN-name crosswalk,
  needed because Prospector masses are keyed by Union3 naming rather than
  Pantheon+ CID. Carries no mass information.
- `stat_tools.py` - small bootstrap-statistics helper used by
  `compute_mass_offset.py`.

## Usage

```bash
pip install numpy matplotlib astropy
python extract_lowz_bias_correction.py
python verify_covariance_unchanged.py
```

This regenerates both `.dat` output files and both `.png` plots in place.

`plt.style.use('science'/'publish')` calls are cosmetic only (from the
optional [SciencePlots](https://github.com/garrettj403/SciencePlots) package
plus a local custom style) and fall back silently if not installed - they
don't affect the corrected data.

## Method details

- Low-z cut: `0.01 < zCMB < 0.15`.
- For each of the low-mass and high-mass populations within that cut, fit a
  cubic polynomial of `biasCor_m_b` vs. SALT2 color `c`.
- Shift every low-z host mass by `delta_mass` dex (see "Mass offset").
- For supernovae whose shifted mass crosses from `<= 10` to `> 10`, replace
  `biasCor_m_b` with the high-mass cubic fit evaluated at that supernova's
  color, and propagate the change through to `m_b_corr`.
- All other supernovae (outside the low-z cut, or not crossing the
  threshold) are unchanged.

## Mass offset

Prospector SED-fit masses are a single, homogeneously derived mass estimate
available across the full Pantheon+SH0ES redshift range. What matters about
them here isn't their absolute value relative to some other reference - it's
that they put every Pantheon+SH0ES supernova's host mass on the same,
standardized system regardless of redshift. They're used as a yardstick to
compare Pantheon+SH0ES host masses against *themselves* at low vs. high z.

Pantheon+SH0ES host masses (`HOST_LOGMASS`) agree with Prospector masses at
high redshift (`zHEL >= 0.15`), but disagree at low redshift (`zHEL < 0.15`).
`delta_mass` is the shift that brings the low-z Pantheon+/Prospector offset
onto the same zero-point as the high-z offset:

```
delta_mass = median(HOST_LOGMASS - prosp_mass)[highz] - median(HOST_LOGMASS - prosp_mass)[lowz]
           = 0.0232 - (-0.5757)
           = 0.5989 dex
```

(`zHEL >= 0.15` / `< 0.15`, both excluding SH0ES calibrators, sigma-clipped
at 4 sigma; see `mass_offset.json` for the full numbers including bootstrap
confidence intervals and sample sizes).

The median (not the mean) is adopted for `delta_mass`. The per-SN
Pantheon+/Prospector differences are long-tailed - a handful of SNe disagree
by more than 1-2 dex - which pulls the mean noticeably further from zero
than the median in both samples:

| | median | mean | N |
|---|---|---|---|
| low-z (`zHEL<0.15`) offset | -0.576 | -0.795 | 414 |
| high-z (`zHEL>=0.15`) offset | +0.023 | +0.005 | 436 |

Both statistics are reported in `mass_offset.json` and plotted in
`mass_compare_lowZ.png` / `mass_compare_highZ.png` so the tail is visible,
but the median is the more robust summary of the bulk offset and is what's
used for `delta_mass`.

`compute_mass_offset.py` reproduces this calculation and regenerates the
comparison figures, but it depends on a per-SN table of Prospector host
masses that is **not included in this repository**: it underlies a still-
unpublished mass analysis (the Union3.1 paper, currently in review) and
will be added here once that paper is accepted. (`panPlus_union3_IDmatch.csv`,
the Pantheon+/Union3 ID crosswalk it also needs, carries no mass information
and is included.) Until the mass table is available, `compute_mass_offset.py`
is included for transparency on the method - it documents exactly how
`delta_mass` and the comparison figures were derived - but isn't runnable
as-is, and `mass_offset.json` plus the two `.png` figures (committed,
derived-numbers/plots only, no per-SN masses) are what's actually used
downstream.

## Covariance matrix

The published Pantheon+SH0ES statistical+systematic covariance matrix does
not need to be modified to use the corrected `m_b_corr` values here. The only
quantity that changes is the bias-correction term (`biasCor_m_b`), and its
uncertainty (`biasCorErr_m_b`) is small and effectively flat across the
low-z sample compared to the peculiar-velocity uncertainty
(`m_b_corr_err_VPEC`) already baked into the covariance matrix - by a median
factor of ~30x for `0.01 < zCMB < 0.15` (see `verify_covariance_unchanged.py`
and its output plot for the full comparison). In other words, the shift
changes the central values of `m_b_corr` for the affected supernovae, but is
well within the uncertainty budget the covariance matrix already accounts
for, so it is left unchanged.
