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
catalog, shifts each low-z host mass by `delta_mass = 0.6 dex`, and
recomputes `biasCor_m_b` (and therefore `m_b_corr`) for any supernova that
crosses the M = 10 threshold once the shift is applied. Supernovae that don't
cross the threshold are left untouched.

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
- Shift every low-z host mass by `delta_mass = 0.6 dex`.
- For supernovae whose shifted mass crosses from `<= 10` to `> 10`, replace
  `biasCor_m_b` with the high-mass cubic fit evaluated at that supernova's
  color, and propagate the change through to `m_b_corr`.
- All other supernovae (outside the low-z cut, or not crossing the
  threshold) are unchanged.

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
