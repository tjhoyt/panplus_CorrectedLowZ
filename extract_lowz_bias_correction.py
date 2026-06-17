"""
Extracted from extract_all_parameters.ipynb: the pipeline that shifts low-z
host masses by delta_mass, recomputes the BS21 bias correction for SNe that
cross the M=10 mass threshold as a result, and writes out a corrected
Pantheon+SH0ES table.

Cells pulled from the notebook (in order): 0, 2, 43, 51, 60, 61.
Cells 52-59 sat between 51 and 60/61 in the notebook but are unrelated
diagnostic plots (one of them isn't even valid Python) - they don't feed
into t_modify, so they're omitted here.
"""
import json

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, MaskedColumn

# 'science'/'publish' are cosmetic styles (from the scienceplots package plus a
# local custom style); fall back silently if they aren't installed so the
# script still runs and produces correct .dat output.
def use_style(name):
    try:
        plt.style.use(name)
    except OSError:
        pass

# --- cell 2: load the Pantheon+SH0ES table ---
old_names = ['1994DRichmond', '2005df_ANU', '2008fv_comb']
new_names = ['1994D', '2005df', '2008fv']

t = Table.read('Pantheon+SH0ES.dat', format='ascii')
for i in range(len(old_names)):
    old, new = old_names[i], new_names[i]
    obj_ind = np.where(t['CID'] == [old])[0][0]
    t['CID'][obj_ind] = new

# Replace -9 host mass errors with a mask
t['HOST_LOGMASS_ERR'] = MaskedColumn(t['HOST_LOGMASS_ERR'], mask=t['HOST_LOGMASS_ERR']==-9)

# --- cell 43: low-z mass cut + cubic bias-correction-vs-color fits ---
use_style('science')
fig=plt.figure(figsize=(8,4))

# Make low redshift cut
zcut = (t['zCMB'] < 0.15) & (t['zCMB'] > 0.01)

x, y = t['c'], t['biasCor_m_b']

lo_cut, hi_cut = t['HOST_LOGMASS'] <= 10, t['HOST_LOGMASS'] > 10
colors = ['C0', 'C3']
labels = ['Low', 'High']
quad_coeffs, mass_biases = [], []
for i, mass_cut in enumerate([lo_cut, hi_cut]):
    cut = mass_cut & zcut
    print(np.sum(cut))
    xp, yp = x[cut], y[cut]
    quad_coeffs.append(np.polyfit(xp, yp, deg=3))
    mass_bias = np.median(yp)
    plt.scatter(xp,yp, fc=colors[i], ec='black', lw=0.2, label=labels[i],s=18)
    mass_biases.append(mass_bias)

x_curve = np.linspace(-0.35, 0.35, 500)
for i, coeffs in enumerate(quad_coeffs):
    func = np.poly1d(coeffs)
    plt.plot(x_curve, func(x_curve), color='black')

plt.xlim(-0.34, 0.34)
plt.ylim(-0.34, 0.34)
plt.axhline(color='black',lw=0.7)
plt.legend()
print(np.diff(mass_biases))
plt.xlabel('SALT $c$')
plt.ylabel('$\\delta_{bias,BS21}$ (mag)')

# --- cell 51: shift low-z host masses, recompute bias correction ---
# low-z table
t_lowz = t[zcut]

# Define new masses. delta_mass is the shift that brings the low-z
# Pantheon+/Prospector mass offset onto the same zero-point as the high-z
# Pantheon+/Prospector offset - see compute_mass_offset.py and mass_offset.json.
with open('mass_offset.json') as f:
    delta_mass = json.load(f)['delta_mass']
t_lowz['new_mass_lowz'] = t_lowz['HOST_LOGMASS'] + delta_mass

# Unpack the cubic bias corr functions
func_lo = np.poly1d(quad_coeffs[0])
func_hi = np.poly1d(quad_coeffs[1])
delta_bias_func = func_hi - func_lo

# Updated calculations
bias_corr_new = np.where((t_lowz['new_mass_lowz']>10) & (t_lowz['HOST_LOGMASS']<10),
                          t_lowz['biasCor_m_b']+delta_bias_func(t_lowz['c']),
                          t_lowz['biasCor_m_b'])
delta_biasCorr = bias_corr_new - t_lowz['biasCor_m_b']
delta_biasCorr_mean = np.mean(delta_biasCorr)
delta_biasCorr_avg = np.average(delta_biasCorr, weights=1/(t_lowz['m_b_corr_err_RAW']**2))

# Mask of SNe actually affected by the mass shift (crossed the M=10 threshold).
# NOTE: `change_cut` gets reassigned to the *unaffected* mask further down in
# this cell (for the "unaffected" scatter plot) - `affected_cut` is the stable
# name to use for anything downstream (cell 60, the correction table).
affected_cut = np.abs(delta_biasCorr) > 0

use_style('science')
use_style('publish')
fig=plt.figure(figsize=(8,4))

# Plot corrected P+ masses (all red/high mass)
change_cut = np.abs(delta_biasCorr) > 0
print(np.sum(change_cut))
x, y = t_lowz['c'], bias_corr_new
lo_cut, hi_cut = t_lowz['new_mass_lowz'] <= 10, t_lowz['new_mass_lowz'] > 10
colors = ['C3']
labels = ['High']
for i, mass_cut in enumerate([hi_cut]):
    cut = mass_cut & change_cut
    xp, yp = x[cut], y[cut]
    print(np.sum(cut))
    plt.scatter(xp,yp, fc=colors[i], ec='black', lw=0.5, label=f'{labels[i]} (Corrected P+ mass)',s=30)

# Plot original P+ masses (all blue/low mass)
change_cut = np.abs(delta_biasCorr) > 0
x, y = t_lowz['c'], t_lowz['biasCor_m_b']
lo_cut, hi_cut = t_lowz['HOST_LOGMASS'] <= 10, t_lowz['HOST_LOGMASS'] > 10
colors = ['C0']
labels = ['Low']
for i, mass_cut in enumerate([lo_cut]):
    cut = mass_cut & change_cut
    xp, yp = x[cut], y[cut]
    plt.scatter(xp,yp, fc='none', ec=colors[i], lw=0.5,s=30, label=f'{labels[i]} (Published P+ mass)')

change_cut = np.abs(delta_biasCorr) == 0
print(np.sum(change_cut))
x, y = t_lowz['c'], t_lowz['biasCor_m_b'] + delta_biasCorr
lo_cut, hi_cut = t_lowz['new_mass_lowz'] <= 10, t_lowz['new_mass_lowz'] > 10
colors = ['C0', 'C3']
labels = ['Low', 'High']
for i, mass_cut in enumerate([lo_cut, hi_cut]):
    cut = mass_cut & change_cut
    xp, yp = x[cut], y[cut]
    plt.scatter(xp,yp, fc=colors[i], ec='black', lw=0.2, label=f'{labels[i]} (unaffected)',s=5, alpha=0.5)

# Plot theory curves
x_curve = np.linspace(-0.35, 0.35, 500)
for i, coeffs in enumerate(quad_coeffs):
    func = np.poly1d(coeffs)
    plt.plot(x_curve, func(x_curve), color='black')

plt.xlim(-0.22, 0.32)
plt.ylim(-0.32, 0.32)
plt.axhline(color='black',lw=0.7)
plt.legend(frameon=True, framealpha=1, ncol=2)
print(np.diff(mass_biases))
plt.xlabel('SALT $c$')
plt.ylabel('$\\delta_{bias,BS21}$ (mag)')
plt.savefig('biascorrVScolor_panplus_correctedMasses.png',dpi=400, bbox_inches='tight')

# --- cell 60: apply the bias-correction shift to m_b_corr ---
t_modify = t.copy()
mask_final = np.zeros_like(zcut, dtype=bool)
mask_final[zcut] = affected_cut
t_modify['m_b_corr'][mask_final] += -delta_biasCorr[affected_cut]
# Check the new mbcorr values were shifted correctly (round to 10th decimal place to avoid float precision issues)
print(np.all(np.round(t[zcut]['m_b_corr'] - t_modify[zcut]['m_b_corr'] - delta_biasCorr, 10) == 0))

# Replace masks with exactly the fill values the DESI SN data text file has
t_modify['HOST_LOGMASS_ERR'] = t_modify['HOST_LOGMASS_ERR'].filled(-9)

# --- cell 61: write the corrected table ---
from astropy.io import ascii

outfile = 'Pantheon+SH0ES_HoytFixLowZbiasCorr.dat'
t_modify.write(outfile, format='ascii', overwrite=True, delimiter=' ')
Table.read(outfile, format='ascii')

# --- correction table: just the SNe whose bias correction actually changed ---
t_corr = Table()
t_corr['CID'] = t_lowz['CID'][affected_cut]
t_corr['old_mass'] = t_lowz['HOST_LOGMASS'][affected_cut]
t_corr['new_mass'] = t_lowz['new_mass_lowz'][affected_cut]
t_corr['old_biasCorr'] = t_lowz['biasCor_m_b'][affected_cut]
t_corr['new_biasCorr'] = bias_corr_new[affected_cut]
# t_lowz['m_b_corr'] is read off of `t` (untouched by t_modify), so it's still
# the original/published value here - safe to use as "old". "new" is derived
# the same way cell 60 derives t_modify's m_b_corr, via affected_cut/delta_biasCorr,
# rather than reading back out of t_modify (avoids relying on mask_final ordering).
t_corr['old_mbCorr'] = t_lowz['m_b_corr'][affected_cut]
t_corr['new_mbCorr'] = t_lowz['m_b_corr'][affected_cut] - delta_biasCorr[affected_cut]

corr_outfile = 'Pantheon+SH0ES_HoytLowZbiasCorr_corrections.dat'
t_corr.write(corr_outfile, format='ascii', overwrite=True, delimiter=' ')
Table.read(corr_outfile, format='ascii')
