"""
Diagnostic backing the standard reply to "did you also update the covariance
matrix?": at low z the bias-correction uncertainty is tiny and essentially
identical between the low/high-mass functional forms, while the peculiar
velocity uncertainty (already baked into the existing covariance matrix)
dwarfs it. Produces verify_covariance_unchanged.png for use in correspondence.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def use_style(name):
    try:
        plt.style.use(name)
    except OSError:
        pass

t = Table.read('Pantheon+SH0ES.dat', format='ascii')
zcut = (t['zCMB'] < 0.15) & (t['zCMB'] > 0.01)
t_lowz = t[zcut]

bias_err = t_lowz['biasCorErr_m_b']
vpec_err = t_lowz['m_b_corr_err_VPEC']
lo_mass = t_lowz['HOST_LOGMASS'] < 10
hi_mass = t_lowz['HOST_LOGMASS'] >= 10

# --- summary stats printed to stdout, for quoting in the email ---
ratio = vpec_err / np.where(bias_err == 0, np.nan, bias_err)
print(f'N low-z SNe (0.01 < z < 0.15): {len(t_lowz)}')
print(f'median bias-correction error: {np.median(bias_err):.4f} mag')
print(f'median peculiar-velocity error: {np.median(vpec_err):.4f} mag')
print(f'median VPEC/bias-correction-error ratio: {np.nanmedian(ratio):.1f}x')
print(f'fraction of bias-correction errors < 0.02 mag: {np.mean(bias_err < 0.02):.3f}')
print(f'median bias-correction error, low mass (logM<10): {np.median(bias_err[lo_mass]):.4f} mag')
print(f'median bias-correction error, high mass (logM>=10): {np.median(bias_err[hi_mass]):.4f} mag')

# --- figure ---
use_style('science')
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))

# Panel A: bias-correction error vs. VPEC error, as a function of redshift.
# VPEC error falls off with z (it's ~constant sigma_v converted to a magnitude
# error at distance cz); the bias-correction error stays small and flat.
ax = axs[0]
ax.scatter(t_lowz['zHEL'], vpec_err, s=12, fc='C0', ec='none', alpha=0.6,
           label='Peculiar velocity error\n(already in covariance)')
ax.scatter(t_lowz['zHEL'], bias_err, s=12, fc='C3', ec='none', alpha=0.6,
           label='Bias correction error\n(mass-step term)')
ax.axhline(0.02, color='black', ls='--', lw=1, label='0.02 mag')
ax.set_yscale('log')
ax.set_xlabel(r'$z_{\rm HEL}$')
ax.set_ylabel('Magnitude uncertainty [mag]')
ax.legend(frameon=True, fontsize=8)
ax.set_title('Bias-correction error vs. peculiar velocity error')

# Panel B: bias-correction error distribution, split by mass bin, to show the
# two functional forms (low-mass / high-mass fits) give essentially the same
# error budget.
ax = axs[1]
bins = np.linspace(0, 0.05, 30)
ax.hist(bias_err[lo_mass], bins=bins, histtype='step', lw=1.5, color='C0',
        density=True, label=f'Low mass (median {np.median(bias_err[lo_mass]):.3f})')
ax.hist(bias_err[hi_mass], bins=bins, histtype='step', lw=1.5, color='C3',
        density=True, label=f'High mass (median {np.median(bias_err[hi_mass]):.3f})')
ax.set_xlabel('Bias correction error [mag]')
ax.set_ylabel('Density')
ax.legend(frameon=True, fontsize=8)
ax.set_title('Bias-correction error: low vs. high mass')

fig.tight_layout()
fig.savefig('verify_covariance_unchanged.png', dpi=300, bbox_inches='tight')
