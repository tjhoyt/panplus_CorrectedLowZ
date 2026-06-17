"""
Computes delta_mass: the host-mass shift applied to low-z SNe in
extract_lowz_bias_correction.py. Also generates the low-z/high-z mass
comparison figures backing that calculation.

Prospector SED-fit masses are a single, homogeneously derived mass estimate
available across the full Pantheon+SH0ES redshift range, so they serve as a
common reference point for comparing Pantheon+SH0ES against *itself* at low
vs. high z, rather than as ground truth in their own right. Pantheon+SH0ES
host masses (HOST_LOGMASS) disagree with Prospector masses at low redshift
(zHEL < 0.15), but agree with them at high redshift (zHEL >= 0.15) - so the
high-z sample's Pantheon+/Prospector offset is taken as the zero-point, and
delta_mass is the shift that brings the low-z Pantheon+/Prospector offset
onto that same zero-point:

    delta_mass = median(HOST_LOGMASS - prosp_mass)[highz] - median(HOST_LOGMASS - prosp_mass)[lowz]

The median (not the mean) is adopted for this offset: the distribution of
per-SN differences is long-tailed (a handful of SNe disagree by >1-2 dex),
which pulls the mean noticeably further from zero than the median in both
the low-z and high-z samples. Both statistics are computed and reported
below so that tail is visible, but the median is the more robust summary of
the bulk offset and is what's used downstream.

Requires PROSPECTOR_MASS_FILE, a per-SN table of Prospector host masses
(column 'prosp_mass', indexed by Union3 SN name) that is NOT included in
this repository - it underlies a still-unpublished mass analysis and will
be added here once that paper (Union3.1) is accepted. Until then this
script is provided for transparency on the method and is not runnable
as-is. ID_MATCH_FILE (a Pantheon+ CID <-> Union3 SN name crosswalk, needed
because Prospector masses are keyed by Union3 naming rather than Pantheon+
CID) is included in this repository - it carries no mass information.
Re-extracted from mass_compare_multipanel.ipynb.
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table

from stat_tools import calc_mean, calc_median

PROSPECTOR_MASS_FILE = 'all_TractorProsp_hosts_full_table.csv'  # not in this repo; see module docstring
ID_MATCH_FILE = 'panPlus_union3_IDmatch.csv'  # included in this repo


def use_style(name):
    try:
        plt.style.use(name)
    except OSError:
        pass


def mass_offset(zcut_mask, prosp_mass, cid_to_sn):
    # one row per CID: HOST_LOGMASS is identical across an SN's duplicate
    # photometry rows, so dedup to avoid double-weighting repeated SNe.
    pantheon_mass = pd.Series(np.array(t['HOST_LOGMASS'][zcut_mask]), index=np.array(t['CID'][zcut_mask]))
    pantheon_mass = pantheon_mass.groupby(pantheon_mass.index).first()

    sn_for_cid = cid_to_sn.reindex(pantheon_mass.index)
    matched = sn_for_cid.notna() & sn_for_cid.isin(prosp_mass.index)
    pantheon_mass = pantheon_mass[matched]
    sn_olap = sn_for_cid[matched]

    delta = pantheon_mass.values - prosp_mass.loc[sn_olap].values
    delta = pd.Series(delta, index=pantheon_mass.index)

    medi, medi_bs, n_medi, _ = calc_median(delta, sigclip=4, confidence_level=0.68)
    mean, mean_bs, n_mean, _ = calc_mean(delta, sigclip=4, confidence_level=0.68)
    return {
        'n_overlap': len(sn_olap),
        'n_clipped': int(n_medi),
        'median': medi,
        'median_ci68_lo': medi_bs.confidence_interval.low,
        'median_ci68_hi': medi_bs.confidence_interval.high,
        'mean': mean,
        'mean_ci68_lo': mean_bs.confidence_interval.low,
        'mean_ci68_hi': mean_bs.confidence_interval.high,
        '_delta': delta,
        '_pantheon_mass': pantheon_mass,
    }


def plot_mass_compare(stats, label, outfile, mass_lims=(8.1, 11.99), delta_lims=(-1.49, 1.49),
                       delta_bins=np.arange(-1.3, 1.3, 0.2)):
    """Re-extracted from mass_compare_multipanel.ipynb: mass-vs-mass scatter
    (top-left), delta-mass-vs-mass scatter (bottom-left) with median (solid)
    and mean (dashed) offset lines, and a delta-mass histogram (bottom-right)."""
    pantheon_mass, delta = stats['_pantheon_mass'], stats['_delta']
    prosp_mass_matched = pantheon_mass - delta

    use_style('publish')
    fig, axs = plt.subplots(2, 2, figsize=(8, 5),
                             gridspec_kw={'height_ratios': [2, 1.2], 'width_ratios': [3, 1.2],
                                          'wspace': 0, 'hspace': 0.02})

    ax = axs[0, 0]
    xline = yline = np.linspace(mass_lims[0] - 0.3, mass_lims[1] + 0.3)
    ax.plot(xline, yline, color='black', lw=1.7)
    ax.scatter(prosp_mass_matched, pantheon_mass, s=18, marker='o', ec='black', fc='C0', alpha=0.6)
    ax.set_ylabel(f'LogMass ({label})')
    ax.set_xlim(*mass_lims)
    ax.set_ylim(*mass_lims)

    ax = axs[1, 0]
    ax.scatter(prosp_mass_matched, delta, s=18, marker='o', ec='black', fc='C0', alpha=0.6)
    ax.set_xlabel('LogMass (Prospector)')
    ax.set_ylabel(r'$\Delta$ LogMass')
    ax.set_xlim(*mass_lims)
    ax.set_ylim(*delta_lims)
    ax.axhline(stats['median'], color='C0', ls=(0, (3, 3)), label=f"median = {stats['median']:.3f}")
    ax.axhline(stats['mean'], color='C3', ls=(0, (1, 1)), label=f"mean = {stats['mean']:.3f}")
    ax.legend(loc='lower right', frameon=True, fontsize=7)

    axs[0, 1].set_yticks([])
    axs[0, 1].set_xticks([])

    ax = axs[1, 1]
    ax.hist(delta, bins=delta_bins, orientation='horizontal', ec='black', fc='C0', alpha=0.6)
    ax.axhline(stats['median'], color='C0', ls=(0, (3, 3)))
    ax.axhline(stats['mean'], color='C3', ls=(0, (1, 1)))
    ax.set_ylim(*delta_lims)
    ax.set_yticklabels([])

    for ax in axs.flatten():
        ax.axhline(0, color='black', lw=0.5)
        ax.grid()

    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close(fig)


t = Table.read('Pantheon+SH0ES.dat', format='ascii')
not_calibrator = t['IS_CALIBRATOR'] == 0

prosp = pd.read_csv(PROSPECTOR_MASS_FILE, index_col='sn')
prosp_mass = prosp['prosp_mass']

id_match = pd.read_csv(ID_MATCH_FILE, index_col='CID')
cid_to_sn = id_match['sn']

lowz = mass_offset((t['zHEL'] < 0.15) & not_calibrator, prosp_mass, cid_to_sn)
highz = mass_offset((t['zHEL'] >= 0.15) & not_calibrator, prosp_mass, cid_to_sn)
delta_mass = highz['median'] - lowz['median']

plot_mass_compare(lowz, 'PanPlus, z<0.15', 'mass_compare_lowZ.png',
                   mass_lims=(8.1, 11.99), delta_lims=(-3.49, 3.49), delta_bins=np.arange(-2.3, 3.3, 0.35))
plot_mass_compare(highz, 'PanPlus, z>=0.15', 'mass_compare_highZ.png')

print(f"low-z (zHEL<0.15):   median = {lowz['median']:.4f} ({lowz['median_ci68_lo']:.4f}, {lowz['median_ci68_hi']:.4f})"
      f"   mean = {lowz['mean']:.4f} ({lowz['mean_ci68_lo']:.4f}, {lowz['mean_ci68_hi']:.4f})   N={lowz['n_clipped']}")
print(f"high-z (zHEL>=0.15): median = {highz['median']:.4f} ({highz['median_ci68_lo']:.4f}, {highz['median_ci68_hi']:.4f})"
      f"   mean = {highz['mean']:.4f} ({highz['mean_ci68_lo']:.4f}, {highz['mean_ci68_hi']:.4f})   N={highz['n_clipped']}")
print(f"delta_mass = highz_median - lowz_median = {delta_mass:.4f}")

result = {
    'description': 'delta_mass = highz_offset - lowz_offset, offset = median(HOST_LOGMASS - prosp_mass)',
    'lowz': {k: v for k, v in lowz.items() if not k.startswith('_')},
    'highz': {k: v for k, v in highz.items() if not k.startswith('_')},
    'delta_mass': delta_mass,
}
with open('mass_offset.json', 'w') as f:
    json.dump(result, f, indent=2)
