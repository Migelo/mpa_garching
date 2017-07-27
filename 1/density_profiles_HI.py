import pygad as pg
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import numpy as np
import pygad.plotting
import glob
import utils
from scipy import stats
from scipy import signal
from multiprocessing import Pool

filename = __file__

bins = np.arange(0, 205, 5)

profiles = np.zeros(9, dtype=object)
HI_profiles = np.zeros(3, dtype=object)
MgII_profiles = np.zeros(3, dtype=object)
OVI_profiles = np.zeros(3, dtype=object)

for i, args in enumerate(utils.combinations):
    halo = args[0]
    definition = 'ism'

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    m_nrec = s.gas['num_recycled'] == -1
    m_rec = np.max(s.gas['T_at_ejection'] > 0, axis=1)
    hot = pg.ExprMask("temp > '10**5.2 K'")
    y_label = ''    
    
    if i == 0:
        profiles[0] = pg.analysis.profile_dens(s.gas, 'mass', r_edges=bins)
        profiles[1] = pg.analysis.profile_dens(s.gas[hot], 'mass', r_edges=bins)
        profiles[2] = pg.analysis.profile_dens(s.gas[~hot], 'mass', r_edges=bins)
        profiles[3] = pg.analysis.profile_dens(s.gas[m_rec], 'mass', r_edges=bins)
        profiles[4] = pg.analysis.profile_dens(s.gas[m_rec][hot], 'mass', r_edges=bins)
        profiles[5] = pg.analysis.profile_dens(s.gas[m_rec][~hot], 'mass', r_edges=bins)
        profiles[6] = pg.analysis.profile_dens(s.gas[m_nrec], 'mass', r_edges=bins)
        profiles[7] = pg.analysis.profile_dens(s.gas[m_nrec][hot], 'mass', r_edges=bins)
        profiles[8] = pg.analysis.profile_dens(s.gas[m_nrec][~hot], 'mass', r_edges=bins)
        HI_profiles[0] = pg.analysis.profile_dens(s.gas, 'HI', r_edges=bins)
        HI_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], 'HI', r_edges=bins)
        HI_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], 'HI', r_edges=bins)
        MgII_profiles[0] = pg.analysis.profile_dens(s.gas, 'MgII', r_edges=bins)
        MgII_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], 'MgII', r_edges=bins)
        MgII_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], 'MgII', r_edges=bins)
        OVI_profiles[0] = pg.analysis.profile_dens(s.gas, 'OVI', r_edges=bins)
        OVI_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], 'OVI', r_edges=bins)
        OVI_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], 'OVI', r_edges=bins)
    else:
        profiles[0] = np.column_stack((profiles[0], pg.analysis.profile_dens(s.gas, 'mass', r_edges=bins)))
        profiles[1] = np.column_stack((profiles[1], pg.analysis.profile_dens(s.gas[hot], 'mass', r_edges=bins)))
        profiles[2] = np.column_stack((profiles[2], pg.analysis.profile_dens(s.gas[~hot], 'mass', r_edges=bins)))
        profiles[3] = np.column_stack((profiles[3], pg.analysis.profile_dens(s.gas[m_rec], 'mass', r_edges=bins)))
        profiles[4] = np.column_stack((profiles[4], pg.analysis.profile_dens(s.gas[m_rec][hot], 'mass', r_edges=bins)))
        profiles[5] = np.column_stack((profiles[5], pg.analysis.profile_dens(s.gas[m_rec][~hot], 'mass', r_edges=bins)))
        profiles[6] = np.column_stack((profiles[6], pg.analysis.profile_dens(s.gas[m_nrec], 'mass', r_edges=bins)))
        profiles[7] = np.column_stack((profiles[7], pg.analysis.profile_dens(s.gas[m_nrec][hot], 'mass', r_edges=bins)))
        profiles[8] = np.column_stack((profiles[8], pg.analysis.profile_dens(s.gas[m_nrec][~hot], 'mass', r_edges=bins)))
        HI_profiles[0] = np.column_stack((HI_profiles[0], pg.analysis.profile_dens(s.gas, 'HI', r_edges=bins)))
        HI_profiles[1] = np.column_stack((HI_profiles[1], pg.analysis.profile_dens(s.gas[m_rec], 'HI', r_edges=bins)))
        HI_profiles[2] = np.column_stack((HI_profiles[2], pg.analysis.profile_dens(s.gas[m_nrec], 'HI', r_edges=bins)))
        MgII_profiles[0] = np.column_stack((MgII_profiles[0], pg.analysis.profile_dens(s.gas, 'MgII', r_edges=bins)))
        MgII_profiles[1] = np.column_stack((MgII_profiles[1], pg.analysis.profile_dens(s.gas[m_rec], 'MgII', r_edges=bins)))
        MgII_profiles[2] = np.column_stack((MgII_profiles[2], pg.analysis.profile_dens(s.gas[m_nrec], 'MgII', r_edges=bins)))
        OVI_profiles[0] = np.column_stack((OVI_profiles[0], pg.analysis.profile_dens(s.gas, 'OVI', r_edges=bins)))
        OVI_profiles[1] = np.column_stack((OVI_profiles[1], pg.analysis.profile_dens(s.gas[m_rec], 'OVI', r_edges=bins)))
        OVI_profiles[2] = np.column_stack((OVI_profiles[2], pg.analysis.profile_dens(s.gas[m_nrec], 'OVI', r_edges=bins)))

y_label = pg.analysis.profile_dens(s.gas, 'mass', r_edges=bins).units.latex()

# log and smooth
for i, item in enumerate(profiles):
    for j in range(item.shape[1]):
        profiles[i][:, j] = signal.savgol_filter(np.log10(profiles[i][:, j]), 5, 3)
        if i < 3:
            HI_profiles[i][:, j] = signal.savgol_filter(np.log10(HI_profiles[i][:, j]), 5, 3)
            MgII_profiles[i][:, j] = signal.savgol_filter(np.log10(MgII_profiles[i][:, j]), 5, 3)
            OVI_profiles[i][:, j] = signal.savgol_filter(np.log10(OVI_profiles[i][:, j]), 5, 3)

f = plt.figure(figsize=utils.figsize[::-1] * 2)
gs = gridspec.GridSpec(3, 4)
ax1 = plt.subplot(gs[:-1, 0])
ax2 = plt.subplot(gs[-1, 0])
ax3 = plt.subplot(gs[:-1, 1])
ax4 = plt.subplot(gs[-1, 1])
ax5 = plt.subplot(gs[:-1, 2])
ax6 = plt.subplot(gs[-1, 2])
ax7 = plt.subplot(gs[:-1, 3])
ax8 = plt.subplot(gs[-1, 3])

ax1.set_ylabel(r'$log_{10}\left(\rho\ [%s]\right)$' % y_label, fontsize=44)
ax2.set_ylabel(r'$\rho_{non\ recycled}/\rho_{recycled}$', fontsize=44)


for a in (ax1, ax3, ax5, ax7):
    a.tick_params(labelbottom='off')
for a in (ax2, ax4, ax6, ax8):
    a.set_ylim((1e-1, 1e1))
    a.set_yscale('log')
    a.plot([-200, 200], [1, 1], color='k')
    a.set_xlabel(r'$r\ [kpc]$', fontsize=44)

ax1.set_ylim((-1, 7))
ax3.set_ylim((-3, 5))
ax5.set_ylim((-6, 2))
ax7.set_ylim((-5, 3))

ax1.set_title('All gas')
ax3.set_title('HI')
ax5.set_title('MgII')
ax7.set_title('OVI')

ax1.plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color='b', label='total')
ax1.plot(bins[:-1], np.percentile(profiles[1], 50, axis=1), color='b', ls='-.', lw=2, label='total T>10^5.2 K')
ax1.plot(bins[:-1], np.percentile(profiles[2], 50, axis=1), color='b', ls='dashed', label='total T<=10^5.2 K')
ax1.plot(bins[:-1], np.percentile(profiles[3], 50, axis=1), color='r', label='recycled')
ax1.plot(bins[:-1], np.percentile(profiles[6], 50, axis=1), color='g', label='non recycled')
ax1.fill_between(bins[:-1], np.percentile(profiles[0], 25, axis=1),
    np.percentile(profiles[0], 75, axis=1), alpha=.25, color='b')
ax1.fill_between(bins[:-1], np.percentile(profiles[3], 25, axis=1),
    np.percentile(profiles[3], 75, axis=1), alpha=.25, color='r')
ax1.fill_between(bins[:-1], np.percentile(profiles[6], 25, axis=1),
    np.percentile(profiles[6], 75, axis=1), alpha=.25, color='g')
ax1.legend(loc='upper right', fontsize=20)

ax3.plot(bins[:-1], np.percentile(HI_profiles[0], 50, axis=1), color='b', label='total')
ax3.plot(bins[:-1], np.percentile(HI_profiles[1], 50, axis=1), color='r', label='recycled')
ax3.plot(bins[:-1], np.percentile(HI_profiles[2], 50, axis=1), color='g', label='non recycled')
ax3.fill_between(bins[:-1], np.percentile(HI_profiles[0], 25, axis=1),
    np.percentile(HI_profiles[0], 75, axis=1), alpha=.25, color='b')
ax3.fill_between(bins[:-1], np.percentile(HI_profiles[1], 25, axis=1),
    np.percentile(HI_profiles[1], 75, axis=1), alpha=.25, color='r')
ax3.fill_between(bins[:-1], np.percentile(HI_profiles[2], 25, axis=1),
    np.percentile(HI_profiles[2], 75, axis=1), alpha=.25, color='g')

ax5.plot(bins[:-1], np.percentile(MgII_profiles[0], 50, axis=1), color='b', label='total')
ax5.plot(bins[:-1], np.percentile(MgII_profiles[1], 50, axis=1), color='r', label='recycled')
ax5.plot(bins[:-1], np.percentile(MgII_profiles[2], 50, axis=1), color='g', label='non recycled')
ax5.fill_between(bins[:-1], np.percentile(MgII_profiles[0], 25, axis=1),
    np.percentile(MgII_profiles[0], 75, axis=1), alpha=.25, color='b')
ax5.fill_between(bins[:-1], np.percentile(MgII_profiles[1], 25, axis=1),
    np.percentile(MgII_profiles[1], 75, axis=1), alpha=.25, color='r')
ax5.fill_between(bins[:-1], np.percentile(MgII_profiles[2], 25, axis=1),
    np.percentile(MgII_profiles[2], 75, axis=1), alpha=.25, color='g')

ax7.plot(bins[:-1], np.percentile(OVI_profiles[0], 50, axis=1), color='b', label='total')
ax7.plot(bins[:-1], np.percentile(OVI_profiles[1], 50, axis=1), color='r', label='recycled')
ax7.plot(bins[:-1], np.percentile(OVI_profiles[2], 50, axis=1), color='g', label='non recycled')
ax7.fill_between(bins[:-1], np.percentile(OVI_profiles[0], 25, axis=1),
    np.percentile(OVI_profiles[0], 75, axis=1), alpha=.25, color='b')
ax7.fill_between(bins[:-1], np.percentile(OVI_profiles[1], 25, axis=1),
    np.percentile(OVI_profiles[1], 75, axis=1), alpha=.25, color='r')
ax7.fill_between(bins[:-1], np.percentile(OVI_profiles[2], 25, axis=1),
    np.percentile(OVI_profiles[2], 75, axis=1), alpha=.25, color='g')

ax2.plot(bins[:-1], 10**(np.percentile(profiles[6], 50, axis=1)) / 10**(np.percentile(profiles[3], 50, axis=1)))
ax4.plot(bins[:-1], 10**(np.percentile(HI_profiles[2], 50, axis=1)) / 10**(np.percentile(HI_profiles[1], 50, axis=1)))
ax6.plot(bins[:-1], 10**(np.percentile(MgII_profiles[2], 50, axis=1)) / 10**(np.percentile(MgII_profiles[1], 50, axis=1)))
ax8.plot(bins[:-1], 10**(np.percentile(OVI_profiles[2], 50, axis=1)) / 10**(np.percentile(OVI_profiles[1], 50, axis=1)))

ax2.set_xlim(ax1.get_xlim())
ax4.set_xlim(ax3.get_xlim())
ax6.set_xlim(ax5.get_xlim())
ax8.set_xlim(ax7.get_xlim())

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
ax.set_ylabel(r'$log_{10}\left(\rho\ [%s]\right)$' % y_label, fontsize=44)
ax.set_xlabel(r'$r\ [kpc]$', fontsize=44)
ax.set_ylim((1.5, 5))
ax.set_title('All gas')
ax.plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color='b', label='total')
ax.plot(bins[:-1], np.percentile(profiles[1], 50, axis=1), color='b', ls='-.', lw=2, label='total T>10^5.2 K')
ax.plot(bins[:-1], np.percentile(profiles[2], 50, axis=1), color='b', ls='dashed', label='total T<=10^5.2 K')
ax.plot(bins[:-1], np.percentile(profiles[3], 50, axis=1), color='r', label='recycled')
ax.plot(bins[:-1], np.percentile(profiles[4], 50, axis=1), color='r', ls='-.', lw=2, label='recycled T>10^5.2 K')
ax.plot(bins[:-1], np.percentile(profiles[5], 50, axis=1), color='r', ls='dashed', label='recycled T<=10^5.2 K')
ax.plot(bins[:-1], np.percentile(profiles[6], 50, axis=1), color='g', label='non recycled')
ax.plot(bins[:-1], np.percentile(profiles[7], 50, axis=1), color='g', ls='-.', lw=2, label='non recycled T>10^5.2 K')
ax.plot(bins[:-1], np.percentile(profiles[8], 50, axis=1), color='g', ls='dashed', label='non recycled T<=10^5.2 K')
ax.fill_between(bins[:-1], np.percentile(profiles[0], 25, axis=1),
    np.percentile(profiles[0], 75, axis=1), alpha=.25, color='b')
ax.fill_between(bins[:-1], np.percentile(profiles[3], 25, axis=1),
    np.percentile(profiles[3], 75, axis=1), alpha=.25, color='r')
ax.fill_between(bins[:-1], np.percentile(profiles[6], 25, axis=1),
    np.percentile(profiles[6], 75, axis=1), alpha=.25, color='g')
ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + "_all_gas_only.png", bbox_inches='tight')

