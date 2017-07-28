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
        MgII_profiles[0] = pg.analysis.radially_binned(s.gas, 'Mg/mass', av='MgII', r_edges=bins, proj=1)
        MgII_profiles[1] = pg.analysis.radially_binned(s.gas[m_rec], 'Mg/mass', av='MgII', r_edges=bins, proj=1)
        MgII_profiles[2] = pg.analysis.radially_binned(s.gas[m_nrec], 'Mg/mass', av='MgII', r_edges=bins, proj=1)
        OVI_profiles[0] = pg.analysis.radially_binned(s.gas, 'Mg/mass', av='OVI', r_edges=bins, proj=1)
        OVI_profiles[1] = pg.analysis.radially_binned(s.gas[m_rec], 'Mg/mass', av='OVI', r_edges=bins, proj=1)
        OVI_profiles[2] = pg.analysis.radially_binned(s.gas[m_nrec], 'Mg/mass', av='OVI', r_edges=bins, proj=1)
    else:
        MgII_profiles[0] = np.column_stack((MgII_profiles[0], pg.analysis.radially_binned(s.gas, 'Mg/mass', av='MgII', r_edges=bins, proj=1)))
        MgII_profiles[1] = np.column_stack((MgII_profiles[1], pg.analysis.radially_binned(s.gas[m_rec], 'Mg/mass', av='MgII', r_edges=bins, proj=1)))
        MgII_profiles[2] = np.column_stack((MgII_profiles[2], pg.analysis.radially_binned(s.gas[m_nrec], 'Mg/mass', av='MgII', r_edges=bins, proj=1)))
        OVI_profiles[0] = np.column_stack((OVI_profiles[0], pg.analysis.radially_binned(s.gas, 'Mg/mass', av='OVI', r_edges=bins, proj=1)))
        OVI_profiles[1] = np.column_stack((OVI_profiles[1], pg.analysis.radially_binned(s.gas[m_rec], 'Mg/mass', av='OVI', r_edges=bins, proj=1)))
        OVI_profiles[2] = np.column_stack((OVI_profiles[2], pg.analysis.radially_binned(s.gas[m_nrec], 'Mg/mass', av='OVI', r_edges=bins, proj=1)))

#y_label = pg.analysis.radially_binned(s.gas, 'metallicity', av='MgII', r_edges=bins, proj=1).units.latex()

# log and smooth
for i, item in enumerate(MgII_profiles):
    for j in range(item.shape[1]):
        MgII_profiles[i][:, j] = signal.savgol_filter(np.log10(MgII_profiles[i][:, j]/pg.solar.Mg()), 5, 3)
        OVI_profiles[i][:, j] = signal.savgol_filter(np.log10(OVI_profiles[i][:, j]/pg.solar.O()), 5, 3)

f = plt.figure(figsize=utils.figsize[::-1] * 2)
gs = gridspec.GridSpec(3, 2)
ax1 = plt.subplot(gs[:-1, 0])
ax2 = plt.subplot(gs[-1, 0])
ax3 = plt.subplot(gs[:-1, 1])
ax4 = plt.subplot(gs[-1, 1])
#ax5 = plt.subplot(gs[:-1, 2])
#ax6 = plt.subplot(gs[-1, 2])
#ax7 = plt.subplot(gs[:-1, 3])
#ax8 = plt.subplot(gs[-1, 3])

ax1.set_ylabel(r'$log_{10}[Z]_\odot$', fontsize=44)
ax2.set_ylabel(r'$Z_{non\ recycled}/Z_{recycled}$', fontsize=44)


for a in (ax1, ax3):
    a.tick_params(labelbottom='off')
for a in (ax2, ax4):
    a.set_ylim((1e-1, 1e1))
    a.set_yscale('log')
    a.plot([-200, 200], [1, 1], color='k')
    a.set_xlabel(r'$r\ [kpc]$', fontsize=44)

#ax1.set_ylim((-1, 7))
#ax3.set_ylim((-3, 5))
#ax5.set_ylim((-6, 2))
#ax7.set_ylim((-5, 3))

#ax1.set_title('All gas')
#ax3.set_title('HI')
ax1.set_title('MgII')
ax3.set_title('OVI')

ax1.plot(bins[:-1], np.percentile(MgII_profiles[0], 50, axis=1), color='b', label='total')
ax1.plot(bins[:-1], np.percentile(MgII_profiles[1], 50, axis=1), color='r', label='recycled')
ax1.plot(bins[:-1], np.percentile(MgII_profiles[2], 50, axis=1), color='g', label='non recycled')
ax1.fill_between(bins[:-1], np.percentile(MgII_profiles[0], 25, axis=1),
    np.percentile(MgII_profiles[0], 75, axis=1), alpha=.25, color='b')
ax1.fill_between(bins[:-1], np.percentile(MgII_profiles[1], 25, axis=1),
    np.percentile(MgII_profiles[1], 75, axis=1), alpha=.25, color='r')
ax1.fill_between(bins[:-1], np.percentile(MgII_profiles[2], 25, axis=1),
    np.percentile(MgII_profiles[2], 75, axis=1), alpha=.25, color='g')
ax1.legend(loc='upper right', fontsize=20)

ax3.plot(bins[:-1], np.percentile(OVI_profiles[0], 50, axis=1), color='b', label='total')
ax3.plot(bins[:-1], np.percentile(OVI_profiles[1], 50, axis=1), color='r', label='recycled')
ax3.plot(bins[:-1], np.percentile(OVI_profiles[2], 50, axis=1), color='g', label='non recycled')
ax3.fill_between(bins[:-1], np.percentile(OVI_profiles[0], 25, axis=1),
    np.percentile(OVI_profiles[0], 75, axis=1), alpha=.25, color='b')
ax3.fill_between(bins[:-1], np.percentile(OVI_profiles[1], 25, axis=1),
    np.percentile(OVI_profiles[1], 75, axis=1), alpha=.25, color='r')
ax3.fill_between(bins[:-1], np.percentile(OVI_profiles[2], 25, axis=1),
    np.percentile(OVI_profiles[2], 75, axis=1), alpha=.25, color='g')

ax2.plot(bins[:-1], 10**(np.percentile(MgII_profiles[2], 50, axis=1)) / 10**(np.percentile(MgII_profiles[1], 50, axis=1)))
ax4.plot(bins[:-1], 10**(np.percentile(OVI_profiles[2], 50, axis=1)) / 10**(np.percentile(OVI_profiles[1], 50, axis=1)))

ax2.set_xlim(ax1.get_xlim())
ax4.set_xlim(ax3.get_xlim())
#ax6.set_xlim(ax5.get_xlim())
#ax8.set_xlim(ax7.get_xlim())

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

#f, ax = plt.subplots(3, figsize=utils.figsize)
#for a in ax:
#    a.set_ylabel(r'$log_{10}\left(\rho\ [%s]\right)$' % y_label, fontsize=15)
#    a.set_ylim((1.5, 5))
#ax[0].set_title('All gas')
#ax[2].set_xlabel(r'$r\ [kpc]$', fontsize=20)
#
#ax[0].plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color='b', label='total')
#ax[0].plot(bins[:-1], np.percentile(profiles[1], 50, axis=1), color='b', ls='-.', lw=2, label=r'total $T>10^{5.2}$ K')
#ax[0].plot(bins[:-1], np.percentile(profiles[2], 50, axis=1), color='b', ls='dashed', label=r'total $T\leq10^{5.2}$ K')
#ax[0].fill_between(bins[:-1], np.percentile(profiles[0], 25, axis=1),
#    np.percentile(profiles[0], 75, axis=1), alpha=.25, color='b')
#ax[1].plot(bins[:-1], np.percentile(profiles[3], 50, axis=1), color='r', label='recycled')
#ax[1].plot(bins[:-1], np.percentile(profiles[4], 50, axis=1), color='r', ls='-.', lw=2, label=r'recycled $T>10^{5.2}$ K')
#ax[1].plot(bins[:-1], np.percentile(profiles[5], 50, axis=1), color='r', ls='dashed', label=r'recycled $T\leq10^{5.2}$ K')
#ax[1].fill_between(bins[:-1], np.percentile(profiles[3], 25, axis=1),
#    np.percentile(profiles[3], 75, axis=1), alpha=.25, color='r')
#ax[2].plot(bins[:-1], np.percentile(profiles[6], 50, axis=1), color='g', label='non recycled')
#ax[2].plot(bins[:-1], np.percentile(profiles[7], 50, axis=1), color='g', ls='-.', lw=2, label=r'non recycled $T>10^{5.2}$ K')
#ax[2].plot(bins[:-1], np.percentile(profiles[8], 50, axis=1), color='g', ls='dashed', label=r'non recycled $T\leq10^{5.2}$ K')
#ax[2].fill_between(bins[:-1], np.percentile(profiles[6], 25, axis=1),
#    np.percentile(profiles[6], 75, axis=1), alpha=.25, color='g')
#
#ax[0].legend(loc=1, fontsize=15)
#ax[1].legend(loc=1, fontsize=15)
#ax[2].legend(loc=1, fontsize=15)
#
#f.tight_layout()
#
#plt.savefig(filename.split("/")[-1][:-3] + "_all_gas_only.png", bbox_inches='tight')

