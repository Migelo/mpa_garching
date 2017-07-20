import pygad as pg
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import numpy as np
import pygad.plotting
import glob
import utils
from scipy import stats
from multiprocessing import Pool

filename = __file__


bins = np.arange(0, 200, 5)

profiles = np.zeros(3, dtype=object)
HI_profiles = np.zeros(3, dtype=object)

for i, args in enumerate(utils.combinations):
    halo = args[0]
    definition = 'ism'

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    m_nrec = s.gas['num_recycled'] == -1
    m_rec = np.max(s.gas['T_at_ejection'] > 0, axis=1)
    
    if i == 0:
        profiles[0] = pg.analysis.profile_dens(s.gas, 'rho', r_edges=bins)
        profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], 'rho', r_edges=bins)
        profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], 'rho', r_edges=bins)
        HI_profiles[0] = pg.analysis.profile_dens(s.gas, 'HI/dV', r_edges=bins)
        HI_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], 'HI/dV', r_edges=bins)
        HI_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], 'HI/dV', r_edges=bins)
    else:
        profiles[0] = np.column_stack((profiles[0], pg.analysis.profile_dens(s.gas, 'rho', r_edges=bins)))
        profiles[1] = np.column_stack((profiles[1], pg.analysis.profile_dens(s.gas[m_rec], 'rho', r_edges=bins)))
        profiles[2] = np.column_stack((profiles[2], pg.analysis.profile_dens(s.gas[m_nrec], 'rho', r_edges=bins)))
        HI_profiles[0] = np.column_stack((HI_profiles[0], pg.analysis.profile_dens(s.gas, 'HI/dV', r_edges=bins)))
        HI_profiles[1] = np.column_stack((HI_profiles[1], pg.analysis.profile_dens(s.gas[m_rec], 'HI/dV', r_edges=bins)))
        HI_profiles[2] = np.column_stack((HI_profiles[2], pg.analysis.profile_dens(s.gas[m_nrec], 'HI/dV', r_edges=bins)))
 
f = plt.figure(figsize=utils.figsize)
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[1, 0])
ax4 = plt.subplot(gs[1, 1])

ax1.set_ylabel(r'$\rho\ [M_\odot\ kpc^{-3}]$', fontsize=20)
ax3.set_ylabel(r'$\rho_{accreted}/\rho_{recycled}$', fontsize=20)

ax1.set_yscale('log')
ax2.set_yscale('log')
ax1.set_ylim((1e-2, 1e8))
ax2.set_ylim((1e-5, 1e8))

for a in [ax1, ax2, ax3, ax4]:
    a.set_yscale('log')
    a.set_xlabel('r [kpc]')

ax1.plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color='b', label='total ')
ax1.plot(bins[:-1], np.percentile(profiles[1], 50, axis=1), color='r', label='recycled')
ax1.plot(bins[:-1], np.percentile(profiles[2], 50, axis=1), color='g', label='non recycled')
ax1.fill_between(bins[:-1], np.percentile(profiles[0], 25, axis=1),
    np.percentile(profiles[0], 75, axis=1), alpha=.25, color='b')
ax1.fill_between(bins[:-1], np.percentile(profiles[1], 25, axis=1),
    np.percentile(profiles[1], 75, axis=1), alpha=.25, color='r')
ax1.fill_between(bins[:-1], np.percentile(profiles[2], 25, axis=1),
    np.percentile(profiles[2], 75, axis=1), alpha=.25, color='g')
ax1.legend(loc='upper right')

ax2.plot(bins[:-1], np.percentile(HI_profiles[0], 50, axis=1), color='b', label='total ')
ax2.plot(bins[:-1], np.percentile(HI_profiles[1], 50, axis=1), color='r', label='recycled')
ax2.plot(bins[:-1], np.percentile(HI_profiles[2], 50, axis=1), color='g', label='non recycled')
ax2.fill_between(bins[:-1], np.percentile(HI_profiles[0], 25, axis=1),
    np.percentile(HI_profiles[0], 75, axis=1), alpha=.25, color='b')
ax2.fill_between(bins[:-1], np.percentile(HI_profiles[1], 25, axis=1),
    np.percentile(HI_profiles[1], 75, axis=1), alpha=.25, color='r')
ax2.fill_between(bins[:-1], np.percentile(HI_profiles[2], 25, axis=1),
    np.percentile(HI_profiles[2], 75, axis=1), alpha=.25, color='g')
ax2.legend(loc='upper right')


ax3.plot(bins[:-1], np.percentile(profiles[2], 50, axis=1) / np.percentile(profiles[1], 50, axis=1))
ax3.plot([0, 200], [1, 1], color='k')
ax3.set_ylim((1e-2, 1e2))
ax4.plot(bins[:-1], np.percentile(HI_profiles[2], 50, axis=1) / np.percentile(HI_profiles[1], 50, axis=1))
ax4.plot([0, 200], [1, 1], color='k')
ax4.set_ylim((1e-2, 1e2))

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

