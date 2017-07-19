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

f = plt.figure(figsize=utils.figsize)
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[1, 0])
ax4 = plt.subplot(gs[1, 1])

ax1.set_ylabel(r'$\rho\ [M_\odot\ kpc^{-3}]$', fontsize=20)
ax3.set_ylabel(r'$\rho\ [M_\odot\ kpc^{-3}\ Gyr^{-1}]$', fontsize=20)

for a in [ax1, ax2, ax3, ax4]:
    a.set_yscale('log')
    a.set_xlabel('r [kpc]')

def binned_data(x, y, bins, mode=0, HI=False):
    data = []
    x = np.array(x)
    y = np.array(y)
    for i, bin in enumerate(bins[1:]):
        if not HI:
            if mode == 0:
                data_in_bins_tot[i] += y[(x >= bins[i]) & (x < bin)].tolist()
            if mode == 1:
                data_in_bins_rec[i] += y[(x >= bins[i]) & (x < bin)].tolist()
            if mode == 2:
                data_in_bins_nrec[i] += y[(x >= bins[i]) & (x < bin)].tolist()
        if HI:
            if mode == 0:
                HI_data_in_bins_tot[i] += y[(x >= bins[i]) & (x < bin)].tolist()
            if mode == 1:
                HI_data_in_bins_rec[i] += y[(x >= bins[i]) & (x < bin)].tolist()
            if mode == 2:
                HI_data_in_bins_nrec[i] += y[(x >= bins[i]) & (x < bin)].tolist()


bins = np.arange(0, 200, 2)
data_in_bins_tot = [[] for bin in bins[1:]]
data_in_bins_rec = [[] for bin in bins[1:]]
data_in_bins_nrec = [[] for bin in bins[1:]]
HI_data_in_bins_tot = [[] for bin in bins[1:]]
HI_data_in_bins_rec = [[] for bin in bins[1:]]
HI_data_in_bins_nrec = [[] for bin in bins[1:]]
for args in utils.combinations:
    halo = args[0]
    definition = 'ism'

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    dens = s.gas['rho']
    m_nrec = s.gas['num_recycled'] == -1
    m_rec = np.max(s.gas['T_at_ejection'] > 0, axis=1)
    binned_data(s.gas['r'], dens, bins, mode=0)
    binned_data(s.gas[m_nrec]['r'], dens[m_nrec], bins, mode=2)
    binned_data(s.gas[m_rec]['r'], dens[m_rec], bins, mode=1)

    HI_dens = s.gas.get('HI/dV')
    binned_data(s.gas['r'], HI_dens, bins, mode=0, HI=True)
    binned_data(s.gas[m_nrec]['r'], HI_dens[m_nrec], bins, mode=2, HI=True)
    binned_data(s.gas[m_rec]['r'], HI_dens[m_rec], bins, mode=1, HI=True)
    

HI_median, HI_percentile25, HI_percentile75 = [], [], []
median_tot = [np.nanpercentile(x, 50) for x in data_in_bins_tot]
percentile25_tot = [np.nanpercentile(x, 25) for x in data_in_bins_tot]
percentile75_tot = [np.nanpercentile(x, 75) for x in data_in_bins_tot]
median_rec = [np.nanpercentile(x, 50) for x in data_in_bins_rec]
percentile25_rec = [np.nanpercentile(x, 25) for x in data_in_bins_rec]
percentile75_rec = [np.nanpercentile(x, 75) for x in data_in_bins_rec]
median_nrec = [np.nanpercentile(x, 50) for x in data_in_bins_nrec]
percentile25_nrec = [np.nanpercentile(x, 25) for x in data_in_bins_nrec]
percentile75_nrec = [np.nanpercentile(x, 75) for x in data_in_bins_nrec]
HI_median_tot = [np.nanpercentile(x, 50) for x in HI_data_in_bins_tot]
HI_percentile25_tot = [np.nanpercentile(x, 25) for x in HI_data_in_bins_tot]
HI_percentile75_tot = [np.nanpercentile(x, 75) for x in HI_data_in_bins_tot]
HI_median_rec = [np.nanpercentile(x, 50) for x in HI_data_in_bins_rec]
HI_percentile25_rec = [np.nanpercentile(x, 25) for x in HI_data_in_bins_rec]
HI_percentile75_rec = [np.nanpercentile(x, 75) for x in HI_data_in_bins_rec]
HI_median_nrec = [np.nanpercentile(x, 50) for x in HI_data_in_bins_nrec]
HI_percentile25_nrec = [np.nanpercentile(x, 25) for x in HI_data_in_bins_nrec]
HI_percentile75_nrec = [np.nanpercentile(x, 75) for x in HI_data_in_bins_nrec]

 
ax1.plot(bins[:-1], median_tot, color='b', label='median')
ax1.plot(bins[:-1], median_rec, color='r', label='recycled')
ax1.plot(bins[:-1], median_nrec, color='g', label='non recycled')
ax1.fill_between(bins[:-1], percentile25_tot, percentile75_tot, alpha=.25, color='b', label='total')
ax1.fill_between(bins[:-1], percentile25_rec, percentile75_rec, alpha=.25, color='r', label='recycled')
ax1.fill_between(bins[:-1], percentile25_nrec, percentile75_nrec, alpha=.25, color='g', label='non recycled')
ax1.legend(loc='upper right')
ax2.plot(bins[:-1], HI_median_tot, color='b', label='median')
ax2.plot(bins[:-1], HI_median_rec, color='r', label='recycled')
ax2.plot(bins[:-1], HI_median_nrec, color='g', label='non recycled')
ax2.fill_between(bins[:-1], HI_percentile25_tot, HI_percentile75_tot, alpha=.25, color='b', label='total')
ax2.fill_between(bins[:-1], HI_percentile25_rec, HI_percentile75_rec, alpha=.25, color='r', label='recycled')
ax2.fill_between(bins[:-1], HI_percentile25_nrec, HI_percentile75_nrec, alpha=.25, color='g', label='non recycled')
ax2.legend(loc='upper right')

ax3.plot(bins[:-1], np.array(median_nrec) / np.array(median_rec))
ax3.plot([0, 200], [1, 1], color='k')
ax4.plot(bins[:-1], np.array(HI_median_nrec) / np.array(HI_median_rec))
ax4.plot([0, 200], [1, 1], color='k')

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

