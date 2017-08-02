import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
import utils
import glob
from multiprocessing import Pool
import matplotlib.gridspec as gridspec
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


filename = __file__

bins = np.linspace(-4, -1, 20)

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    mask = s.gas['mass_at_ejection'] > 0
    s_m = s.gas[mask]
    metallicity_at_infall = s_m.gas['metals_at_infall'] / s_m.gas['mass_at_infall']
    metallicity_at_ejection = s_m.gas['metals_at_ejection'] / s_m.gas['mass_at_ejection']
    metals_infall =  s_m.gas['metals_at_infall']
    metals_ejection = s_m.gas['metals_at_ejection']
    mass_infall =  s_m.gas['metals_at_infall']
    mass_ejection = s_m.gas['metals_at_ejection']

    metallicity_in_hist, _, _ = stats.binned_statistic(np.log10(metallicity_at_infall),
        mass_infall, statistic='sum', bins=bins)
    metallicity_ej_hist, _, _ = stats.binned_statistic(np.log10(metallicity_at_ejection),
        mass_ejection, statistic='sum', bins=bins)

    f = plt.figure(figsize=utils.figsize)
    gs = gridspec.GridSpec(3, 1)
    ax1 = plt.subplot(gs[:2, 0])
    ax2 = plt.subplot(gs[2, 0])

    _, _, _, cbar = pg.plotting.scatter_map(np.log10(metallicity_at_infall),
        np.log10(metallicity_at_ejection), qty=s_m['mass_at_infall'],
        colors=s_m['T_at_ejection'], colors_av=s_m['mass_at_infall'],
        logscale=True, extent=[[-4, 0], [-4, 0]], clogscale=True,
        clim=[10**3.2, 1e8], zero_is_white=True, ax=ax1)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r'$\log_{10} ( T_\mathrm{ejection} [%s]) $' % s.gas['T_at_ejection'].units.latex(),
        fontsize=23)
    plt.setp(cbar_ax.xaxis.get_majorticklabels(), fontsize=20)
    cbar_ax.set_xticklabels(cbar_ax.xaxis.get_majorticklabels()[:-1])
    ax1.plot([-4, 0], [-4, 0], color='r')
    ax1.set_xlim((-4, -1))
    ax1.set_ylim((-3.5, -.5))
    ax1.set_xlabel(r"$\log_{10}(Z_\mathrm{at\ infall})$", fontsize=23)
    ax1.set_ylabel(r"$\log_{10}(Z_\mathrm{at\ ejection})$", fontsize=23)
    ax1.tick_params(labelbottom='off')
    plt.setp(ax1.yaxis.get_majorticklabels(), fontsize=20)

    ax2.set_ylabel(r"$Mass\ \mathrm{[10^{8}M_\odot]}$", fontsize=23)
    ax2.set_xlabel(r"$\mathrm{\log_{10}}\mathrm{(}Z\mathrm{)}$", fontsize=23)
    ax2.set_xlim((-4, -1))
    ax2.bar(bins[:-1], metallicity_in_hist / 1e8, width=np.diff(bins), alpha=.5, label='infall', align='edge')
    ax2.bar(bins[:-1], metallicity_ej_hist / 1e8, np.diff(bins), alpha=.5, label='ejection', align='edge')
    plt.setp(ax2.xaxis.get_majorticklabels(), fontsize=20)
    plt.setp(ax2.yaxis.get_majorticklabels(), fontsize=20)
    leg2 = ax2.legend(loc='best', fontsize=14)

    f.tight_layout()
    plt.subplots_adjust(top=0.95)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(7)
p.map(plot, utils.combinations)

