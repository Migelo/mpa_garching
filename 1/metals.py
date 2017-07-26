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

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    mask = s.gas['mass_at_ejection'] > 0
    metallicity_at_infall = s.gas['metals_at_infall'][mask] / s.gas['mass_at_infall'][mask]
    metallicity_at_ejection = s.gas['metals_at_ejection'][mask] / s.gas['mass_at_ejection'][mask]
    metals_infall =  s.gas['metals_at_infall'][mask]
    metals_ejection = s.gas['metals_at_ejection'][mask]

    f = plt.figure(figsize=utils.figsize)
    gs = gridspec.GridSpec(4, 1)
    ax1 = plt.subplot(gs[:2, 0])
    ax2 = plt.subplot(gs[2, 0])
    ax3 = plt.subplot(gs[3, 0])

    _, _, _, cbar = pg.plotting.scatter_map(np.log10(metallicity_at_infall),
        np.log10(metallicity_at_ejection), qty=s.gas['mass_at_infall'][mask],
        colors=s.gas['T_at_ejection'][mask], colors_av=s.gas['mass_at_infall'][mask],
        logscale=True, extent=[[-4, 0], [-4, 0]], clogscale=True, clim=[10**3.2, 1e8],
        zero_is_white=True, ax=ax1)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r'$log_{10} ( Temperature\ [K] ) $')
    ax1.plot([-4, 0], [-4, 0], color='r')
    ax1.set_xlim((-4, -1))
    ax1.set_ylim((-3.5, -.5))
    ax1.set_xlabel(r"$log_{10}(z_{at\ infall})$")
    ax1.set_ylabel(r"$log_{10}(z_{at\ ejection})$")

    ax2.set_ylabel("count")
    ax2.set_xlabel("metallicity")
    ax2.set_xlim((0, .02))
    ax2.hist(metallicity_at_infall, bins=np.linspace(0, .1, 101), alpha=.5, label='infall')
    ax2.hist(metallicity_at_ejection, bins=np.linspace(0, .1, 101), alpha=.5, label='ejection')
    leg1 = ax2.legend(loc='best', fontsize=14)
    plt.setp(ax2.get_xticklabels(), rotation=45)

    ax3.set_ylabel("count")
    ax3.set_xlabel("metallicity")
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.hist(metallicity_at_infall, bins=np.logspace(-4, -1, 31), alpha=.5, label='infall')
    ax3.hist(metallicity_at_ejection, bins=np.logspace(-4, -1, 31), alpha=.5, label='ejection')
    leg2 = ax3.legend(loc='best', fontsize=14)

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(7)
p.map(plot, utils.combinations)

