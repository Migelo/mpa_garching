import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import utils
import glob
from scipy import stats
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    s_m = s.gas[s.gas['num_recycled'] > -1]
    mask = s_m['mass_at_ejection'] > 0
    bins = np.linspace(0, s.cosmic_time(), 14*2)

    z_ej = s_m['metals_at_ejection'][mask] / s_m['mass_at_ejection'][mask]
    r_max = s_m['cycle_r_max'][mask]

    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    _, _, _, cbar = pg.plotting.scatter_map(np.log10(r_max), np.log10(z_ej),
        qty=s_m['mass_at_ejection'][mask], logscale=True,
        extent=[[.5, 4], [-4, -.5]], bins=25, vlim=[10**6, 10**9.49],
        ax=ax, zero_is_white=True)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r'$\log_{10}\left( Mass\ \mathrm{[M_\odot]} \right)$', fontsize=utils.axes_labelsize*.8)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=utils.tick_labelsize*.8)
    text = cbar_ax.xaxis.get_label().set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    cbar_ax.xaxis.set_label(text)
    ax.grid(False)
    ax.set_ylabel(r'$\log_{10}\left( Z_\mathrm{ejection} \right)$',
        fontsize=utils.axes_labelsize)
    ax.set_xlabel(r'$\log_{10}\left( r_\mathrm{max} \right)$',
        fontsize=utils.axes_labelsize)
    plt.setp(ax.get_xticklabels(), fontsize=utils.tick_labelsize)    
    plt.setp(ax.get_yticklabels(), fontsize=utils.tick_labelsize)    

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=32)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png")

p = Pool(7)
p.map(plot, utils.combinations)
    
