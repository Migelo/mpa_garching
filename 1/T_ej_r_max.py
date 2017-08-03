import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils
import glob
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    
    mask = np.max(s.gas['T_at_ejection'], axis=-1) > 0
    s_masked = s.gas[mask]
    
    f, ax = plt.subplots(2, figsize=utils.figsize)

    _, _, _, cbar = pg.plotting.scatter_map(np.log10(s_masked['cycle_r_max'].flatten()),
        np.log10(s_masked['T_at_ejection'].flatten()), s=s_masked,
        extent=[[.5, 3.5], [3, 8]], logscale=True,
        qty=s_masked['mass_at_ejection'].flatten(),
        colors='ejection_time.flatten()',
        colors_av=s_masked['mass_at_ejection'].flatten(),
        clogscale=False, clim=[2, s.cosmic_time()],
        zero_is_white=True, ax=ax[0])
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel('Ejection time [Gyr]', fontsize=16)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)    
    ax[0].set_ylabel(r'$T_\mathrm{ejection}\ \mathrm{[K]}$')
    ax[0].tick_params(labelbottom='off')
    
    _, _, _, cbar = pg.plotting.scatter_map(np.log10(s_masked['cycle_r_max'].flatten()),
        np.log10(s_masked['T_at_ejection'].flatten()), s=s_masked,
        extent=[[.5, 3.5], [3, 8]], logscale=True,
        qty=s_masked['mass_at_ejection'].flatten(),
        colors=s_masked['metals_at_ejection'].flatten()/s_masked['mass_at_ejection'].flatten()/pg.solar.Z(),
        colors_av=s_masked['mass_at_ejection'].flatten(),
        clogscale=True, clim=[10**-1.5, 10**.5],
        zero_is_white=True, ax=ax[1])
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r'$\log_{10} [ Z ]_\odot $', fontsize=16)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)    
    ax[1].set_xlabel(r'$\log_{10}\left(\mathrm{cycle }r_\mathrm{max}\ \mathrm{[kpc]}\right)$')
    ax[1].set_ylabel(r'$\log_{10}\left(T_\mathrm{ejection}\ \mathrm{[K]}\right)$')
    
    f.tight_layout()
    plt.subplots_adjust(top=0.95)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(7)
p.map(plot, utils.combinations)

