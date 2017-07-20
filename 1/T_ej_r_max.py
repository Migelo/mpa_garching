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
    
    s_masked = s.gas[np.max(s.gas['T_at_ejection'], axis=-1) > 0]
    
    f, ax = plt.subplots(2, figsize=utils.figsize)
    pg.plotting.scatter_map(s_masked['cycle_r_max'].flatten(),
        np.log10(s_masked['T_at_ejection'].flatten()),
        extent=[[0, 2000], [3, 8]], logscale=True, ax=ax[0],
        vlim=[1e0, 1e2])
    ax[0].set_xlabel('cycle $r_{max}$ [kpc]')
    ax[0].set_ylabel(r'$log_{10}(T\ at\ ejection)$ [K]')
    
    pg.plotting.scatter_map(np.log10(s_masked['cycle_r_max'].flatten()),
        np.log10(s_masked['T_at_ejection'].flatten()),
        extent=[[0.5, 3.5], [3, 8]], logscale=True, ax=ax[1])
    ax[1].set_xlabel(r'$log_{10}\left(cycle\ r_{max}\right)$ [kpc]')
    ax[1].set_ylabel(r'$log_{10}\left(T\ at\ ejection\right)$ [K]')
    
    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

