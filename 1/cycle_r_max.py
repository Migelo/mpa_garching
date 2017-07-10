import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
import glob
from multiprocessing import Pool
from scipy import stats
import utils

filename = __file__

def cycle_r_max(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    cycle_r_max = s.gas['cycle_r_max'] #[s.gas['num_recycled'] > -1]
    particle_mass = s.gas['mass_at_ejection'] #[s.gas['num_recycled'] > -1]
    mask = np.max(s.gas['cycle_r_max'], axis=1) > '100 kpc'
    cycle_r_max_lim = cycle_r_max[mask]
    ID_mask = pg.IDMask(s.gas['ID'][mask])


    mass_in_bins, edges2, count2 = stats.binned_statistic(cycle_r_max.flatten(),
                                     particle_mass.flatten(), statistic='sum', bins=np.logspace(-1.1, 4, 90))
    mass_in_bins /= 1e10
    
    fig, ax = plt.subplots(4, figsize=(16,12))
    ax[0].set_xlabel("$r_{max}$ [kpc]")
    ax[0].set_ylabel("mass [$10^{10}\ M_{\odot}$]")
    ax[0].set_xscale('log')
    ax[0].set_xlim((1e-2, 1e4))
    ax[0].set_yscale('log')
    ax[0].grid(True)
    ax[0].bar(edges2[:-1], mass_in_bins, width=np.diff(edges2))

    ax[1].set_xlabel("$r_{max}$ [kpc]")
    ax[1].set_ylabel("count")
    ax[1].grid(True)
    ax[1].set_yscale('log')
    ax[1].hist(cycle_r_max.flatten(), bins=np.linspace(0, 300, 301))
    
    ax[2].set_title("same as above but logscale")
    ax[2].set_xlabel("$r_{max}$ [kpc]")
    ax[2].set_ylabel("count")
    ax[2].grid(True)
    ax[2].set_yscale('log')
    ax[2].set_xscale('log')
    ax[2].hist(cycle_r_max.flatten(), bins=np.logspace(-1.1, 4, 90), log=True)

    ax[3].set_title("same as above but just particles with max(cycle_r_max) > 100)")
    ax[3].set_xlabel("$r_{max}$ [kpc]")
    ax[3].set_ylabel("count")
    ax[3].grid(True)
    ax[3].set_yscale('log')
    ax[3].set_xscale('log')
    ax[3].hist(cycle_r_max_lim.flatten(), bins=np.logspace(-1.1, 4, 90), log=True)

    fig.tight_layout()
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')
    plt.clf()
    plt.close()
    #for snap in sorted(glob.glob(path))[70::20]:
    #    s, h, g = pg.prepare_zoom(snap, gas_trace=None, star_form=None)
    #    fig, ax = plt.subplots(1, 3)
    #    fig.suptitle(s.cosmic_time())
    #    pltargs = dict(clim=(10**3.5, 10**6.5))
    #    pg.plotting.image(s.gas, extent='400 kpc', ax=ax[0], **pltargs)
    #    pg.plotting.image(s.gas[ID_mask], extent='400 kpc', ax=ax[1], **pltargs)
    #    pg.plotting.image(s.gas[~ID_mask], extent='400 kpc', ax=ax[2], **pltargs)
    #    plt.tight_layout()
    #    print '%s/%s/%s_%s.png' % (halo, definition, definition, snap.split('/')[-1])
    #    plt.savefig('%s/%s/%s_%s.png' % (halo, definition, definition, snap.split('/')[-1]))
    #    plt.clf()
    #    plt.close()

p = Pool(4)
p.map(plot, utils.combinations)

