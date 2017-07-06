import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import utils
import glob
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    r = s.gas['r'][s.gas['num_recycled'] > -1] 
    nrc = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    z = s.gas['metallicity'][s.gas['num_recycled'] > -1]

    nrc_avg, edges_nrc, count_nrc = stats.binned_statistic(r, nrc, statistic='mean', bins=np.linspace(0, 3000, 3000))
    z_avg, edges_z, count_z = stats.binned_statistic(r, z, statistic='mean', bins=np.linspace(0, 3000, 3000))
    
    nrc_avg = utils.prepare_step(nrc_avg)
    z_avg = utils.prepare_step(z_avg)

    fig, (ax1, ax2) = plt.subplots(2)
    plt.tight_layout()
    ax1.grid(True)
    ax1.set_xlim((0, 100))
    ax1.set_ylim((0, np.max(nrc[r < ax1.get_xlim()[1]])))
    ax1.set_xlabel('r [kpc]')
    ax1.set_ylabel('number of recycles')
    ax1.scatter(r, nrc)
    ax1.step(edges_nrc, nrc_avg, c='r')

    ax2.grid(True)
    ax2.set_xlim((0, 100))
    ax2.set_ylim((5e-4, 1e-1))
    ax2.set_yscale(('log'))
    ax2.set_xlabel('r [kpc]')
    ax2.set_ylabel('z')
    ax2.scatter(r, z)
    ax2.step(edges_z, z_avg, c='r')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

