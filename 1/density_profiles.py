import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import pygad.plotting
import glob
import utils
from multiprocessing import Pool

filename = __file__
def plot(args):
    halo = args[0]
    type = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, type), star_form=None)
    
    f, ax = plt.subplots(3, 2, figsize=(20,12))
    
    for axes in ax.flatten():
        axes.grid(True)

    pg.plotting.profile(s.gas, '200 kpc', 'mass', ax=ax[0,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'mass', ax=ax[0,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'mass', ax=ax[0,0], label='recycled')
    ax[0,0].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'MgII', ax=ax[1,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'MgII', ax=ax[1,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'MgII', ax=ax[1,0], label='recycled')
    ax[1,0].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'CIV', ax=ax[0,1], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'CIV', ax=ax[0,1], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'CIV', ax=ax[0,1], label='recycled')
    ax[0,1].set_ylim((1e-4, 1e0))
    ax[0,1].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'HI', ax=ax[1,1], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'HI', ax=ax[1,1], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'HI', ax=ax[1,1], label='recycled')
    ax[1,1].legend(loc='best')
    
    pg.plotting.profile(s.gas, '200 kpc', 'OVI', ax=ax[2,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'OVI', ax=ax[2,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'OVI', ax=ax[2,0], label='recycled')
    ax[2,0].set_ylim((1e-3, 1e-1))
    ax[2,0].legend(loc='best')
    
    f.tight_layout()
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')
    
p = Pool(4)
p.map(plot, utils.combinations)
