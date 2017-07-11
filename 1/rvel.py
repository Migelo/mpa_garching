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
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    
    f, ax = plt.subplots(3, 2, figsize=(20, 12))
    
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], log=False, label='recycled')
    #ax[0,0].set_ylim((1e4, 1e8))
    ax[0,0].legend(loc='upper right')

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], log=False, label='recycled')
    #ax[0,1].set_ylim((1e0, 1e8))
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], log=False, label='recycled')
    #ax[1,0].set_ylim((1e-4, 1e5))

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], log=False, label='recycled')
    #ax[1,1].set_ylim((1e-1, 1e2))

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], log=False, label='recycled')
    #ax[2,0].set_ylim((1e0, 1e3))
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], log=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], log=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], log=False, label='recycled')
    #ax[2,1].set_ylim((1e-1, 1e3))

    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
#    f.tight_layout()
    
    for axes in ax.flatten():
        axes.grid(True)
        axes.set_ylim((-150, 150))
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".pdf", bbox_inches='tight')
    
p = Pool(8)
p.map(plot, utils.combinations)

