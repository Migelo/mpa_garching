import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import pygad.plotting

filename = __file__
for halo in ("M0977", "M1192"):
    for type in ('disc-Uebler', 'disc', 'ball', 'ism'):
        s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_470' % (halo, halo), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, type))
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
        
