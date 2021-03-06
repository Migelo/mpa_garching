import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import utils
import glob
from multiprocessing import Pool

filename = __file__
def plot(args):
    halo = args[0]
    type = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, type), star_form=None)

    mask = s.gas['ejection_time'][:,0] > '8 Gyr'
    ej_hot = s.gas['mass_at_ejection'][mask][ s.gas['T_at_ejection'][:,0][mask]>1e5 , 0 ].sum()
    ej_cold = s.gas['mass_at_ejection'][mask][ s.gas['T_at_ejection'][:,0][mask]<=1e5 , 0 ].sum()
    print ej_hot / (ej_cold+ej_hot)
    ej_hot = s.gas['mass_at_ejection'][ s.gas['T_at_ejection']>1e5 ].sum()
    ej_cold = s.gas['mass_at_ejection'][ s.gas['T_at_ejection']<=1e5 ].sum()
    print ej_hot / (ej_cold+ej_hot)

    T_ejection = s.gas['T_at_ejection'][s.gas['num_recycled'] > -1]
    T_infall = s.gas['T_at_infall'][s.gas['num_recycled'] > -1]

    for i, temp in enumerate(T_ejection):
        if len(temp) < len(T_infall[i]):
            T_infall[i] = T_infall[i][:-1]
        elif len(temp) > len(T_infall[i]):
            print 'hmm'

    fig, ax = plt.subplots(3)
    plt.tight_layout()
    ax[0].set_xlabel("$T_{infall}$")
    ax[0].set_ylabel("$T_{ejection}$")
    ax[0].set_xlim((1e2, 1e9))
    ax[0].set_ylim((1e2, 1e9))
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].scatter(T_infall, T_ejection, alpha=.1, edgecolor='none')
    ax[0].plot([0, 1e80], [0, 1e80], color='r')

    ax[1].set_ylabel("count")
    ax[1].set_xlabel("$T$ [K]")
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].hist(T_infall.flatten(), bins=np.logspace(2, 9, 100/2), alpha=.5, label='infall')
    ax[1].hist(T_ejection.flatten(), bins=np.logspace(2, 9, 100/2), alpha=.5, label='ejection')
    lgd1 = ax[1].legend(loc='best')

    ax[2].set_title('Initial')
    ax[2].set_ylabel("count")
    ax[2].set_xlabel("$T$ [K]")
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].hist(T_infall[:,0], bins=np.logspace(2, 9, 100/2), alpha=.5, label='initial infall')
    ax[2].hist(T_ejection[:,0], bins=np.logspace(2, 9, 100/2), alpha=.5, label='initial ejection')
    lgd2 = ax[2].legend(loc='best')


    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

