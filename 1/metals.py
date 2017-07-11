import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
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

    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    metals_infall_full = copy(metals_infall)

    for i, temp in enumerate(metals_ejection):
        if len(temp) < len(metals_infall[i]):
            metals_infall[i] = metals_infall[i][:-1]
        elif len(temp) > len(metals_infall[i]):
            print 'hmm'

    metals_ejection = [item for sublist in metals_ejection for item in sublist]
    metals_infall = [item for sublist in metals_infall for item in sublist]
    metals_infall_full = [item for sublist in metals_infall_full for item in sublist]

    fig, ax = plt.subplots(3)
    ax[0].set_xlabel("z at infall")
    ax[0].set_ylabel("z at ejection")
    ax[0].set_xlim((0, .1))
    ax[0].set_ylim((0, .1))
    pg.plotting.scatter_map(metals_infall, metals_ejection, ax=ax[0], logscale=True, extent=[[0, .1], [0, .1]])
    ax[0].plot([0, 1], [0, 1], color='r')

    ax[1].set_ylabel("count")
    ax[1].set_xlabel("z")
    ax[1].set_xlim((0, .1))
    ax[1].hist(metals_infall_full, bins=np.linspace(0, .1, 101), alpha=.5, label='infall')
    ax[1].hist(metals_ejection, bins=np.linspace(0, .1, 101), alpha=.5, label='ejection')
    leg1 = plt.legend(loc='best')

    ax[2].set_ylabel("count")
    ax[2].set_xlabel("z")
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].hist(metals_infall_full, bins=np.logspace(-4, -1, 31), alpha=.5, label='infall')
    ax[2].hist(metals_ejection, bins=np.logspace(-4, -1, 31), alpha=.5, label='ejection')
    leg2 = plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

