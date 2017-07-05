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

    x = range(max(s.gas['num_recycled']) + 1)
    y = [s.gas['mass'][s.gas['num_recycled'] == i].sum() for i in range(max(s.gas['num_recycled']) + 1)]
    y = np.array(y)
    fig, ax = plt.subplots(2)
    plt.tight_layout()

    ax[0].set_xlabel("$N_{recycled}$")
    ax[0].set_ylabel("$M_{\odot}$")
    ax[1].set_xlabel("$N_{recycled}$")
    ax[1].set_ylabel("$M_{\odot}$")
    ax[1].set_xlim((0, len(y[y > .1 * max(y)])))
    ax[1].set_xticks(range(int(ax[1].get_xlim()[-1])))

    ax[0].bar(x, y, width=1)
    ax[1].bar(x, y, width=1)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

