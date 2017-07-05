import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import utils
import glob
from multiprocessing import Pool


def time_in(s, array):
    array = np.array(array[array > 0])
    return pg.UnitArr(np.r_[array[1:], float(s.cosmic_time())] - array, 'Gyr')[:-1]


filename = __file__
def plot(args):
    halo = args[0]
    type = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, type), star_form=None)

    data = [time_in(s, item) for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
    data = [item for sublist in data for item in sublist]

    fig, ax = plt.subplots(1)
    ax.set_xlabel("infall time [Gyr]")
    ax.set_ylabel("count")
    ax.hist(data, bins=np.linspace(0, 1.5, 50))

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

