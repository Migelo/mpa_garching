import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
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

    bins = np.linspace(0, s.cosmic_time(), 14*2)

    infall_time = s.gas['infall_time']
    ejection_time = s.gas['ejection_time']
    time_in = ejection_time - infall_time
    time_in = time_in[s.gas['num_recycled'] > -1]
    
    time_in_non_rcy = time_in[:, 0][time_in[:,0] < 0] + s.cosmic_time()
    time_in_non_rcy = time_in_non_rcy[time_in_non_rcy > 0]

    time_in_rcy = time_in[:, 1:]
    time_in_rcy[time_in_rcy < 0] += s.cosmic_time()
    time_in_rcy = time_in_rcy[time_in_rcy > 0]
    
    f, ax = plt.subplots(1)
    ax.set_xlabel("time in the galaxy [Gyr]")
    ax.set_ylabel("count")
#    ax.set_ylim((1e0, ax.get_ylim()[-1]))
    ax.hist(time_in_rcy, bins=bins, alpha=.3, label="recycled", log=True)
    ax.hist(time_in_non_rcy, bins=bins, alpha=.3, label="non recycled", log=True)
    ax.legend(loc="upper right")

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

