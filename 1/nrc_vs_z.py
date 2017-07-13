import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils
import glob
from multiprocessing import Pool
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    infall_time = [item[item > 0] for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
    ejection_time = [item[item > 0] for item in s.gas['ejection_time'][s.gas['num_recycled'] > -1]]
    number_of_recycles = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    most_recent_z = []

    for i, item in enumerate(metals_infall):
        metals = np.vstack([np.column_stack([item, infall_time[i]]), np.column_stack([metals_ejection[i], ejection_time[i]])])
        most_recent_z.append(metals[np.argmax(metals[:, 1])][0])

    average_nor, edges_nor, count_nor = stats.binned_statistic(most_recent_z, number_of_recycles, statistic='mean', bins=np.linspace(0, .1, 100))
    
    upper_lim = np.percentile(number_of_recycles, 99.98) + 2
    bins, y_edge, x_edge = np.histogram2d(number_of_recycles, most_recent_z,
        bins=[np.logspace(-4, -1, 31), np.arange(-.5, upper_lim, 1)][::-1])
    vmin, vmax = np.percentile(utils.finite(np.log10(bins.flatten())), [0, 95])

    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    abc = ax.pcolormesh(x_edge, y_edge, np.log10(bins), vmin=vmin, vmax=vmax)
    cax = inset_axes(ax, width="70%", height="3%", loc=1)
    cbar = f.colorbar(abc, cax=cax, orientation='horizontal')
    cbar.set_label('$log_{10}(count)$', color='w')
    cbar_obj = plt.getp(cbar.ax.axes, 'xticklabels')
    plt.setp(cbar_obj, color='w')
    ax.set_ylim((y_edge[0], y_edge[-1]))
    ax.set_xlabel("z")
    ax.set_ylabel("number of recycles")
    ax.set_xscale('log')
    ax.plot(edges_nor[: -1], average_nor, color='k')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)    

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

