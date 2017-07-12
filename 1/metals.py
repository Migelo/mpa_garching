import pygad as pg
import pygad.plotting
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
import utils
import glob
from multiprocessing import Pool
import matplotlib.gridspec as gridspec
from scipy import stats
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
    metals_infall_full = copy(metals_infall)

    for i, temp in enumerate(metals_ejection):
        if len(temp) < len(metals_infall[i]):
            metals_infall[i] = metals_infall[i][:-1]
        elif len(temp) > len(metals_infall[i]):
            print 'hmm'

    metals_ejection = [item for sublist in metals_ejection for item in sublist]
    metals_infall = [item for sublist in metals_infall for item in sublist]
    metals_infall_full = [item for sublist in metals_infall_full for item in sublist]

    count_infall, edges, foo = stats.binned_statistic(metals_infall,
        metals_ejection, statistic='count', bins=np.logspace(-4, -1, 31))
    count_ejection, edges, foo = stats.binned_statistic(metals_ejection,
        metals_infall, statistic='count', bins=np.logspace(-4, -1, 31))

#    normalisation = .1/float(np.max(count_infall))
    count_infall = np.array(count_infall).astype(float) * 1e-1 / np.max(count_infall)
    count_ejection = np.array(count_ejection).astype(float) * 1e-1 / np.max(count_ejection)


    f = plt.figure(figsize=(8.268, 8.268*2**.5))
    gs = gridspec.GridSpec(4, 1)
    ax1 = plt.subplot(gs[:2, 0])
    ax2 = plt.subplot(gs[2, 0])
    ax3 = plt.subplot(gs[3, 0])

    bins, x_edge, y_edge = np.histogram2d(metals_ejection, metals_infall,
        bins=[np.logspace(-4, 0, 61), np.logspace(-4, 0, 61)])
    vmin, vmax = np.percentile(utils.finite(np.log10(bins.flatten())), [0, 95])
    
    abc = ax1.pcolormesh(x_edge, y_edge, np.log10(bins), vmin=vmin, vmax=vmax)
    cax = inset_axes(ax1, width="70%", height="3%", loc=1)
    cbar = f.colorbar(abc, cax=cax, orientation='horizontal')
    cbar.set_label('$log_{10}(count)$', color='w')
    cbar_obj = plt.getp(cbar.ax.axes, 'xticklabels')
    plt.setp(cbar_obj, color='w')
    
    ax1.plot([0, 1], [0, 1], color='r')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim((1e-4, 1e0))
    ax1.set_ylim((1e-4, 1e0))
    ax1.set_xlabel("z at infall")
    ax1.set_ylabel("z at ejection")
    #ax1.bar(edges[:-1], count_infall, width=np.diff(edges), alpha=.3)
    #ax1.barh(edges[:-1], count_ejection, height=np.diff(edges), alpha=.3)
    #ax1.scatter(metals_infall, metals_ejection, alpha=.1, edgecolor=None)

    ax2.set_ylabel("count")
    ax2.set_xlabel("z")
    ax2.set_xlim((0, .1))
    ax2.hist(metals_infall, bins=np.linspace(0, .1, 101), alpha=.5, label='infall')
    ax2.hist(metals_ejection, bins=np.linspace(0, .1, 101), alpha=.5, label='ejection')
    leg1 = ax2.legend(loc='best')

    ax3.set_ylabel("count")
    ax3.set_xlabel("z")
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.hist(metals_infall, bins=np.logspace(-4, -1, 31), alpha=.5, label='infall')
    ax3.hist(metals_ejection, bins=np.logspace(-4, -1, 31), alpha=.5, label='ejection')
    leg2 = ax3.legend(loc='best')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

