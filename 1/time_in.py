import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import utils
import glob
from scipy import stats
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    s_m = s.gas[s.gas['num_recycled'] > -1]
    bins = np.linspace(0, s.cosmic_time(), 14*2)

    infall_time = s_m['infall_time']
    ejection_time = s_m['ejection_time']
    time_in = ejection_time - infall_time
    infall_mass = s_m['mass_at_infall']
   
    one_infall_mask = time_in[:, 0] < 0
    time_one_infall = time_in[one_infall_mask][:, 0] + s.cosmic_time()
    infall_mass_one = infall_mass[one_infall_mask]
    infall_mass_one = infall_mass_one[infall_mass_one > 0]

    time_multiple_infalls = time_in[~one_infall_mask]
    time_multiple_infalls = time_multiple_infalls[time_multiple_infalls != 0]
    time_multiple_infalls[time_multiple_infalls < 0] += s.cosmic_time()
    infall_mass_multiple = infall_mass[~one_infall_mask]
    infall_mass_multiple = infall_mass_multiple[infall_mass_multiple > 0]

    infall_mass_one_binned, _, _ = stats.binned_statistic(time_one_infall, infall_mass_one, bins=bins, statistic='sum')
    infall_mass_multiple_binned, _, _ = stats.binned_statistic(time_multiple_infalls, infall_mass_multiple, bins=bins, statistic='sum')

    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.set_xlabel(r"$Time$ in the galaxy [Gyr]")
    ax.set_ylabel(r"$\log_{10}\left(Mass\ \left[\mathrm{M_\odot}\right]\right)$")
    ax.bar(bins[:-1], np.log10(infall_mass_one_binned), width=np.diff(bins),
        alpha=.3, label="one infall", color='b', align='edge')
    ax.bar(bins[:-1], np.log10(infall_mass_multiple_binned),
        width=np.diff(bins), alpha=.3, label="multiple infalls", color='r',
        align='edge')
    ax.set_xlim((0, 10.2))
    ax.set_ylim((5, 11))
    ax.legend(loc="upper right")

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

