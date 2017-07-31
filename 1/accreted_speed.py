import pygad as pg
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import numpy as np
import pygad.plotting
import glob
import utils
from scipy import stats
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo, definition = args
    print halo    
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    
    timerange = np.arange(0, s.cosmic_time(), 0.5)

    mask = s.gas['num_recycled'] > -1
    s_m = s.gas[mask]
    mass_infall = s_m['mass_at_infall']
    mass_infall_initial = mass_infall[:, 0]
    mass_infall_reac = mass_infall[:, 1:]
    mass_ejection = s_m['mass_at_ejection'] 
    mass_ejection_initial = mass_ejection[:, 0]
    mass_ejection_reac = mass_ejection[:, 1:]
    
    infall_time = s_m['infall_time']
    infall_time_initial = infall_time[:, 0]
    infall_time_reac = infall_time[:, 1:]
    ejection_time = s_m['ejection_time']
    ejection_time_initial = ejection_time[:, 0]
    ejection_time_reac = ejection_time[:, 1:]
    
    mass_infall, _, _ = stats.binned_statistic(np.concatenate(infall_time), np.concatenate(mass_infall), statistic='sum', bins=timerange)
    mass_infall_initial, _, _ = stats.binned_statistic(infall_time_initial, mass_infall_initial, statistic='sum', bins=timerange)
    mass_infall_reac, _, _ = stats.binned_statistic(np.concatenate(infall_time_reac), np.concatenate(mass_infall_reac), statistic='sum', bins=timerange)
    mass_ejection, _, _ = stats.binned_statistic(np.concatenate(ejection_time), np.concatenate(mass_ejection), statistic='sum', bins=timerange)
    mass_ejection_initial, _, _ = stats.binned_statistic(ejection_time_initial, mass_ejection_initial, statistic='sum', bins=timerange)
    mass_ejection_reac, _, _ = stats.binned_statistic(np.concatenate(ejection_time_reac), np.concatenate(mass_ejection_reac), statistic='sum', bins=timerange)

    dt = timerange[1] * 1e9

    mass_infall /= dt
    mass_infall_initial /= dt
    mass_infall_reac /= dt
    mass_ejection /= dt

    total_infall = mass_infall.sum()
    halftime = 0
    for i, x in enumerate(timerange):
        if mass_infall[:i].sum() / total_infall >= .5:
            halftime = x
            break

    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.grid(True)
    ax.set_ylabel(r'Rate [$M_{\odot}\ yr^{-1}$]')
    ax.set_yscale('log')
    ax.set_ylim((1e-1, 2e2))
    ax.set_xlabel('Time [Gyr]')
    ax.step(timerange[:-1], mass_infall, label='Total infall')
    ax.step(timerange[:-1], mass_infall_initial, label='First infall')
    ax.step(timerange[:-1], mass_infall_reac, label='Reacreation')
    ax.plot([halftime, halftime], [0, np.max(mass_infall)], color='k')
    spacing, ha = 1, 'left'
    if halftime >= s.cosmic_time() / 2:
        spacing = -1
        ha='right'
    ax.text(halftime + spacing * .1, np.max(mass_infall) / 2, 'half of integrated total infall', fontsize=17, ha=ha)
    ax.legend(loc='upper right')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
arguments = [(halo, 'halo') for halo in utils.halos]
p.map(plot, arguments)

    



