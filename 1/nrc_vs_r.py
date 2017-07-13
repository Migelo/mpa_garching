import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils
import glob
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    cycle_r_max = [item[item > 0] for item in s.gas['cycle_r_max'][s.gas['num_recycled'] > -1]]
    number_of_recycles = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    extended_nor = []
    
    sm = s.gas[np.max(s.gas['T_at_ejection'], axis=1) > 0]

    time_bins = [4, 3, 2, 1, .5, 0][::-1]
    mask = sm['ejection_a'].flatten() > 0
    ejection_a = sm['ejection_a'].flatten()[mask]
    ejection_z = pg.cosmology.a2z(ejection_a)
    metallicity_ejection = np.concatenate(sm['metals_at_ejection']/sm['mass_at_ejection'])[mask]


    average_metallicity, edges_met, _ = stats.binned_statistic(
        ejection_z, metallicity_ejection,
        bins=time_bins, statistic='mean')
    median_metallicity, edges_met, _ = stats.binned_statistic(
        ejection_z, metallicity_ejection,
        bins=time_bins, statistic='median')


    for i, particle in enumerate(cycle_r_max):
        for line in particle:
            extended_nor.append(number_of_recycles[i])

    average_nor, edges_nor, count_nor = stats.binned_statistic(np.concatenate(cycle_r_max), extended_nor, statistic='mean', bins=np.linspace(0, 200, 201))
    average_z, edges, count = stats.binned_statistic(np.concatenate(cycle_r_max), np.concatenate(metals_ejection), statistic='mean', bins=np.linspace(0, 200, 201))

    average_z = utils.prepare_step(average_z)
    average_nor = utils.prepare_step(average_nor)

    f, ax = plt.subplots(2, figsize=utils.figsize)
    ax[0].grid(True)
    ax[0].set_xlabel(('cycle r max [kpc]'))
    ax[0].set_ylabel('z', color='b')
    ax[0].step(edges, average_z)
    ax[0].tick_params('y', colors='b')

    ax2 = ax[0].twinx()
    ax2.step(edges_nor, average_nor, color='r')
    ax2.set_ylabel('average number of recycles', color='r')
    ax2.tick_params('y', colors='r')
    
    ax[1].grid(True)
    ax[1].set_xlim((0, 4))
    ax[1].set_xlabel('redshift')
    ax[1].set_ylabel('metallicity')
    ax[1].step(edges_met, utils.prepare_step(average_metallicity), label='average')
    ax[1].step(edges_met, utils.prepare_step(median_metallicity), label='median')
    ax[1].legend(loc='upper right')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

