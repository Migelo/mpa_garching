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
    
    time_at_r_max = pg.cosmology.a2z(s.gas['cycle_r_max_at'])
    nrc = s.gas['num_recycled']
    
    start, stop, step = 0, 200, 5
    bins=np.arange(start, stop + step, step)
    time_bins = [4, 3, 2, 1, .5, 0]
        
    f, ax = plt.subplots(7, figsize=(8, 8/2*2**.5*7))

    y_min = []
    y_max = []
    y2_min = []
    y2_max = []
    second_axes = []
    for i, item in enumerate(time_bins[1:]):
        bin_start = item
        bin_end = time_bins[i]
        
        mask = (time_at_r_max <= bin_end) & (time_at_r_max > bin_start)
        mask2 = np.sum(mask, axis=1)
        mask3 = mask2 * 1.
        mask3[mask3 > .1] = 1
        mask3 = np.array(mask3, dtype=np.bool)
    
        extended_nrc = []
        for j, n in enumerate(nrc[mask3]):
            extended_nrc.append([[n] * mask2[mask3][j]])
        extended_nrc = utils.flatten_list(utils.flatten_list(extended_nrc))
        extended_nrc = np.array(extended_nrc)

        metallicity_at_infall = s.gas['metals_at_infall'][mask] / s.gas['mass_at_infall'][mask]
        metallicity_at_ejection = s.gas['metals_at_ejection'][mask] / s.gas['mass_at_ejection'][mask]
        metals_infall =  s.gas['metals_at_infall'][mask]
        metals_ejection = s.gas['metals_at_ejection'][mask]
        cycle_r_max = s.gas['cycle_r_max'][mask]
        ejection_z = pg.cosmology.a2z(s.gas['ejection_a'])[mask]
        
        if len(cycle_r_max) == 0:
            continue
        average_nor, _, _ = stats.binned_statistic(cycle_r_max,
            extended_nrc, statistic='median', bins=bins)
        average_z, _, _ = stats.binned_statistic(cycle_r_max,
            metallicity_at_ejection, statistic='median', bins=bins)
        
        ax[i].set_ylabel('z', color='b')
        ax[i].step(bins[:-1], average_z)
        ax[i].tick_params('y', colors='b')
        ax[i].set_title('z: %s-%s' % (bin_end, bin_start))

        ax[i].tick_params(labelbottom='off') 

        second_axes.append(ax[i].twinx())
        second_axes[-1].grid(False)
        second_axes[-1].step(bins[:-1], average_nor, color='r')
        second_axes[-1].set_ylabel('median num of rec', color='r')
        second_axes[-1].tick_params('y', colors='r')
        y_min.append(ax[i].get_ylim()[0])
        y_max.append(ax[i].get_ylim()[1])
        y2_min.append(second_axes[-1].get_ylim()[0])
        y2_max.append(second_axes[-1].get_ylim()[1])

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

    average_nor, _, _ = stats.binned_statistic(np.concatenate(cycle_r_max),
        extended_nor, statistic='median', bins=bins)
    average_z, _, _ = stats.binned_statistic(np.concatenate(cycle_r_max),
        np.concatenate(metals_ejection), statistic='median', bins=bins)

    ax[-2].set_title('All together')
    ax[-2].set_xlabel(('cycle r max [kpc]'))
    ax[-2].set_ylabel('z', color='b')
    ax[-2].step(bins[:-1], average_z)
    ax[-2].tick_params('y', colors='b')
    y_min.append(ax[-2].get_ylim()[0])
    y_max.append(ax[-2].get_ylim()[1])
    for a in ax[:-1]:
        a.set_ylim(np.min(y_min), np.max(y_max))
        a.set_xlim((0, 200))

    ax2 = ax[-2].twinx()
    ax2.grid(False)
    ax2.step(bins[:-1], average_nor, color='r')
    ax2.set_ylabel('median num of rec', color='r')
    ax2.tick_params('y', colors='r')
    y2_min.append(ax2.get_ylim()[0])
    y2_max.append(ax2.get_ylim()[1])
    for a in second_axes:
        a.set_ylim(np.min(y2_min), np.max(y2_max))
    ax2.set_ylim(np.min(y2_min), np.max(y2_max))
     
    ax[-1].set_xlim((0, 4))
    ax[-1].set_xlabel('redshift')
    ax[-1].set_ylabel('metallicity')
    ax[-1].step(edges_met[:-1], average_metallicity, label='average')
    ax[-1].step(edges_met[:-1], median_metallicity, label='median')
    ax[-1].legend(loc='upper right')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition))
        
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(7)
p.map(plot, utils.combinations)

