import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils
import glob
from multiprocessing import Pool

filename = __file__

start, stop, step = 0, 200, 5
bins=np.arange(start, stop + step, step)
time_bins = [4, 3, 2, 1, .5, 0]

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    
    time_at_r_max = pg.cosmology.a2z(s.gas['cycle_r_max_at'])
    nrc = s.gas['num_recycled']
        
    f, ax = plt.subplots(5, figsize=(8, 8/2*2**.5*5))

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
#        number_of_recycles = s.gas['num_recycled'][mask]
    
#        average_metallicity, _, _ = stats.binned_statistic(
#            ejection_z, metallicity_at_ejection,
#            bins=time_bins[::-1], statistic='mean')
#        median_metallicity, _, _ = stats.binned_statistic(
#            ejection_z, metallicity_at_ejection,
#            bins=time_bins[::-1], statistic='median')
        
        if len(cycle_r_max) == 0:
            continue
        average_nor, _, _ = stats.binned_statistic(cycle_r_max,
            extended_nrc, statistic='mean', bins=bins)
        average_z, _, _ = stats.binned_statistic(cycle_r_max,
            metallicity_at_ejection, statistic='mean', bins=bins)
        
        ax[i].set_ylabel('z', color='b')
        ax[i].step(bins[:-1], average_z)
        ax[i].tick_params('y', colors='b')
        ax[i].set_title('z: %s-%s' % (bin_end, bin_start))

        if item == time_bins[-1]:
            ax[i].set_xlabel(('cycle r max [kpc]'))
        else:
            ax[i].tick_params(labelbottom='off') 


        ax2 = ax[i].twinx()
        ax2.step(bins[:-1], average_nor, color='r')
        ax2.set_ylabel('average number of recycles', color='r')
        ax2.tick_params('y', colors='r')
        
#        ax[i].set_xlim((0, 4))
#        ax[i].set_xlabel('redshift')
#        ax[i].set_ylabel('metallicity')
#        ax[i].step(edges_met, utils.prepare_step(average_metallicity), label='average')
#        ax[i].step(edges_met, utils.prepare_step(median_metallicity), label='median')
#        ax[i].legend(loc='upper right')

    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
        
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)
#plot(('M1196', 'ism'))

#for item in utils.combinations:
#    plot(item)
