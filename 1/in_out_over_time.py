import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import glob
from multiprocessing import Pool
import utils

filename = __file__

def plot(args):
    halo, definition = args

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    timerange = np.linspace(0, s.cosmic_time(), 13 * 8)
    rec = s.gas[s.gas['num_recycled'] > -1]
    metals_infall = rec['metals_at_infall']
    metals_infall_initial = metals_infall[:, 0]
    metals_infall_reac = metals_infall[:, 1:]
    metals_ejection = rec['metals_at_ejection'] 
    metals_ejection_initial = metals_ejection[:, 0]
    metals_ejection_reac = metals_ejection[:, 1:]

    mass_infall = rec['mass_at_infall']
    mass_infall_initial = mass_infall[:, 0]
    mass_infall_reac = mass_infall[:, 1:]
    mass_ejection = rec['mass_at_ejection'] 
    mass_ejection_initial = mass_ejection[:, 0]
    mass_ejection_reac = mass_ejection[:, 1:]
    
    infall_time = rec['infall_time']
    infall_time_initial = infall_time[:, 0]
    infall_time_reac = infall_time[:, 1:]
    ejection_time = rec['ejection_time']
    ejection_time_initial = ejection_time[:, 0]
    ejection_time_reac = ejection_time[:, 1:]

    metals_infall, edges, count = stats.binned_statistic(np.concatenate(infall_time), np.concatenate(metals_infall), statistic='sum', bins=timerange)
    metals_ejection, edges, count = stats.binned_statistic(np.concatenate(ejection_time), np.concatenate(metals_ejection), statistic='sum', bins=timerange)
    metals_infall_initial, edges_init, count = stats.binned_statistic(infall_time_initial, metals_infall_initial, statistic='sum', bins=timerange)
    metals_ejection_initial, edges_init, count = stats.binned_statistic(ejection_time_initial, metals_ejection_initial, statistic='sum', bins=timerange)
    metals_infall_reac, edges_reac, count = stats.binned_statistic(np.concatenate(infall_time_reac), np.concatenate(metals_infall_reac), statistic='sum', bins=timerange)
    metals_ejection_reac, edges_reac, count = stats.binned_statistic(np.concatenate(ejection_time_reac), np.concatenate(metals_ejection_reac), statistic='sum', bins=timerange)
    mass_infall, edges, count = stats.binned_statistic(np.concatenate(infall_time), np.concatenate(mass_infall), statistic='sum', bins=timerange)
    mass_ejection, edges, count = stats.binned_statistic(np.concatenate(ejection_time), np.concatenate(mass_ejection), statistic='sum', bins=timerange)
    mass_infall_initial, edges_init, count = stats.binned_statistic(infall_time_initial, mass_infall_initial, statistic='sum', bins=timerange)
    mass_ejection_initial, edges_init, count = stats.binned_statistic(ejection_time_initial, mass_ejection_initial, statistic='sum', bins=timerange)
    mass_infall_reac, edges_reac, count = stats.binned_statistic(np.concatenate(infall_time_reac), np.concatenate(mass_infall_reac), statistic='sum', bins=timerange)
    mass_ejection_reac, edges_reac, count = stats.binned_statistic(np.concatenate(ejection_time_reac), np.concatenate(mass_ejection_reac), statistic='sum', bins=timerange)

    dt = timerange[1]*1e9

    metals_infall /=  dt
    metals_ejection /= dt
    metals_infall_initial /= dt
    metals_ejection_initial /= dt
    metals_infall_reac /= dt
    metals_ejection_reac /= dt
    mass_infall /= dt
    mass_ejection /= dt
    mass_infall_initial /= dt
    mass_ejection_initial /= dt
    mass_infall_reac /= dt
    mass_ejection_reac /= dt


    metals_infall = utils.prepare_step(metals_infall)
    metals_ejection = utils.prepare_step(metals_ejection)
    metals_infall_initial = utils.prepare_step(metals_infall_initial)
    metals_ejection_initial = utils.prepare_step(metals_ejection_initial)
    metals_infall_reac = utils.prepare_step(metals_infall_reac)
    metals_ejection_reac = utils.prepare_step(metals_ejection_reac)
    mass_infall = utils.prepare_step(mass_infall)
    mass_ejection = utils.prepare_step(mass_ejection)
    mass_infall_initial = utils.prepare_step(mass_infall_initial)
    mass_ejection_initial = utils.prepare_step(mass_ejection_initial)
    mass_infall_reac = utils.prepare_step(mass_infall_reac)
    mass_ejection_reac = utils.prepare_step(mass_ejection_reac)


    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(25, 20))
    plt.tight_layout()
    for ax in (ax1, ax2, ax3, ax4):
        ax.grid(True)
        ax.set_xlabel("Time [yr]")
        ax.set_xlim((0, s.cosmic_time()))

    ax1.set_ylabel("mass rate [$M_{\odot}/yr$]")
    ax2.set_ylabel("mass rate [$M_{\odot}/yr$]")
    ax3.set_ylabel("metal rate [$M_{\odot}/yr$]")
    ax4.set_ylabel("metal rate [$M_{\odot}/yr$]")

    ax1.set_title("Mass")
    ax1.step(edges, mass_infall, label='Total infall')
    ax1.step(edges_init, mass_infall_initial, label='First infall')
    ax1.step(edges_reac, mass_infall_reac, label='Subsequent infall')
    lgd1 = ax1.legend(loc='upper left')
    ax2.step(edges, mass_ejection, label='Total ejection')
    ax2.step(edges_init, mass_ejection_initial, label='First ejection')
    ax2.step(edges_reac, mass_ejection_reac, label='Subsequent ejection')
    lgd2 = ax2.legend(loc='upper left')

    ax3.set_title("Metals")
    ax3.step(edges, metals_infall, label='Total infall')
    ax3.step(edges_init, metals_infall_initial, label='First infall')
    ax3.step(edges_reac, metals_infall_reac, label='Subsequent ejection')
    lgd3 = ax3.legend(loc='upper left')
    ax4.step(edges, metals_ejection, label='Total ejection')
    ax4.step(edges_init, metals_ejection_initial, label='First ejection')
    ax4.step(edges_reac, metals_ejection_reac, label='Subsequent ejection')
    lgd4 = ax4.legend(loc='upper left')

    utils.mkdir('./%s' % (halo))
    plt.savefig('./%s/' % (halo) + ("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

