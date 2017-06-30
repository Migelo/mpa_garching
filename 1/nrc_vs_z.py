import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils

filename = __file__

for type in ('disc-Uebler', 'disc', 'ball', 'ism'):
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/u/mihac/data/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))


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

    fix, ax = plt.subplots(1)
    ax.grid(True)
    ax.set_xlabel("z")
    ax.set_ylabel("number of recycles")
    ax.set_xlim((0, .08))
#ax.set_ylim((0, 80))
    ax.scatter(most_recent_z, number_of_recycles, label='number of recycles')
    ax.plot(edges_nor[: -1], average_nor, c='r', label='average')
    plt.legend(loc='best')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".png", bbox_inches='tight')

