import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils

filename = __file__

type = ('disc-Uebler', 'disc', 'ball')[2]
s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))

metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
cycle_r_max = [item[item > 0] for item in s.gas['cycle_r_max'][s.gas['num_recycled'] > -1]]
number_of_recycles = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
extended_nor = []

for i, particle in enumerate(cycle_r_max):
    for line in particle:
        extended_nor.append(number_of_recycles[i])

average_nor, edges_nor, count_nor = stats.binned_statistic(np.concatenate(cycle_r_max), extended_nor, statistic='mean', bins=np.linspace(0, 40, 41))
average_z, edges, count = stats.binned_statistic(np.concatenate(cycle_r_max), np.concatenate(metals_ejection), statistic='sum', bins=np.linspace(0, 40, 41))

edges, average_z = utils.prepare_step(edges, average_z)
edges_nor, average_nor = utils.prepare_step(edges_nor, average_nor)

fig, ax = plt.subplots(1)
ax.grid(True)
ax.set_xlabel(('cycle r max [kpc]'))
ax.set_ylabel('z', color='b')
ax.step(edges, average_z)
ax.tick_params('y', colors='b')

ax2 = ax.twinx()
ax2.step(edges_nor, average_nor, color='r')
ax2.set_ylabel('average number of recycles', color='r')
ax2.tick_params('y', colors='r')

plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".png", bbox_inches='tight')
