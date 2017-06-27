import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import utils

filename = __file__

type = ('disc-Uebler', 'disc', 'ball')[2]
s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))

timerange = np.linspace(0, s.cosmic_time(), 26*4)

metals_infall = s.gas['metals_at_infall'][s.gas['num_recycled'] > -1]
metals_infall_initial = [item[0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1]]
metals_infall_reac = [item[item > 0][1:] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]
metals_ejection = s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1]
metals_ejection_initial = [item[item > 0][0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 0]
metals_ejection_reac = [item[item > 0][1:] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]

mass_infall = s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]
mass_ejection = s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]
mass_ejection_initial = [item[item > 0][0] for item in s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 0]
mass_infall_initial = [item[0] for item in s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
mass_infall_reac = [item[item > 0][1:] for item in s.gas['mass_at_infall'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]
mass_ejection_reac = [item[item > 0][1:] for item in s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]

infall_time = s.gas['infall_time'][s.gas['num_recycled'] > -1]
ejection_time = s.gas['ejection_time'][s.gas['num_recycled'] > -1]
ejection_time_initial = [item[item > 0][0] for item in s.gas['ejection_time'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 0]
infall_time_initial = [item[0] for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
infall_time_reac = [item[item > 0][1:] for item in s.gas['infall_time'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]
ejection_time_reac = [item[item > 0][1:] for item in s.gas['ejection_time'][s.gas['num_recycled'] > -1] if len(item[item > 0]) > 1]

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

dt = timerange[1]

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


edges, metals_infall = utils.prepare_step(edges, metals_infall)
foo, metals_ejection = utils.prepare_step(timerange, metals_ejection)
edges_init, metals_infall_initial = utils.prepare_step(edges_init, metals_infall_initial)
foo, metals_ejection_initial = utils.prepare_step(timerange, metals_ejection_initial)
edges_reac, metals_infall_reac = utils.prepare_step(edges_reac, metals_infall_reac)
foo, metals_ejection_reac = utils.prepare_step(timerange, metals_ejection_reac)
foo, mass_infall = utils.prepare_step(timerange, mass_infall)
foo, mass_ejection = utils.prepare_step(timerange, mass_ejection)
foo, mass_infall_initial = utils.prepare_step(timerange, mass_infall_initial)
foo, mass_ejection_initial = utils.prepare_step(timerange, mass_ejection_initial)
foo, mass_infall_reac = utils.prepare_step(timerange, mass_infall_reac)
foo, mass_ejection_reac = utils.prepare_step(timerange, mass_ejection_reac)


fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
plt.tight_layout()
for ax in (ax1, ax2, ax3, ax4):
    ax.grid(True)
    ax.set_xlabel("Time [Gyr]")
    ax.set_xlim((0, s.cosmic_time()))

ax1.set_ylabel("mass rate [$M_{\odot}/Gyr$]")
ax2.set_ylabel("mass rate [$M_{\odot}/Gyr$]")
ax3.set_ylabel("metal rate [$M_{\odot}/Gyr$]")
ax4.set_ylabel("metal rate [$M_{\odot}/Gyr$]")

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

plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".pdf", bbox_inches='tight')

