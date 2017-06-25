import pygad as pg
import matplotlib.pyplot as plt
import numpy as np


def time_in(s, array):
    array = np.array(array[array > 0])
    return pg.UnitArr(np.r_[array[1:], float(s.cosmic_time())] - array, 'Gyr')[:-1]


filename = __file__
type = ('disc-Uebler', 'disc', 'ball')[2]
s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))

data = [time_in(s, item) for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
data = [item for sublist in data for item in sublist]

fig, ax = plt.subplots(1)
ax.set_xlabel("infall time [Gyr]")
ax.set_ylabel("count")
ax.hist(data, bins=np.linspace(0, 1.5, 50))

plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".png", bbox_inches='tight')

