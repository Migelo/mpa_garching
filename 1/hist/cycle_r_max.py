import pygad as pg
import matplotlib.pyplot as plt
import numpy as np

filename = __file__

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_disc-Uebler_070_470.dat')

cycle_r_max = [item[item > 0] for item in s.gas['cycle_r_max'][s.gas['num_recycled'] > 0]]
cycle_r_max = [item for sublist in cycle_r_max for item in sublist]

fig, ax = plt.subplots(1)
plt.tight_layout()
ax.set_xlabel("$r_{max}$ [kpc]")
ax.set_ylabel("count")
ax.hist(cycle_r_max, bins=np.linspace(0, 40, 30))
plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')
