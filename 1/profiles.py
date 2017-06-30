import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import utils

filename = __file__

for type in ('disc-Uebler', 'disc', 'ball', 'ism'):
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/u/mihac/data/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))

r = np.array([np.sqrt(x**2 + y**2 + z**2) for x, y, z in s.gas['pos'][s.gas['num_recycled'] > -1]])
nrc = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
z = s.gas['metals'][s.gas['num_recycled'] > -1]/s.gas['mass'][s.gas['num_recycled'] > -1]

nrc_avg, edges_nrc, count_nrc = stats.binned_statistic(r, nrc, statistic='mean', bins=np.linspace(0, 3000, 3000))
z_avg, edges_z, count_z = stats.binned_statistic(r, z, statistic='mean', bins=np.linspace(0, 3000, 3000))
edges_nrc, nrc_avg = utils.prepare_step(edges_nrc, nrc_avg)
edges_z, z_avg = utils.prepare_step(edges_z, z_avg)

fig, (ax1, ax2) = plt.subplots(2)
plt.tight_layout()
ax1.grid(True)
ax1.set_xlim((0, 100))
ax1.set_ylim((0, max(nrc[r < ax1.get_xlim()[1]])))
ax1.set_xlabel('r [kpc]')
ax1.set_ylabel('number of recycles')
ax1.scatter(r, nrc)
ax1.step(edges_nrc, nrc_avg, c='r')

ax2.grid(True)
ax2.set_xlim((0, 100))
ax2.set_ylim((5e-4, 1e-1))
ax2.set_yscale(('log'))
ax2.set_xlabel('r [kpc]')
ax2.set_ylabel('z')
ax2.scatter(r, z)
ax2.step(edges_z, z_avg, c='r')

plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".png", bbox_inches='tight')

