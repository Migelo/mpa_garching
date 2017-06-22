import pygad as pg
import matplotlib.pyplot as plt
import numpy as np

filename = __file__

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_disc-Uebler_070_470.dat')

metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > 0]]
metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > 0]]

for i, temp in enumerate(metals_ejection):
    if len(temp) < len(metals_infall[i]):
        metals_infall[i] = metals_infall[i][:-1]
    elif len(temp) > len(metals_infall[i]):
        print 'hmm'

metals_ejection = [item for sublist in metals_ejection for item in sublist]
metals_infall = [item for sublist in metals_infall for item in sublist]

fig, ax = plt.subplots(2)
plt.tight_layout()
ax[0].set_xlabel("$T_{infall}$")
ax[0].set_ylabel("$T_{ejection}$")
ax[0].set_xlim((-1e6, 2.1e7))
ax[0].set_ylim((-1e6, 2e7))
ax[0].scatter(s.gas['metals_at_infall'][s.gas['num_recycled'] > -1], s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1])

ax[1].set_ylabel("count")
ax[1].set_xlabel("z")
ax[1].hist(metals_infall, bins=np.linspace(0, 1.2e4, 100), alpha=.5, label='infall')
ax[1].hist(metals_ejection, bins=np.linspace(0, 1.2e4, 100), alpha=.5, label='ejection')
plt.legend(loc='best')
plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

