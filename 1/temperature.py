import pygad as pg
import matplotlib.pyplot as plt
import numpy as np

filename = __file__
type = ('disc-Uebler', 'disc', 'ball')[0]
s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_%s_070_470.dat' % (type))

T_ejection = [item[item > 0] for item in s.gas['T_at_ejection'][s.gas['num_recycled'] > 0]]
T_infall = [item[item > 0] for item in s.gas['T_at_infall'][s.gas['num_recycled'] > 0]]

for i, temp in enumerate(T_ejection):
    if len(temp) < len(T_infall[i]):
        T_infall[i] = T_infall[i][:-1]
    elif len(temp) > len(T_infall[i]):
        print 'hmm'

T_ejection = [item for sublist in T_ejection for item in sublist]
T_infall = [item for sublist in T_infall for item in sublist]

fig, ax = plt.subplots(2)
plt.tight_layout()
ax[0].set_xlabel("$T_{infall}$")
ax[0].set_ylabel("$T_{ejection}$")
ax[0].set_xlim((1e2, 1e9))
ax[0].set_ylim((1e2, 1e9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].scatter(T_infall, T_ejection)
ax[0].plot([0, 1e80], [0, 1e80], color='r')

ax[1].set_ylabel("count")
ax[1].set_xlabel("$T$ [K]")
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].hist(T_infall, bins=np.logspace(2, 9, 100/2), alpha=.5, label='infall')
ax[1].hist(T_ejection, bins=np.logspace(2, 9, 100/2), alpha=.5, label='ejection')
plt.legend(loc='best')

plt.savefig(filename.split("/")[-1][:-3] + '_' + type + ".png", bbox_inches='tight')

