import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import glob
from multiprocessing import Pool
import utils

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M1196/SF_X/4x-2phase/out/snap_M1196_4x_470', gas_trace='/u/mihac/data/M1196/4x-2phase/gastrace_disc', star_form=None)

R200_frac, Tcrit, rhocrit = [.15, '2e4 K', '1e-2 u/cm**3']
R200, M200 = pg.analysis.virial_info(s)
ism = s.gas[pg.BallMask(R200_frac*R200) & \
                       pg.ExprMask('(temp < "%s") & (rho > "%s")' % (Tcrit, rhocrit))]
halo = s.gas[pg.BallMask(R200)]
cgm = s.gas[pg.IDMask(set(halo['ID']) - set(ism['ID']))]
mask = np.max(cgm['T_at_ejection'] > 0, axis=-1)
cgm_recycled = cgm[mask]
last_T = np.array([x[x > 0][-1] for x in cgm_recycled['T_at_ejection']])

fig, ax = plt.subplots(1)
ax.grid(True)
ax.set_xlim((1e3, 1e9))
ax.set_ylim((1e3, 1e9))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("T at z=0")
ax.set_ylabel("T at last ejection")
fig.tight_layout()
ax.scatter(cgm_recycled['temp'], last_T, alpha=.1, edgecolor='none')
#plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')
plt.savefig('T_z0_vs_ejection.png', bbox_inches='tight')

