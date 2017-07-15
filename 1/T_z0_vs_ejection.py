import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import glob
from multiprocessing import Pool
import utils

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    R200_frac, Tcrit, rhocrit = [.15, '2e4 K', '1e-2 u/cm**3']
    R200, M200 = pg.analysis.virial_info(s)
    ism = pg.BallMask(R200_frac*R200) & \
                           pg.ExprMask('(temp < "%s") & (rho > "%s")' % (Tcrit, rhocrit))
    cgm = h.gas[~ism]
    last_T = np.array([x[x>0][-1] for x in cgm['T_at_ejection'] if len(x[x>0])>0])

    x_h, edges = np.histogram(last_T, bins=np.logspace(3, 9, 31))
    y_h, _ = np.histogram(cgm['temp'][np.max(cgm['T_at_ejection'], axis=-1) > 0], bins=np.logspace(3, 9, 31))

    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.grid(True)
    ax.set_xlim((1e3, 1e9))
    ax.set_ylim((1e3, 1e9))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("T at z=0")
    ax.set_ylabel("T at last ejection")
    ax.scatter(cgm['temp'][np.max(cgm['T_at_ejection'], axis=-1) > 0], last_T, alpha=.2, edgecolor='none')

    ax1 = ax.twinx()
    ax1.set_ylabel('count', color='r')
    ax1.tick_params('y', colors='r')
    ax1.bar(edges[:-1], x_h, width=np.diff(edges), alpha=.1, color='r')

    ax2 = ax.twiny()
    ax2.set_xlabel('count', color='g')
    ax2.tick_params('x', colors='g')
    ax2.barh(edges[:-1], y_h, height=np.diff(edges), alpha=.1, color='g')

    f.tight_layout()
    plt.subplots_adjust(top=0.88)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

