import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

    f = plt.figure(figsize=utils.figsize[::-1] * 1.3)
    gs = gridspec.GridSpec(3, 3)
    ax1 = plt.subplot(gs[:2, 1:])
    ax2 = plt.subplot(gs[:2, 0])
    ax3 = plt.subplot(gs[2, 1:])

    ax1.grid(True)
    ax1.set_xlim((1e3, 1e9))
    ax1.set_ylim((1e3, 1e9))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("T at z=0")
    ax1.set_ylabel("T at last ejection")
    ax1.scatter(cgm['temp'][np.max(cgm['T_at_ejection'], axis=-1) > 0], last_T, alpha=.2, edgecolor='none')

    ax3.set_ylabel('count')
    ax3.set_xlim((1e3, 1e9))
    ax3.set_xscale('log')
    ax3.bar(edges[:-1], x_h, width=np.diff(edges))
    ax3.set_ylim(ax3.get_ylim()[::-1])

    ax2.set_xlabel('count')
    ax2.set_yscale('log')
    ax2.barh(edges[:-1], y_h, height=np.diff(edges))
    ax2.set_xlim(ax2.get_xlim()[::-1])

    f.tight_layout()
    plt.subplots_adjust(top=0.88)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

