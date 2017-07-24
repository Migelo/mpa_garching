import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pygad as pg
import pygad.plotting
from scipy import stats
import glob
from multiprocessing import Pool
import utils

plt.style.use('general')

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
    mask = np.max(cgm['T_at_ejection'], axis=-1) > 0
    last_T = np.array([x[x>0][-1] for x in cgm['T_at_ejection'] if len(x[x>0])>0])
    
    bins = np.arange(3, 8.1, .1)
    x_h, edges, _ = stats.binned_statistic(np.log10(last_T),
        np.log10(cgm[mask]['mass']), bins=bins, statistic='sum')
    y_h, _, _ = stats.binned_statistic(np.log10(cgm['temp'][mask]),
        np.log10(cgm[mask]['mass']), bins=bins, statistic='sum')

    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])
    
    pg.plotting.scatter_map(np.log10(cgm['temp'][mask]),
        np.log10(last_T), s=cgm[mask], qty='mass',
        logscale=True, ax=ax1, bins=[bins, bins],
        extent=[[3, 8], [3, 8]])
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0})$")
    ax1.set_ylabel(r"$log_{10}(T_{last\ ej})$")
    limits = ax1.get_xlim()

    ax3.set_ylabel(r'$log_{10}(mass\ [M_\odot])$')
    ax3.bar(edges[:-1], x_h, width=np.diff(edges))
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)

    ax2.set_xlabel(r'$log_{10}(mass\ [M_\odot])$')
    ax2.barh(edges[:-1], y_h, height=np.diff(edges))
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    f.tight_layout()
    plt.subplots_adjust(top=0.95)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

