import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
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
    
    mask = np.max(s.gas['T_at_ejection'], axis=-1) > 0
    last_ids = s.gas[mask]['ID']

    T = []
    t = []
    id = []
    snaps = sorted(glob.glob(path))
    for snapshot in snaps[70::]:
        s1, h1, g1 = pg.prepare_zoom(snapshot, gas_trace=None, star_form=None)
        T.append(s1.gas['temp'])
        t.append(s1.cosmic_time())
        id.append(s1.gas['ID'])

    t = np.array(t)
    id_set = set.intersection(*map(set, id))
    id_set.intersection_update(last_ids)
    id_mask = pg.IDMask(id_set)
    
    T_id = []
    for snapshot in snaps[70::]:
        s2, h2, g2 = pg.prepare_zoom(snapshot, gas_trace=None, star_form=None)
        s_m = s2[id_mask]
        T_id.append(s_m.gas['temp'])
    
    T_id = np.array(T_id).T
    last_t = np.array([x[x>0][-1] for x in s.gas[id_mask]['ejection_time']])
    to_plot = []
    for i, item in enumerate(T_id):
        to_plot.append(np.column_stack(((item[t >= last_t[i]]), t[t >= last_t[i]])))

    colors = s.gas[id_mask]['metallicity']/pg.solar.Z()
    norm = mpl.colors.LogNorm(vmin=colors.min(), vmax=colors.max())
    scalarMap = cm.ScalarMappable(norm=norm, cmap='viridis')
    
#    Z = [[0,0],[0,0]]
#    levels = [scalarMap.to_rgba(x) for x in colors[:10]]
#    CS3 = plt.contourf(Z, levels, cmap='viridis')
#    plt.clf()
    
    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.grid(True)
    ax.set_ylim((1e2, 1e7))
    ax.set_yscale('log')
    ax.set_xlabel("time [Gyr]")
    ax.set_ylabel("T [K]")
    for i, item in enumerate(to_plot):
        colorVal = scalarMap.to_rgba(colors[i])
        ax.plot(item[:, 1], item[:, 0])

    f.tight_layout()
    plt.subplots_adjust(top=0.88)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)


