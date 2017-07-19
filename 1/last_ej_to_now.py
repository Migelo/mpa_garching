import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import pygad as pg
import pygad.plotting
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
    
    # get IDs with at least 1 ejection
    mask = np.max(s.gas['T_at_ejection'], axis=-1) > 0
    last_ids = s.gas[mask]['ID']
    id_mask = pg.IDMask(last_ids)
    s_idx = np.argsort(last_ids)

    # iterate over all snapshots and get the cosmic time and
    # IDs of all gas particles
    T_id = []
    t = []
    snaps = sorted(glob.glob(path))
    start, step = 70, 1
    for snapshot in snaps[start::step]:
        s1, h1, g1 = pg.prepare_zoom(snapshot, gas_trace=None, star_form=None)
        t.append(s1.cosmic_time())
        idx = np.argsort(s1.gas[id_mask]['ID'])
        T_id.append(s1.gas[pg.IDMask(last_ids)]['temp'][idx])
    t = np.array(t)
    T_id = np.array(T_id).T

    last_t = np.array([x[x>0][-1] for x in s.gas[id_mask]['ejection_time'][s_idx]])
    to_plot = []
    for i, item in enumerate(T_id):
        to_plot.append(np.column_stack((t[t >= last_t[i]], item[t >= last_t[i]])))

    colors = np.log10(s.gas[id_mask]['metallicity'][s_idx]/pg.solar.Z())
    norm = mpl.colors.Normalize(vmin=-2, vmax=0.5)
    scalarMap = cm.ScalarMappable(norm=norm, cmap=pg.plotting.cm_my_viridis)
    
    f, ax = plt.subplots(2, figsize=utils.figsize[::1])
    for axis in ax:
        axis.grid(True)
        axis.set_ylim((1e2, 1e8))
        axis.set_yscale('log')
        axis.set_xlabel("time [Gyr]")
        axis.set_ylabel("T [K]")
    
        ax[0].set_title('All data')
    for i, item in enumerate(to_plot):
        colorVal = scalarMap.to_rgba(colors[i])
        ax[0].plot(item[:, 0], item[:, 1], color=colorVal)
        
        ax[1].set_title('Every 1000th line')
    for i, item in enumerate(to_plot[::1000]):
        colorVal = scalarMap.to_rgba(colors[1000*i])
        ax[1].plot(item[:, 0], item[:, 1], color=colorVal)

    f.tight_layout()
    plt.subplots_adjust(top=0.88)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)


