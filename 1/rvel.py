import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import pygad.plotting
import glob
import utils
import os
import pdb
from multiprocessing import Pool

filename = __file__

def plot(args):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    
    x_limits = []
    y_limits = []
    f, ax = plt.subplots(3, 2, figsize=(20, 12))
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="mass", ax=ax[0,0], logbin=False, ylog=False, label='recycled')
    ax[0,0].legend(loc='upper right')

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="HI", ax=ax[0,1], logbin=False, ylog=False, label='recycled')
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="MgII", ax=ax[1,0], logbin=False, ylog=False, label='recycled')

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="CIV", ax=ax[1,1], logbin=False, ylog=False, label='recycled')

    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="OVI", ax=ax[2,0], logbin=False, ylog=False, label='recycled')
    
    pg.plotting.profile(s.gas, '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], logbin=False, ylog=False, label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], logbin=False, ylog=False, label='non-recycled')
    pg.plotting.profile(s.gas[np.max(s.gas['T_at_ejection'] > 0, axis=1)], '200 kpc', 'vrad', dens=False, av="metallicity", ax=ax[2,1], logbin=False, ylog=False, label='recycled')

    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=20)
    
    for axes in ax.flatten():
        axes.grid(True)
        x_limits.append(axes.get_xlim())
        y_limits.append(axes.get_ylim())

    x_limits = np.array(x_limits)
    y_limits = np.array(y_limits)
    x_limits = np.min(x_limits[:, 0]), np.max(x_limits[:, 1])
    y_limits = np.min(y_limits[:, 0]), np.max(y_limits[:, 1])
    
    if args[2] == '0.0':
        return (x_limits, y_limits)
    for axes in ax.flatten():
        axes.set_xlim(args[2][0])
        axes.set_ylim(args[2][1])

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".pdf", bbox_inches='tight')
p = Pool(4)
limits = p.map(plot, np.column_stack((utils.combinations, np.zeros(len(utils.combinations)))))
x_lim = np.array([x[0] for x in limits])
y_lim = np.array([x[1] for x in limits])

limits = [[x_lim[:,0].min(), x_lim[:,1].max()], [y_lim[:,0].min(), y_lim[:,1].max()]]
arguments = []
for item in utils.combinations:
    item = list(item)
    arguments.append(item + [limits])
p.map(plot, arguments)

