import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
    modification = ''
    if '-' in halo:
        modification = halo[5:]
        halo = halo[:5]
    if definition.split('_')[-1] != definition:
        definition_filename = definition.split('_')[-1]
        definition = definition.split('_')[0]
    else:
        definition_filename = ''
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_???' % (halo, modification, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    print halo, definition
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_%s' % (halo, modification, halo, max),
        gas_trace=None, star_form=None)
    
    f, ax = plt.subplots(2, figsize=(10, 5))
    f.suptitle(halo, fontsize=44)
    pg.plotting.image(s.gas, qty='HI', extent='306 kpc', ax=ax[0])
    pg.plotting.image(s.gas, qty='HI', extent='306 kpc', ax=ax[1], xaxis=0, yaxis=2)
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '.png', bbox_inches='tight')

p = Pool(4)
combinations = utils.combinations
#combinations = [('M0858', 'ism')]
p.map(plot, combinations)
