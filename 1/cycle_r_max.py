import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
import glob
from multiprocessing import Pool
import utils

filename = __file__
def plot(args):
    halo = args[0]
    type = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, type), star_form=None)

    cycle_r_max = [item[item > 0] for item in s.gas['cycle_r_max'][s.gas['num_recycled'] > -1]]
    particle_mass = [np.average(item[item > 0]) for item in s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    bins = np.logspace(-1, 6, 15)
    mass_in_bins = np.zeros(len(bins) - 1)

    for i, particle in enumerate(cycle_r_max):
        for r in particle:
            for j in range(len(mass_in_bins)):
                if bins[j] < r <= bins[j+1]:
                    mass_in_bins[j] += particle_mass[i]

    mass_in_bins /= 1e10

    y = mass_in_bins.tolist()
    x = bins.tolist()
    y.insert(0, y[0])
    x.insert(0, x[0])
    y.insert(0, 0)
    x.append(x[-1])
    y.append(0)


    cycle_r_max = [item for sublist in cycle_r_max for item in sublist]
    fig, (ax, ay) = plt.subplots(2)
    ax.set_xlabel("$r_{max}$ [kpc]")
    ax.set_ylabel("mass [$10^{10}\ M_{\odot}$]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
    ax.grid(True)
    ax.step(x, y)

    ay.set_xlabel("$r_{max}$ [kpc]")
    ay.set_ylabel("count")
    ay.grid(True)
    ay.hist(cycle_r_max, bins=np.linspace(min(cycle_r_max)*.9, np.median(cycle_r_max)*4, 25))
    plt.tight_layout()

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

p = Pool(4)
p.map(plot, utils.combinations)

