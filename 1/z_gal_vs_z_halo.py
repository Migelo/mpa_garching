import pygad as pg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats
import utils
import glob
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from multiprocessing import Pool

filename = __file__

n, start, end = 50, -6, -1
bins=np.logspace(start, end, n)

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    
    su, hu, gu = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    sb, hb, gb = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, 'halo'), star_form=None)

    ball_IDs = set(sb.gas['ID'][sb.gas['num_recycled'] > -1])
    disk_IDs = set(su.gas['ID'][su.gas['num_recycled'] > -1])
    ID_overlapp = ball_IDs.intersection(disk_IDs)
    overlapp_mask = pg.IDMask(list(ID_overlapp))
    tracked_ball = sb[overlapp_mask]
    tracked_disk = su[overlapp_mask]

    halo_z = tracked_ball.gas['metals_at_infall'] / tracked_ball.gas['mass_at_infall'] 
    disk_z = tracked_disk.gas['metals_at_infall'] / tracked_disk.gas['mass_at_infall'] 

    halo_z_1st_infall = [x[0] for x in halo_z]
    disk_z_1st_infall = [x[0] for x in disk_z]

    mass_at_halo_infall, edges_h, count = stats.binned_statistic(halo_z_1st_infall,
        tracked_ball.gas['mass_at_infall'][:, 0], statistic='sum', bins=bins)
    mass_at_disc_infall, edges_d, count = stats.binned_statistic(disk_z_1st_infall,
        tracked_disk.gas['mass_at_infall'][:, 0], statistic='sum', bins=bins)

    disk_infall_avg, edges_h, count = stats.binned_statistic(halo_z_1st_infall,
        disk_z_1st_infall, statistic='mean', bins=bins)

    x, y = bins, bins
    X, Y = np.meshgrid(x[:-1], y[:-1])
    Z = np.histogram2d(disk_z_1st_infall, halo_z_1st_infall, bins=bins)[0]
   
    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    ax1.scatter(halo_z_1st_infall, disk_z_1st_infall, alpha=.1, edgecolor=None)
    ax1.set_xlabel(r'$z_{halo\ infall}$')
    ax1.set_ylabel(r'$z_{galaxy\ infall}$')
    ax1.set_xlim((1e-6, 1e-1))
    ax1.set_ylim((1e-4, 1e-1))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.plot([0, 1e10], [0, 1e10], c='r', label='1:1')
    cont = ax1.contour(X, Y, Z, np.linspace(0, np.percentile(Z, 99), 10))
    ax1.step(edges_h[:-1], disk_infall_avg, where='mid',
        color='k', lw=2.5, label='average')
    ax1.legend(loc='upper right')
   
    ax2.barh(edges_d[:-1], mass_at_disc_infall / 1e9,
        height=np.diff(edges_d), log=True, align='edge')
    ax2.set_yscale('log')
    ax2.set_xscale('linear')
    ax2_xlim = ax2.get_xlim()[::-1]
    ax2_xlim[1] = 0
    ax2.set_xlim(ax2_xlim)
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_xlabel(r'$Mass\ [10^{9}\ M_\odot]$')
    ax2.tick_params(labelleft='off')    


    ax3.bar(edges_h[:-1], mass_at_halo_infall / 1e9,
        width=np.diff(edges_h), align='edge', log=True)
    ax3.set_xscale('log')    
    ax3.set_yscale('linear')
    ax3_ylim = ax3.get_ylim()[::-1]
    ax3_ylim[1] = 0
    ax3.set_ylim(ax3_ylim)
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_ylabel(r'$Mass\ [10^{9}\ M_\odot]$')
    ax3.tick_params(labelbottom='off')    


    f.tight_layout()
    plt.subplots_adjust(top=0.93)
    f.suptitle('%s - %s' % (halo, definition), fontsize=32)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

