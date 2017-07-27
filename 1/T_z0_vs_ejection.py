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
    cgm_m = cgm[mask]
    last_T = np.array([x[x>0][-1] for x in cgm['T_at_ejection'] if len(x[x>0])>0])
    last_t = np.array([x[x>0][-1] for x in cgm['ejection_time'] if len(x[x>0])>0])

    bins = np.arange(3, 8.1, .1)
    mass_last_T, edges, _ = stats.binned_statistic(np.log10(last_T),
        cgm_m['mass'], bins=bins, statistic='sum')
    mass_z0, _, _ = stats.binned_statistic(np.log10(cgm_m['temp']),
        cgm_m['mass'], bins=bins, statistic='sum')

    # calculate the sum of all mass in the 4 quadrants we have
    u_l_m = (np.log10(cgm_m['temp']) <= 5.2) & (np.log10(last_T) > 5.2) #upper left mask
    u_r_m = (np.log10(cgm_m['temp']) > 5.2) & (np.log10(last_T) > 5.2) 
    b_l_m = (np.log10(cgm_m['temp']) <= 5.2) & (np.log10(last_T) <= 5.2) 
    b_r_m = (np.log10(cgm_m['temp']) > 5.2) & (np.log10(last_T) <= 5.2) 
    u_l = cgm_m[u_l_m]['mass'].sum()
    u_r = cgm_m[u_r_m]['mass'].sum()
    b_l = cgm_m[b_l_m]['mass'].sum()
    b_r = cgm_m[b_r_m]['mass'].sum()
    tot_mass = cgm_m['mass'].sum()

    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map('log10(temp)', np.log10(last_T),
        #s=cgm_m, qty='mass', colors=cgm_m['metallicity']/pg.solar.Z(),
        s=cgm_m, qty='mass', colors=last_t,
        logscale=True, bins=[bins, bins], extent=[[3, 8], [3, 8]],
#        zero_is_white=True, clim=[10**-1., 10**.5], colors_av='mass', 
        zero_is_white=True, clim=[2, s.cosmic_time()], colors_av='mass', 
#        clogscale=True, ax=ax1)
        clogscale=False, ax=ax1)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel('last ejection time [Gyr]', fontsize=30)
    #cbar_ax.set_xlabel(r'$log_{10} [ Z ] $', fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c='k')
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c='k')
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$")
    ax1.set_ylabel(r"$log_{10}(T_{last\ ej}\ [K])$")
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += .05
    right -= .05
    (bottom, top) = (5.2 + .05, 5.2 - .05)
    ax1.text(left, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_l/tot_mass * 100, utils.str_fmt(u_l)), fontsize=23,
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(right, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_r/tot_mass * 100, utils.str_fmt(u_r)), fontsize=23, ha='right',
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(left, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_l/tot_mass * 100, utils.str_fmt(b_l)), fontsize=23,
             va='top', bbox=dict(facecolor='white'))
    ax1.text(right, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_r/tot_mass * 100, utils.str_fmt(b_r)), fontsize=23, ha='right',
             va='top', bbox=dict(facecolor='white'))


    ax2.set_xlabel(r'$Mass\ [10^8\ M_\odot]$')
    ax2.barh(edges[:-1], mass_last_T / 1e8, height=np.diff(edges), align='edge')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft='off')    
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r'$Mass\ [10^8\ M_\odot]$')
    ax3.bar(edges[:-1], mass_z0 / 1e8, width=np.diff(edges), align='edge')
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom='off')

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, definition), fontsize=44)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

p = Pool(8)
p.map(plot, utils.combinations)

