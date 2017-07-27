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

bins = np.arange(3, 8.1, .1)

def plot(args):
    halo = args[0]
    definition = args[1]
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s_i, h_i, g_i = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)
    s_h, h_h, g_h = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, 'halo'), star_form=None)
    time_bins = np.arange(0, s_i.cosmic_time() + .1, .25)

    R200_frac, Tcrit, rhocrit = [.15, '2e4 K', '1e-2 u/cm**3']
    R200, M200 = pg.analysis.virial_info(s_i)
    ism = pg.BallMask(R200_frac*R200) & \
                           pg.ExprMask('(temp < "%s") & (rho > "%s")' % (Tcrit, rhocrit))

    cgm_i_ej = h_i.gas[~ism][~pg.BallMask(R200_frac*R200)]
    mask_i = np.max(cgm_i_ej['T_at_ejection'], axis=-1) > 0
    cgm_i_m = cgm_i_ej[mask_i]
    ism_IDs = s_i.gas['ID'][s_i.gas['num_recycled'] > -1]
    
    cgm_h_in = h_h.gas[~ism][~pg.BallMask(R200_frac*R200)][~pg.IDMask(ism_IDs)]
    mask_h = cgm_h_in['infall_time'][:, 0] != cgm_h_in['infall_time'][cgm_h_in['infall_time'] > .1].flatten().min()
    cgm_h_m = cgm_h_in[mask_h]
    mask_h2 = np.max(cgm_h_m['T_at_infall'], axis=1) > 0
    cgm_h_m = cgm_h_m[mask_h2]

    last_T_i = np.array([x[x>0][-1] for x in cgm_i_m['T_at_ejection'] if len(x[x>0])>0])
    last_mass_i = np.array([x[x>0][-1] for x in cgm_i_m['mass_at_ejection'] if len(x[x>0])>0])
    last_t_i = np.array([x[x>0][-1] for x in cgm_i_m['ejection_time'] if len(x[x>0])>0])

    last_T_h = np.array([x[x>0][-1] for x in cgm_h_m['T_at_infall'] if len(x[x>0])>0])
    last_mass_h = np.array([x[x>0][-1] for x in cgm_h_m['mass_at_infall'] if len(x[x>0])>0])
    last_t_h = np.array([x[x>0][-1] for x in cgm_h_m['infall_time'] if len(x[x>0])>0])

    mass_last_T_i, edges, _ = stats.binned_statistic(np.log10(last_T_i),
        cgm_i_m['mass'], bins=bins, statistic='sum')
    mass_z0_i, _, _ = stats.binned_statistic(np.log10(cgm_i_m['temp']),
        cgm_i_m['mass'], bins=bins, statistic='sum')

    mass_last_T_h, edges, _ = stats.binned_statistic(np.log10(last_T_h),
        cgm_h_m['mass'], bins=bins, statistic='sum')
    mass_z0_h, _, _ = stats.binned_statistic(np.log10(cgm_h_m['temp']),
        cgm_h_m['mass'], bins=bins, statistic='sum')

    last_mass_i_hist, _, _ = stats.binned_statistic(last_t_i, last_mass_i,
        bins=time_bins, statistic='sum')
    last_mass_h_hist, _, _ = stats.binned_statistic(last_t_h, last_mass_h,
        bins=time_bins, statistic='sum')
    
    # calculate the sum of all mass in the 4 quadrants for both
    tot_mass = h_i.gas[~ism][~pg.BallMask(R200_frac*R200)]['mass'].sum()
    u_l_i_m = (np.log10(cgm_i_m['temp']) <= 5.2) & (np.log10(last_T_i) > 5.2)
    u_r_i_m = (np.log10(cgm_i_m['temp']) > 5.2) & (np.log10(last_T_i) > 5.2) 
    b_l_i_m = (np.log10(cgm_i_m['temp']) <= 5.2) & (np.log10(last_T_i) <= 5.2) 
    b_r_i_m = (np.log10(cgm_i_m['temp']) > 5.2) & (np.log10(last_T_i) <= 5.2) 
    u_l_i = cgm_i_m[u_l_i_m]['mass'].sum()
    u_r_i = cgm_i_m[u_r_i_m]['mass'].sum()
    b_l_i = cgm_i_m[b_l_i_m]['mass'].sum()
    b_r_i = cgm_i_m[b_r_i_m]['mass'].sum()

    u_l_h_m = (np.log10(cgm_h_m['temp']) <= 5.2) & (np.log10(last_T_h) > 5.2)
    u_r_h_m = (np.log10(cgm_h_m['temp']) > 5.2) & (np.log10(last_T_h) > 5.2) 
    b_l_h_m = (np.log10(cgm_h_m['temp']) <= 5.2) & (np.log10(last_T_h) <= 5.2) 
    b_r_h_m = (np.log10(cgm_h_m['temp']) > 5.2) & (np.log10(last_T_h) <= 5.2) 
    u_l_h = cgm_h_m[u_l_h_m]['mass'].sum()
    u_r_h = cgm_h_m[u_r_h_m]['mass'].sum()
    b_l_h = cgm_h_m[b_l_h_m]['mass'].sum()
    b_r_h = cgm_h_m[b_r_h_m]['mass'].sum()
    
    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map('log10(temp)', np.log10(last_T_i),
        s=cgm_i_m, qty='mass', colors=np.log10(cgm_i_m['metallicity']/pg.solar.Z()),
        #s=cgm_m, qty='mass', colors=last_t,
        logscale=True, bins=[bins, bins], extent=[[3, 8], [3, 8]],
        zero_is_white=True, clim=[-1., .5], colors_av='mass', 
        #zero_is_white=True, clim=[2, s.cosmic_time()], colors_av='mass', 
        clogscale=False, ax=ax1)
#        clogscale=False, ax=ax1)
    cbar_ax = cbar.ax
    #cbar_ax.set_xlabel('last ejection time [Gyr]', fontsize=30)
    cbar_ax.set_xlabel(r'$log_{10} [ Z ] $', fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c='k')
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c='k')
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ ej}\ [K])$", fontsize=30)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += .05
    right -= .05
    (bottom, top) = (5.2 + .05, 5.2 - .05)
    ax1.text(left, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_l_i/tot_mass * 100, utils.str_fmt(u_l_i)), fontsize=23,
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(right, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_r_i/tot_mass * 100, utils.str_fmt(u_r_i)), fontsize=23, ha='right',
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(left, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_l_i/tot_mass * 100, utils.str_fmt(b_l_i)), fontsize=23,
             va='top', bbox=dict(facecolor='white'))
    ax1.text(right, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_r_i/tot_mass * 100, utils.str_fmt(b_r_i)), fontsize=23, ha='right',
             va='top', bbox=dict(facecolor='white'))
    ax1.text(5.2, 5.2, r'$%.0f\%%$' % ((u_l_i + u_r_i + b_l_i + b_r_i) * 100 / tot_mass),
        va='center', ha='center', bbox=dict(facecolor='white'), fontsize=23)


    ax2.set_xlabel(r'$Mass\ [10^8\ M_\odot]$')
    ax2.barh(edges[:-1], mass_last_T_i / 1e8, height=np.diff(edges), align='edge')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft='off')    
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r'$Mass\ [10^8\ M_\odot]$')
    ax3.bar(edges[:-1], mass_z0_i / 1e8, width=np.diff(edges), align='edge')
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom='off')

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, 'ism ejected'), fontsize=44)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_ism_ejection_metallicity.png', bbox_inches='tight')





    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map('log10(temp)', np.log10(last_T_h),
        s=cgm_h_m, qty='mass', colors=np.log10(cgm_h_m['metallicity']/pg.solar.Z()),
        #s=cgm_i_m, qty='mass', colors=last_t,
        logscale=True, bins=[bins, bins], extent=[[3, 8], [3, 8]],
        zero_is_white=True, clim=[-1., .5], colors_av='mass', 
        #zero_is_white=True, clim=[2, s.cosmic_time()], colors_av='mass', 
        clogscale=False, ax=ax1)
#        clogscale=False, ax=ax1)
    cbar_ax = cbar.ax
    #cbar_ax.set_xlabel('last ejection time [Gyr]', fontsize=30)
    cbar_ax.set_xlabel(r'$log_{10} [ Z ] $', fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c='k')
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c='k')
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ infall}\ [K])$", fontsize=30)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += .05
    right -= .05
    (bottom, top) = (5.2 + .05, 5.2 - .05)
    ax1.text(left, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_l_h/tot_mass * 100, utils.str_fmt(u_l_h)), fontsize=23,
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(right, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_r_h/tot_mass * 100, utils.str_fmt(u_r_h)), fontsize=23, ha='right',
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(left, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_l_h/tot_mass * 100, utils.str_fmt(b_l_h)), fontsize=23,
             va='top', bbox=dict(facecolor='white'))
    ax1.text(right, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_r_h/tot_mass * 100, utils.str_fmt(b_r_h)), fontsize=23, ha='right',
             va='top', bbox=dict(facecolor='white'))
    ax1.text(5.2, 5.2, r'$%.0f\%%$' % ((u_l_h + u_r_h + b_l_h + b_r_h) * 100 / tot_mass),
        va='center', ha='center', bbox=dict(facecolor='white'), fontsize=23)


    ax2.set_xlabel(r'$Mass\ [10^8\ M_\odot]$')
    ax2.barh(edges[:-1], mass_last_T_h / 1e8, height=np.diff(edges), align='edge')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft='off')    
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r'$Mass\ [10^8\ M_\odot]$')
    ax3.bar(edges[:-1], mass_z0_h / 1e8, width=np.diff(edges), align='edge')
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom='off')

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, 'halo accreted'), fontsize=44)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_hallo_infall_metallicity.png', bbox_inches='tight')






    # colors by last ej time
    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map('log10(temp)', np.log10(last_T_i),
        s=cgm_i_m, qty='mass', colors=last_t_i,
        logscale=True, bins=[bins, bins], extent=[[3, 8], [3, 8]],
        zero_is_white=True, clim=[2, s_i.cosmic_time()], colors_av='mass', 
        clogscale=False, ax=ax1)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel('last ejection time [Gyr]', fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c='k')
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c='k')
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ ej}\ [K])$", fontsize=30)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += .05
    right -= .05
    (bottom, top) = (5.2 + .05, 5.2 - .05)
    ax1.text(left, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_l_i/tot_mass * 100, utils.str_fmt(u_l_i)), fontsize=23,
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(right, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_r_i/tot_mass * 100, utils.str_fmt(u_r_i)), fontsize=23, ha='right',
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(left, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_l_i/tot_mass * 100, utils.str_fmt(b_l_i)), fontsize=23,
             va='top', bbox=dict(facecolor='white'))
    ax1.text(right, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_r_i/tot_mass * 100, utils.str_fmt(b_r_i)), fontsize=23, ha='right',
             va='top', bbox=dict(facecolor='white'))
    ax1.text(5.2, 5.2, r'$%.0f\%%$' % ((u_l_i + u_r_i + b_l_i + b_r_i) * 100 / tot_mass),
        va='center', ha='center', bbox=dict(facecolor='white'), fontsize=23)


    ax2.set_xlabel(r'$Mass\ [10^8\ M_\odot]$')
    ax2.barh(edges[:-1], mass_last_T_i / 1e8, height=np.diff(edges), align='edge')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft='off')    
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r'$Mass\ [10^8\ M_\odot]$')
    ax3.bar(edges[:-1], mass_z0_i / 1e8, width=np.diff(edges), align='edge')
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom='off')

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, 'ism ejected'), fontsize=44)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_ism_ejection_ej_time.png', bbox_inches='tight')





    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map('log10(temp)', np.log10(last_T_h),
        s=cgm_h_m, qty='mass', colors=last_t_h,
        logscale=True, bins=[bins, bins], extent=[[3, 8], [3, 8]],
        zero_is_white=True, clim=[2, s_i.cosmic_time()], colors_av='mass', 
        clogscale=False, ax=ax1)
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel('last infall time [Gyr]', fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c='k')
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c='k')
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ infall}\ [K])$", fontsize=30)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += .05
    right -= .05
    (bottom, top) = (5.2 + .05, 5.2 - .05)
    ax1.text(left, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_l_h/tot_mass * 100, utils.str_fmt(u_l_h)), fontsize=23,
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(right, bottom, r'$%.0f\%%, \ %s\ M_\odot$' % (u_r_h/tot_mass * 100, utils.str_fmt(u_r_h)), fontsize=23, ha='right',
             va='bottom', bbox=dict(facecolor='white'))
    ax1.text(left, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_l_h/tot_mass * 100, utils.str_fmt(b_l_h)), fontsize=23,
             va='top', bbox=dict(facecolor='white'))
    ax1.text(right, top, r'$%.0f\%%, \ %s\ M_\odot$' % (b_r_h/tot_mass * 100, utils.str_fmt(b_r_h)), fontsize=23, ha='right',
             va='top', bbox=dict(facecolor='white'))
    ax1.text(5.2, 5.2, r'$%.0f\%%$' % ((u_l_h + u_r_h + b_l_h + b_r_h) * 100 / tot_mass),
        va='center', ha='center', bbox=dict(facecolor='white'), fontsize=23)


    ax2.set_xlabel(r'$Mass\ [10^8\ M_\odot]$')
    ax2.barh(edges[:-1], mass_last_T_h / 1e8, height=np.diff(edges), align='edge')
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft='off')    
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r'$Mass\ [10^8\ M_\odot]$')
    ax3.bar(edges[:-1], mass_z0_h / 1e8, width=np.diff(edges), align='edge')
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom='off')

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, 'halo accreted'), fontsize=44)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_hallo_infall_ej_time.png', bbox_inches='tight')





    

    # histogram
    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.set_xlabel('Cosmic time [Gyr]')
    ax.set_ylabel(r'Mass $[10^9M_\odot]$')
    ax.set_xlim((0, 14))
    ax.bar(time_bins[:-1], last_mass_i_hist / 1e9, np.diff(time_bins),
        label='ism ejected')
    ax.bar(time_bins[:-1], last_mass_h_hist / 1e9, np.diff(time_bins),
        label='halo accreted', alpha=.5)
    ax.legend(loc='upper right')
    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    f.suptitle('%s - %s' % (halo, definition), fontsize=44)
    

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_hist.png', bbox_inches='tight')
    










p = Pool(8)
p.map(plot, utils.combinations)

