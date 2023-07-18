import glob
from multiprocessing import Pool

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pygad as pg
import pygad.plotting
import utils
import utils_miha
from scipy import stats

filename = __file__

bins = np.arange(2, 8.1, 0.1)


def plot(args):
    halo = args[0]
    definition = args[1]
    modification = ""
    if "-" in halo:
        modification = halo[5:]
        halo = halo[:5]
    if definition.split("_")[-1] != definition:
        definition_filename = definition.split("_")[-1]
        definition = definition.split("_")[0]
    else:
        definition_filename = ""
    # print args
    path = "/ptmp/mpa/naab/REFINED/{}/SF_X/4x-2phase{}/out/snap_{}_4x_???".format(
        halo,
        modification,
        halo,
    )
    max = int(sorted(glob.glob(path))[-1][-3:])
    s_i, h_i, g_i = pg.prepare_zoom(
        "/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_%s"
        % (halo, modification, halo, max),
        gas_trace="/u/mihac/data/%s/4x-2phase%s/gastrace_%s"
        % (halo, modification, definition),
        star_form=None,
    )
    s_h, h_h, g_h = pg.prepare_zoom(
        "/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_%s"
        % (halo, modification, halo, max),
        gas_trace="/u/mihac/data/%s/4x-2phase%s/gastrace_%s"
        % (halo, modification, "halo"),
        star_form=None,
    )
    time_bins = np.arange(0, s_i.cosmic_time() + 0.1, 0.25)

    R200_frac, Tcrit, rhocrit = [0.15, "2e4 K", "1e-2 u/cm**3"]
    R200, M200 = pg.analysis.virial_info(s_i)
    ism = pg.BallMask(R200_frac * R200) & pg.ExprMask(
        '(temp < "{}") & (rho > "{}")'.format(Tcrit, rhocrit)
    )

    s_i.gas[pg.BallMask(R200)]["mass"].sum() - s_i.gas[ism]["mass"].sum()

    # print("%s %4.2f %s" % (halo, Mcgm / 1e11, Mcgm.units))

    cgm_i_ej = h_i.gas[~ism][~pg.BallMask(R200_frac * R200)]  # cgm gas
    mask_i = np.max(cgm_i_ej["T_at_ejection"], axis=-1) > 0  # ejected at least once
    cgm_i_m = cgm_i_ej[mask_i]  # cgm gas which was ejected at least once
    ism_IDs = s_i.gas["ID"][
        s_i.gas["num_recycled"] > -1
    ]  # IDs of ISM gas which was ejected at least once

    cgm_h_in = h_h.gas[~ism][~pg.BallMask(R200_frac * R200)][~pg.IDMask(ism_IDs)]
    mask_h = (
        cgm_h_in["infall_time"][:, 0]
        != cgm_h_in["infall_time"][cgm_h_in["infall_time"] > 0.1].flatten().min()
    )
    cgm_h_m = cgm_h_in[mask_h]
    mask_h2 = np.max(cgm_h_m["T_at_infall"], axis=1) > 0
    cgm_h_m = cgm_h_m[mask_h2]

    last_T_i = np.array(
        [x[x > 0][-1] for x in cgm_i_m["T_at_ejection"] if len(x[x > 0]) > 0]
    )
    last_mass_i = np.array(
        [x[x > 0][-1] for x in cgm_i_m["mass_at_ejection"] if len(x[x > 0]) > 0]
    )
    last_t_i = np.array(
        [x[x > 0][-1] for x in cgm_i_m["ejection_time"] if len(x[x > 0]) > 0]
    )

    last_T_h = np.array(
        [x[x > 0][-1] for x in cgm_h_m["T_at_infall"] if len(x[x > 0]) > 0]
    )
    last_mass_h = np.array(
        [x[x > 0][-1] for x in cgm_h_m["mass_at_infall"] if len(x[x > 0]) > 0]
    )
    last_t_h = np.array(
        [x[x > 0][-1] for x in cgm_h_m["infall_time"] if len(x[x > 0]) > 0]
    )

    mass_last_T_i, edges, _ = stats.binned_statistic(
        np.log10(last_T_i), cgm_i_m["mass"], bins=bins, statistic="sum"
    )
    mass_z0_i, _, _ = stats.binned_statistic(
        np.log10(cgm_i_m["temp"]), cgm_i_m["mass"], bins=bins, statistic="sum"
    )

    mass_last_T_h, edges, _ = stats.binned_statistic(
        np.log10(last_T_h), cgm_h_m["mass"], bins=bins, statistic="sum"
    )
    mass_z0_h, _, _ = stats.binned_statistic(
        np.log10(cgm_h_m["temp"]), cgm_h_m["mass"], bins=bins, statistic="sum"
    )

    last_mass_i_hist, _, _ = stats.binned_statistic(
        last_t_i, last_mass_i, bins=time_bins, statistic="sum"
    )
    last_mass_h_hist, _, _ = stats.binned_statistic(
        last_t_h, last_mass_h, bins=time_bins, statistic="sum"
    )

    # calculate the sum of all mass in the 4 quadrants for both
    tot_mass = h_i.gas[~ism][~pg.BallMask(R200_frac * R200)]["mass"].sum()
    u_l_i_m = (np.log10(cgm_i_m["temp"]) <= 5.2) & (np.log10(last_T_i) > 5.2)
    u_r_i_m = (np.log10(cgm_i_m["temp"]) > 5.2) & (np.log10(last_T_i) > 5.2)
    b_l_i_m = (np.log10(cgm_i_m["temp"]) <= 5.2) & (np.log10(last_T_i) <= 5.2)
    b_r_i_m = (np.log10(cgm_i_m["temp"]) > 5.2) & (np.log10(last_T_i) <= 5.2)
    u_l_i = cgm_i_m[u_l_i_m]["mass"].sum()
    u_r_i = cgm_i_m[u_r_i_m]["mass"].sum()
    b_l_i = cgm_i_m[b_l_i_m]["mass"].sum()
    b_r_i = cgm_i_m[b_r_i_m]["mass"].sum()

    u_l_h_m = (np.log10(cgm_h_m["temp"]) <= 5.2) & (np.log10(last_T_h) > 5.2)
    u_r_h_m = (np.log10(cgm_h_m["temp"]) > 5.2) & (np.log10(last_T_h) > 5.2)
    b_l_h_m = (np.log10(cgm_h_m["temp"]) <= 5.2) & (np.log10(last_T_h) <= 5.2)
    b_r_h_m = (np.log10(cgm_h_m["temp"]) > 5.2) & (np.log10(last_T_h) <= 5.2)
    u_l_h = cgm_h_m[u_l_h_m]["mass"].sum()
    u_r_h = cgm_h_m[u_r_h_m]["mass"].sum()
    b_l_h = cgm_h_m[b_l_h_m]["mass"].sum()
    b_r_h = cgm_h_m[b_r_h_m]["mass"].sum()
    print(
        "%s %4.2f %4.2f"
        % (
            halo,
            (u_l_h + u_r_h + b_l_h + b_r_h) * 100 / tot_mass,
            (u_l_i + u_r_i + b_l_i + b_r_i) * 100 / tot_mass,
        )
    )
    return

    """Plotting"""
    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map(
        "log10(temp)",
        np.log10(last_T_i),
        s=cgm_i_m,
        qty="mass",
        colors=np.log10(cgm_i_m["metallicity"] / pg.solar.Z()),
        logscale=True,
        bins=[bins, bins],
        extent=[[2, 8], [2, 8]],
        zero_is_white=True,
        clim=[-1.0, 0.49],
        colors_av="mass",
        clogscale=False,
        ax=ax1,
    )
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r"$log_{10} [ Z ] $", fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c="k")
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c="k")
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ ej}\ [K])$", fontsize=30)
    plt.setp(ax1.get_xticklabels(), fontsize=25)
    plt.setp(ax1.get_yticklabels(), fontsize=25)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += 0.05
    right -= 0.05
    (bottom, top) = (5.2 + 0.05, 5.2 - 0.05)
    ax1.text(
        left,
        bottom,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            u_l_i / tot_mass * 100, utils.str_fmt(u_l_i)
        ),
        fontsize=23,
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        bottom,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            u_r_i / tot_mass * 100, utils.str_fmt(u_r_i)
        ),
        fontsize=23,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        left,
        top,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            b_l_i / tot_mass * 100, utils.str_fmt(b_l_i)
        ),
        fontsize=23,
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        top,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            b_r_i / tot_mass * 100, utils.str_fmt(b_r_i)
        ),
        fontsize=23,
        ha="right",
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        5.2,
        5.2,
        r"$%.0f\%%$" % ((u_l_i + u_r_i + b_l_i + b_r_i) * 100 / tot_mass),
        va="center",
        ha="center",
        bbox=dict(facecolor="white"),
        fontsize=23,
    )

    ax2.set_xlabel(r"$Mass\ [10^8\ M_\odot]$")
    ax2.barh(edges[:-1], mass_last_T_i / 1e8, height=np.diff(edges), align="edge")
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft="off")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r"$Mass\ [10^8\ M_\odot]$")
    ax3.bar(edges[:-1], mass_z0_i / 1e8, width=np.diff(edges), align="edge")
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom="off")

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    # f.suptitle('%s%s - %s' % (halo, modification, 'ism ejected'), fontsize=44)

    plt.savefig(
        filename.split("/")[-1][:-3]
        + "_"
        + halo
        + modification
        + definition_filename
        + "_ism_ejection_metallicity.png",
        bbox_inches="tight",
    )

    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map(
        "log10(temp)",
        np.log10(last_T_h),
        s=cgm_h_m,
        qty="mass",
        colors=np.log10(cgm_h_m["metallicity"] / pg.solar.Z()),
        logscale=True,
        bins=[bins, bins],
        extent=[[2, 8], [2, 8]],
        zero_is_white=True,
        clim=[-1.0, 0.49],
        colors_av="mass",
        clogscale=False,
        ax=ax1,
    )
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel(r"$log_{10} [ Z ] $", fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c="k")
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c="k")
    ax1.grid(False)
    ax1.set_xlabel(r"$log_{10}(T_{z=0}\ [K])$", fontsize=30)
    ax1.set_ylabel(r"$log_{10}(T_{last\ in}\ [K])$", fontsize=30)
    plt.setp(ax1.get_xticklabels(), fontsize=25)
    plt.setp(ax1.get_yticklabels(), fontsize=25)
    ax1.set_xlim([2, 8])
    ax1.set_ylim([2, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += 0.05
    right -= 0.05
    (bottom, top) = (5.2 + 0.05, 5.2 - 0.05)
    ax1.text(
        left,
        bottom,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            u_l_h / tot_mass * 100, utils.str_fmt(u_l_h)
        ),
        fontsize=23,
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        bottom,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            u_r_h / tot_mass * 100, utils.str_fmt(u_r_h)
        ),
        fontsize=23,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        left,
        top,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            b_l_h / tot_mass * 100, utils.str_fmt(b_l_h)
        ),
        fontsize=23,
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        top,
        r"${:.0f}\%, \ {}\ M_\odot$".format(
            b_r_h / tot_mass * 100, utils.str_fmt(b_r_h)
        ),
        fontsize=23,
        ha="right",
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        5.2,
        5.2,
        r"$%.0f\%%$" % ((u_l_h + u_r_h + b_l_h + b_r_h) * 100 / tot_mass),
        va="center",
        ha="center",
        bbox=dict(facecolor="white"),
        fontsize=23,
    )

    ax2.set_xlabel(r"$Mass\ [10^8\ M_\odot]$")
    ax2.barh(edges[:-1], mass_last_T_h / 1e8, height=np.diff(edges), align="edge")
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft="off")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r"$Mass\ [10^8\ M_\odot]$")
    ax3.bar(edges[:-1], mass_z0_h / 1e8, width=np.diff(edges), align="edge")
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom="off")

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    # f.suptitle('%s%s - %s' % (halo, modification, 'halo accreted'), fontsize=44)
    plt.savefig(
        filename.split("/")[-1][:-3]
        + "_"
        + halo
        + modification
        + definition_filename
        + "_hallo_infall_metallicity.png",
        bbox_inches="tight",
    )

    # colors by last ej time
    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map(
        "log10(temp)",
        np.log10(last_T_i),
        s=cgm_i_m,
        qty="mass",
        colors=last_t_i,
        logscale=True,
        bins=[bins, bins],
        extent=[[2, 8], [2, 8]],
        zero_is_white=True,
        clim=[2, s_i.cosmic_time()],
        colors_av="mass",
        clogscale=False,
        ax=ax1,
    )
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel("last ISM ejection time [Gyr]", fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c="k")
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c="k")
    ax1.grid(False)
    ax1.set_xlabel(r"$\log_{10}(T_\mathrm{z=0}$ [K])", fontsize=30)
    ax1.set_ylabel(r"$\log_{10}(T_\mathrm{last\ ej}$ [K])", fontsize=30)
    plt.setp(ax1.get_xticklabels(), fontsize=25)
    plt.setp(ax1.get_yticklabels(), fontsize=25)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += 0.05
    right -= 0.05
    (bottom, top) = (5.2 + 0.05, 5.2 - 0.05)
    ax1.text(
        left,
        bottom,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (u_l_i / tot_mass * 100, utils.str_fmt(u_l_i)),
        fontsize=23,
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        bottom,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (u_r_i / tot_mass * 100, utils.str_fmt(u_r_i)),
        fontsize=23,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        left,
        top,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (b_l_i / tot_mass * 100, utils.str_fmt(b_l_i)),
        fontsize=23,
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        top,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (b_r_i / tot_mass * 100, utils.str_fmt(b_r_i)),
        fontsize=23,
        ha="right",
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        5.2,
        5.2,
        r"$%.0f\%%$" % ((u_l_i + u_r_i + b_l_i + b_r_i) * 100 / tot_mass),
        va="center",
        ha="center",
        bbox=dict(facecolor="white"),
        fontsize=23,
    )

    ax2.set_xlabel(r"$\mathrm{Mass\ [10^8\ M_\odot]}$")
    ax2.barh(edges[:-1], mass_last_T_i / 1e8, height=np.diff(edges), align="edge")
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft="off")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r"$\mathrm{Mass\ [10^8\ M_\odot]}$")
    ax3.bar(edges[:-1], mass_z0_i / 1e8, width=np.diff(edges), align="edge")
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom="off")

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    # ax3.text(3, 7.6, 'ISM ejected ', fontsize=44)

    plt.savefig(
        filename.split("/")[-1][:-3]
        + "_"
        + halo
        + modification
        + definition_filename
        + "_ism_ejection_ej_time.png",
        bbox_inches="tight",
    )

    f = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(4, 4)
    ax1 = plt.subplot(gs[:-1, 1:])
    ax2 = plt.subplot(gs[:-1, 0])
    ax3 = plt.subplot(gs[-1, 1:])

    _, _, _, cbar = pg.plotting.scatter_map(
        "log10(temp)",
        np.log10(last_T_h),
        s=cgm_h_m,
        qty="mass",
        colors=last_t_h,
        logscale=True,
        bins=[bins, bins],
        extent=[[2, 8], [2, 8]],
        zero_is_white=True,
        clim=[2, s_i.cosmic_time()],
        colors_av="mass",
        clogscale=False,
        ax=ax1,
    )
    cbar_ax = cbar.ax
    cbar_ax.set_xlabel("last halo infall time [Gyr]", fontsize=30)
    plt.setp(cbar_ax.get_xticklabels(), fontsize=14)
    ax1.plot([5.2, 5.2], [0, 1e100], lw=3, c="k")
    ax1.plot([0, 1e100], [5.2, 5.2], lw=3, c="k")
    ax1.grid(False)
    ax1.set_xlabel(r"$\log_{10}(T_\mathrm{z=0}$ [K])", fontsize=30)
    ax1.set_ylabel(r"$\log_{10}(T_\mathrm{last\ in}$ [K])", fontsize=30)
    plt.setp(ax1.get_xticklabels(), fontsize=25)
    plt.setp(ax1.get_yticklabels(), fontsize=25)
    ax1.set_xlim([3, 8])
    ax1.set_ylim([3, 8])
    limits = ax1.get_xlim()
    (left, right) = limits
    left += 0.05
    right -= 0.05
    (bottom, top) = (5.2 + 0.05, 5.2 - 0.05)
    ax1.text(
        left,
        bottom,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (u_l_h / tot_mass * 100, utils.str_fmt(u_l_h)),
        fontsize=23,
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        bottom,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (u_r_h / tot_mass * 100, utils.str_fmt(u_r_h)),
        fontsize=23,
        ha="right",
        va="bottom",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        left,
        top,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (b_l_h / tot_mass * 100, utils.str_fmt(b_l_h)),
        fontsize=23,
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        right,
        top,
        r"$%.0f\%%, \ %s\ \mathrm{M_\odot}$"
        % (b_r_h / tot_mass * 100, utils.str_fmt(b_r_h)),
        fontsize=23,
        ha="right",
        va="top",
        bbox=dict(facecolor="white"),
    )
    ax1.text(
        5.2,
        5.2,
        r"$%.0f\%%$" % ((u_l_h + u_r_h + b_l_h + b_r_h) * 100 / tot_mass),
        va="center",
        ha="center",
        bbox=dict(facecolor="white"),
        fontsize=23,
    )

    ax2.set_xlabel(r"Mass [$\mathrm{10^8\ M_\odot}$]")
    ax2.barh(edges[:-1], mass_last_T_h / 1e8, height=np.diff(edges), align="edge")
    ax2.set_xlim(ax2.get_xlim()[::-1])
    ax2.set_ylim(limits)
    ax2.tick_params(labelleft="off")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)

    ax3.set_ylabel(r"Mass [$\mathrm{10^8\ M_\odot}$]")
    ax3.bar(edges[:-1], mass_z0_h / 1e8, width=np.diff(edges), align="edge")
    ax3.set_ylim(ax3.get_ylim()[::-1])
    ax3.set_xlim(limits)
    ax3.tick_params(labelbottom="off")

    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    # f.suptitle('%s%s - %s' % (halo, modification, 'halo accreted'), fontsize=44)

    plt.savefig(
        filename.split("/")[-1][:-3]
        + "_"
        + halo
        + modification
        + definition_filename
        + "_hallo_infall_ej_time.png",
        bbox_inches="tight",
    )

    # histogram
    f, ax = plt.subplots(1, figsize=utils.figsize[::-1])
    ax.set_xlabel("Cosmic time [Gyr]")
    ax.set_ylabel(r"Accretion rate $\mathrm{[M_\odot\ yr^{-1}]}$")
    ax.set_xlim(time_bins[0], time_bins[-1] - 0.25 / 2)
    ax.plot([], [], " ", label="CGM gas")
    ax.bar(
        time_bins[:-1],
        last_mass_h_hist / 1e0 / (1e9 * np.diff(time_bins)),
        np.diff(time_bins),
        label="halo accreted",
    )
    ax.bar(
        time_bins[:-1],
        last_mass_i_hist / 1e0 / (1e9 * np.diff(time_bins)),
        np.diff(time_bins),
        label="ISM ejected",
    )
    ax.legend(loc="upper left", fontsize=20)
    f.tight_layout()
    plt.subplots_adjust(top=0.92)
    # ax.text(, 'cgm'), fontsize=44)

    plt.savefig(
        filename.split("/")[-1][:-3]
        + "_"
        + halo
        + modification
        + definition_filename
        + "_hist.png",
        bbox_inches="tight",
    )


p = Pool(1)
combinations = utils_miha.combinations
# combinations = [('M0858', 'ism')]
# combinations.insert(0, ('M1859-weakFB', 'ism'))
p.map(plot, combinations)
# plot(combinations[0])
