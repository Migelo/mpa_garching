import glob

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pygad as pg
import pygad.plotting
import utils

filename = __file__

bins = np.arange(0, 205, 5)

profiles = np.zeros(9, dtype=object)
HI_profiles = np.zeros(3, dtype=object)
MgII_profiles = np.zeros(3, dtype=object)
OVI_profiles = np.zeros(3, dtype=object)
profiles_surface = np.zeros(9, dtype=object)
HI_profiles_surface = np.zeros(3, dtype=object)
MgII_profiles_surface = np.zeros(3, dtype=object)
OVI_profiles_surface = np.zeros(3, dtype=object)

for i, args in enumerate(utils.combinations):
    halo = args[0]
    definition = "ism"

    path = "/ptmp/mpa/naab/REFINED/{}/SF_X/4x-2phase/out/snap_{}_4x_???".format(
        halo, halo
    )
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom(
        "/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s"
        % (halo, halo, max),
        gas_trace="/u/mihac/data/{}/4x-2phase/gastrace_{}".format(halo, definition),
        star_form=None,
    )

    m_nrec = s.gas["num_recycled"] == -1
    m_rec = np.max(s.gas["T_at_ejection"] > 0, axis=1)
    hot = pg.ExprMask("temp > '10**5.2 K'")

    if i == 0:
        profiles[0] = pg.analysis.profile_dens(s.gas, "mass", r_edges=bins)
        profiles[1] = pg.analysis.profile_dens(s.gas[hot], "mass", r_edges=bins)
        profiles[2] = pg.analysis.profile_dens(s.gas[~hot], "mass", r_edges=bins)
        profiles[3] = pg.analysis.profile_dens(s.gas[m_rec], "mass", r_edges=bins)
        profiles[4] = pg.analysis.profile_dens(s.gas[m_rec][hot], "mass", r_edges=bins)
        profiles[5] = pg.analysis.profile_dens(s.gas[m_rec][~hot], "mass", r_edges=bins)
        profiles[6] = pg.analysis.profile_dens(s.gas[m_nrec], "mass", r_edges=bins)
        profiles[7] = pg.analysis.profile_dens(s.gas[m_nrec][hot], "mass", r_edges=bins)
        profiles[8] = pg.analysis.profile_dens(
            s.gas[m_nrec][~hot], "mass", r_edges=bins
        )
        HI_profiles[0] = pg.analysis.profile_dens(s.gas, "HI", r_edges=bins)
        HI_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], "HI", r_edges=bins)
        HI_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], "HI", r_edges=bins)
        MgII_profiles[0] = pg.analysis.profile_dens(s.gas, "MgII", r_edges=bins)
        MgII_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], "MgII", r_edges=bins)
        MgII_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], "MgII", r_edges=bins)
        OVI_profiles[0] = pg.analysis.profile_dens(s.gas, "OVI", r_edges=bins)
        OVI_profiles[1] = pg.analysis.profile_dens(s.gas[m_rec], "OVI", r_edges=bins)
        OVI_profiles[2] = pg.analysis.profile_dens(s.gas[m_nrec], "OVI", r_edges=bins)
        profiles_surface[0] = pg.analysis.profile_dens(
            s.gas, "mass", r_edges=bins, proj=1
        )
        profiles_surface[1] = pg.analysis.profile_dens(
            s.gas[hot], "mass", r_edges=bins, proj=1
        )
        profiles_surface[2] = pg.analysis.profile_dens(
            s.gas[~hot], "mass", r_edges=bins, proj=1
        )
        profiles_surface[3] = pg.analysis.profile_dens(
            s.gas[m_rec], "mass", r_edges=bins, proj=1
        )
        profiles_surface[4] = pg.analysis.profile_dens(
            s.gas[m_rec][hot], "mass", r_edges=bins, proj=1
        )
        profiles_surface[5] = pg.analysis.profile_dens(
            s.gas[m_rec][~hot], "mass", r_edges=bins, proj=1
        )
        profiles_surface[6] = pg.analysis.profile_dens(
            s.gas[m_nrec], "mass", r_edges=bins, proj=1
        )
        profiles_surface[7] = pg.analysis.profile_dens(
            s.gas[m_nrec][hot], "mass", r_edges=bins, proj=1
        )
        profiles_surface[8] = pg.analysis.profile_dens(
            s.gas[m_nrec][~hot], "mass", r_edges=bins, proj=1
        )
        HI_profiles_surface[0] = pg.analysis.profile_dens(
            s.gas, "HI", r_edges=bins, proj=1
        )
        HI_profiles_surface[1] = pg.analysis.profile_dens(
            s.gas[m_rec], "HI", r_edges=bins, proj=1
        )
        HI_profiles_surface[2] = pg.analysis.profile_dens(
            s.gas[m_nrec], "HI", r_edges=bins, proj=1
        )
        MgII_profiles_surface[0] = pg.analysis.profile_dens(
            s.gas, "MgII", r_edges=bins, proj=1
        )
        MgII_profiles_surface[1] = pg.analysis.profile_dens(
            s.gas[m_rec], "MgII", r_edges=bins, proj=1
        )
        MgII_profiles_surface[2] = pg.analysis.profile_dens(
            s.gas[m_nrec], "MgII", r_edges=bins, proj=1
        )
        OVI_profiles_surface[0] = pg.analysis.profile_dens(
            s.gas, "OVI", r_edges=bins, proj=1
        )
        OVI_profiles_surface[1] = pg.analysis.profile_dens(
            s.gas[m_rec], "OVI", r_edges=bins, proj=1
        )
        OVI_profiles_surface[2] = pg.analysis.profile_dens(
            s.gas[m_nrec], "OVI", r_edges=bins, proj=1
        )
    else:
        profiles[0] = np.column_stack(
            (profiles[0], pg.analysis.profile_dens(s.gas, "mass", r_edges=bins))
        )
        profiles[1] = np.column_stack(
            (profiles[1], pg.analysis.profile_dens(s.gas[hot], "mass", r_edges=bins))
        )
        profiles[2] = np.column_stack(
            (profiles[2], pg.analysis.profile_dens(s.gas[~hot], "mass", r_edges=bins))
        )
        profiles[3] = np.column_stack(
            (profiles[3], pg.analysis.profile_dens(s.gas[m_rec], "mass", r_edges=bins))
        )
        profiles[4] = np.column_stack(
            (
                profiles[4],
                pg.analysis.profile_dens(s.gas[m_rec][hot], "mass", r_edges=bins),
            )
        )
        profiles[5] = np.column_stack(
            (
                profiles[5],
                pg.analysis.profile_dens(s.gas[m_rec][~hot], "mass", r_edges=bins),
            )
        )
        profiles[6] = np.column_stack(
            (profiles[6], pg.analysis.profile_dens(s.gas[m_nrec], "mass", r_edges=bins))
        )
        profiles[7] = np.column_stack(
            (
                profiles[7],
                pg.analysis.profile_dens(s.gas[m_nrec][hot], "mass", r_edges=bins),
            )
        )
        profiles[8] = np.column_stack(
            (
                profiles[8],
                pg.analysis.profile_dens(s.gas[m_nrec][~hot], "mass", r_edges=bins),
            )
        )
        HI_profiles[0] = np.column_stack(
            (HI_profiles[0], pg.analysis.profile_dens(s.gas, "HI", r_edges=bins))
        )
        HI_profiles[1] = np.column_stack(
            (HI_profiles[1], pg.analysis.profile_dens(s.gas[m_rec], "HI", r_edges=bins))
        )
        HI_profiles[2] = np.column_stack(
            (
                HI_profiles[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "HI", r_edges=bins),
            )
        )
        MgII_profiles[0] = np.column_stack(
            (MgII_profiles[0], pg.analysis.profile_dens(s.gas, "MgII", r_edges=bins))
        )
        MgII_profiles[1] = np.column_stack(
            (
                MgII_profiles[1],
                pg.analysis.profile_dens(s.gas[m_rec], "MgII", r_edges=bins),
            )
        )
        MgII_profiles[2] = np.column_stack(
            (
                MgII_profiles[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "MgII", r_edges=bins),
            )
        )
        OVI_profiles[0] = np.column_stack(
            (OVI_profiles[0], pg.analysis.profile_dens(s.gas, "OVI", r_edges=bins))
        )
        OVI_profiles[1] = np.column_stack(
            (
                OVI_profiles[1],
                pg.analysis.profile_dens(s.gas[m_rec], "OVI", r_edges=bins),
            )
        )
        OVI_profiles[2] = np.column_stack(
            (
                OVI_profiles[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "OVI", r_edges=bins),
            )
        )
        profiles_surface[0] = np.column_stack(
            (
                profiles_surface[0],
                pg.analysis.profile_dens(s.gas, "mass", r_edges=bins, proj=1),
            )
        )
        profiles_surface[1] = np.column_stack(
            (
                profiles_surface[1],
                pg.analysis.profile_dens(s.gas[hot], "mass", r_edges=bins, proj=1),
            )
        )
        profiles_surface[2] = np.column_stack(
            (
                profiles_surface[2],
                pg.analysis.profile_dens(s.gas[~hot], "mass", r_edges=bins, proj=1),
            )
        )
        profiles_surface[3] = np.column_stack(
            (
                profiles_surface[3],
                pg.analysis.profile_dens(s.gas[m_rec], "mass", r_edges=bins, proj=1),
            )
        )
        profiles_surface[4] = np.column_stack(
            (
                profiles_surface[4],
                pg.analysis.profile_dens(
                    s.gas[m_rec][hot], "mass", r_edges=bins, proj=1
                ),
            )
        )
        profiles_surface[5] = np.column_stack(
            (
                profiles_surface[5],
                pg.analysis.profile_dens(
                    s.gas[m_rec][~hot], "mass", r_edges=bins, proj=1
                ),
            )
        )
        profiles_surface[6] = np.column_stack(
            (
                profiles_surface[6],
                pg.analysis.profile_dens(s.gas[m_nrec], "mass", r_edges=bins, proj=1),
            )
        )
        profiles_surface[7] = np.column_stack(
            (
                profiles_surface[7],
                pg.analysis.profile_dens(
                    s.gas[m_nrec][hot], "mass", r_edges=bins, proj=1
                ),
            )
        )
        profiles_surface[8] = np.column_stack(
            (
                profiles_surface[8],
                pg.analysis.profile_dens(
                    s.gas[m_nrec][~hot], "mass", r_edges=bins, proj=1
                ),
            )
        )
        HI_profiles_surface[0] = np.column_stack(
            (
                HI_profiles_surface[0],
                pg.analysis.profile_dens(s.gas, "HI", r_edges=bins, proj=1),
            )
        )
        HI_profiles_surface[1] = np.column_stack(
            (
                HI_profiles_surface[1],
                pg.analysis.profile_dens(s.gas[m_rec], "HI", r_edges=bins, proj=1),
            )
        )
        HI_profiles_surface[2] = np.column_stack(
            (
                HI_profiles_surface[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "HI", r_edges=bins, proj=1),
            )
        )
        MgII_profiles_surface[0] = np.column_stack(
            (
                MgII_profiles_surface[0],
                pg.analysis.profile_dens(s.gas, "MgII", r_edges=bins, proj=1),
            )
        )
        MgII_profiles_surface[1] = np.column_stack(
            (
                MgII_profiles_surface[1],
                pg.analysis.profile_dens(s.gas[m_rec], "MgII", r_edges=bins, proj=1),
            )
        )
        MgII_profiles_surface[2] = np.column_stack(
            (
                MgII_profiles_surface[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "MgII", r_edges=bins, proj=1),
            )
        )
        OVI_profiles_surface[0] = np.column_stack(
            (
                OVI_profiles_surface[0],
                pg.analysis.profile_dens(s.gas, "OVI", r_edges=bins, proj=1),
            )
        )
        OVI_profiles_surface[1] = np.column_stack(
            (
                OVI_profiles_surface[1],
                pg.analysis.profile_dens(s.gas[m_rec], "OVI", r_edges=bins, proj=1),
            )
        )
        OVI_profiles_surface[2] = np.column_stack(
            (
                OVI_profiles_surface[2],
                pg.analysis.profile_dens(s.gas[m_nrec], "OVI", r_edges=bins, proj=1),
            )
        )

y_label_3d = pg.analysis.profile_dens(s.gas, "mass", r_edges=bins).units.latex()
y_label_2d = pg.analysis.profile_dens(s.gas, "mass", r_edges=bins, proj=1).units.latex()

for i, item in enumerate(profiles):
    for j in range(item.shape[1]):
        profiles[i][:, j] = np.log10(profiles[i][:, j])
        profiles_surface[i][:, j] = np.log10(profiles_surface[i][:, j])
    if i < 3:
        HI_profiles_surface[i][:, j] = np.log10(HI_profiles_surface[i][:, j])
        MgII_profiles_surface[i][:, j] = np.log10(MgII_profiles_surface[i][:, j])
        OVI_profiles_surface[i][:, j] = np.log10(OVI_profiles_surface[i][:, j])
        HI_profiles[i][:, j] = np.log10(HI_profiles[i][:, j])
        MgII_profiles[i][:, j] = np.log10(MgII_profiles[i][:, j])
        OVI_profiles[i][:, j] = np.log10(OVI_profiles[i][:, j])


# 3D density plots
f = plt.figure(figsize=utils.figsize[::-1] * 2)
gs = gridspec.GridSpec(3, 4)
ax1 = plt.subplot(gs[:-1, 0])
ax2 = plt.subplot(gs[-1, 0])
ax3 = plt.subplot(gs[:-1, 1])
ax4 = plt.subplot(gs[-1, 1])
ax5 = plt.subplot(gs[:-1, 2])
ax6 = plt.subplot(gs[-1, 2])
ax7 = plt.subplot(gs[:-1, 3])
ax8 = plt.subplot(gs[-1, 3])

ax1.set_ylabel(r"$log_{10}\left(\rho\ \mathrm{[%s]}\right)$" % y_label_3d, fontsize=30)
ax2.set_ylabel(r"$\rho_{non\ recycled}/\rho_{recycled}$", fontsize=28)

for a in (ax1, ax3, ax5, ax7):
    a.tick_params(labelbottom="off")
    a.set_xlim((0, 200))
for a in (ax2, ax4, ax6, ax8):
    a.set_ylim((1e-1, 1e1))
    a.set_yscale("log")
    a.plot([-200, 200], [1, 1], color="k")
    a.set_xlabel(r"$r\ [kpc]$", fontsize=44)

ax1.set_ylim((-1, 7))
ax3.set_ylim((-3, 5))
ax5.set_ylim((-6, 2))
ax7.set_ylim((-5, 3))

ax1.set_title("All gas")
ax3.set_title("HI")
ax5.set_title("MgII")
ax7.set_title("OVI")

ax1.plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color="b", label="total")
ax1.plot(
    bins[:-1],
    np.percentile(profiles[1], 50, axis=1),
    color="b",
    ls="-.",
    lw=2,
    label="total T>10^5.2 K",
)
ax1.plot(
    bins[:-1],
    np.percentile(profiles[2], 50, axis=1),
    color="b",
    ls="dashed",
    label="total T<=10^5.2 K",
)
ax1.plot(bins[:-1], np.percentile(profiles[3], 50, axis=1), color="r", label="recycled")
ax1.plot(
    bins[:-1], np.percentile(profiles[6], 50, axis=1), color="g", label="non recycled"
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles[0], 25, axis=1),
    np.percentile(profiles[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles[3], 25, axis=1),
    np.percentile(profiles[3], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles[6], 25, axis=1),
    np.percentile(profiles[6], 75, axis=1),
    alpha=0.25,
    color="g",
)
ax1.legend(loc="upper right", fontsize=20)

ax3.plot(bins[:-1], np.percentile(HI_profiles[0], 50, axis=1), color="b", label="total")
ax3.plot(
    bins[:-1], np.percentile(HI_profiles[1], 50, axis=1), color="r", label="recycled"
)
ax3.plot(
    bins[:-1],
    np.percentile(HI_profiles[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles[0], 25, axis=1),
    np.percentile(HI_profiles[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles[1], 25, axis=1),
    np.percentile(HI_profiles[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles[2], 25, axis=1),
    np.percentile(HI_profiles[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax5.plot(
    bins[:-1], np.percentile(MgII_profiles[0], 50, axis=1), color="b", label="total"
)
ax5.plot(
    bins[:-1], np.percentile(MgII_profiles[1], 50, axis=1), color="r", label="recycled"
)
ax5.plot(
    bins[:-1],
    np.percentile(MgII_profiles[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles[0], 25, axis=1),
    np.percentile(MgII_profiles[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles[1], 25, axis=1),
    np.percentile(MgII_profiles[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles[2], 25, axis=1),
    np.percentile(MgII_profiles[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax7.plot(
    bins[:-1], np.percentile(OVI_profiles[0], 50, axis=1), color="b", label="total"
)
ax7.plot(
    bins[:-1], np.percentile(OVI_profiles[1], 50, axis=1), color="r", label="recycled"
)
ax7.plot(
    bins[:-1],
    np.percentile(OVI_profiles[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles[0], 25, axis=1),
    np.percentile(OVI_profiles[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles[1], 25, axis=1),
    np.percentile(OVI_profiles[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles[2], 25, axis=1),
    np.percentile(OVI_profiles[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax2.plot(
    bins[:-1],
    10 ** np.percentile(profiles[6], 50, axis=1)
    / 10 ** np.percentile(profiles[3], 50, axis=1),
)
ax4.plot(
    bins[:-1],
    10 ** np.percentile(HI_profiles[2], 50, axis=1)
    / 10 ** np.percentile(HI_profiles[1], 50, axis=1),
)
ax6.plot(
    bins[:-1],
    10 ** np.percentile(MgII_profiles[2], 50, axis=1)
    / 10 ** np.percentile(MgII_profiles[1], 50, axis=1),
)
ax8.plot(
    bins[:-1],
    10 ** np.percentile(OVI_profiles[2], 50, axis=1)
    / 10 ** np.percentile(OVI_profiles[1], 50, axis=1),
)

ax2.set_xlim(ax1.get_xlim())
ax4.set_xlim(ax3.get_xlim())
ax6.set_xlim(ax5.get_xlim())
ax8.set_xlim(ax7.get_xlim())

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches="tight")

f, ax = plt.subplots(3, figsize=utils.figsize, dpi=220)
for a in ax:
    a.set_ylabel(
        r"$\log_{10}\left(\rho\ \mathrm{[%s]}\right)$" % y_label_3d, fontsize=20
    )
    a.set_ylim((1, 6))
    a.set_xlim((0, 200))
# ax[0].set_title('All gas')
ax[2].set_xlabel(r"$r\ \mathrm{[kpc]}$", fontsize=20)
ax[0].tick_params(labelbottom="off")
ax[1].tick_params(labelbottom="off")

ax[0].plot(bins[:-1], np.percentile(profiles[0], 50, axis=1), color="b", label="total")
ax[0].plot(
    bins[:-1],
    np.percentile(profiles[1], 50, axis=1),
    color="b",
    ls="-.",
    lw=2,
    label=r"total $T>10^{5.2}$ K",
)
ax[0].plot(
    bins[:-1],
    np.percentile(profiles[2], 50, axis=1),
    color="b",
    ls="dashed",
    label=r"total $T\leq10^{5.2}$ K",
)
ax[0].fill_between(
    bins[:-1],
    np.percentile(profiles[0], 25, axis=1),
    np.percentile(profiles[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax[1].plot(
    bins[:-1], np.percentile(profiles[3], 50, axis=1), color="r", label="recycled"
)
ax[1].plot(
    bins[:-1],
    np.percentile(profiles[4], 50, axis=1),
    color="r",
    ls="-.",
    lw=2,
    label=r"recycled $T>10^{5.2}$ K",
)
ax[1].plot(
    bins[:-1],
    np.percentile(profiles[5], 50, axis=1),
    color="r",
    ls="dashed",
    label=r"recycled $T\leq10^{5.2}$ K",
)
ax[1].fill_between(
    bins[:-1],
    np.percentile(profiles[3], 25, axis=1),
    np.percentile(profiles[3], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax[2].plot(
    bins[:-1], np.percentile(profiles[6], 50, axis=1), color="g", label="non recycled"
)
ax[2].plot(
    bins[:-1],
    np.percentile(profiles[7], 50, axis=1),
    color="g",
    ls="-.",
    lw=2,
    label=r"non recycled $T>10^{5.2}$ K",
)
ax[2].plot(
    bins[:-1],
    np.percentile(profiles[8], 50, axis=1),
    color="g",
    ls="dashed",
    label=r"non recycled $T\leq10^{5.2}$ K",
)
ax[2].fill_between(
    bins[:-1],
    np.percentile(profiles[6], 25, axis=1),
    np.percentile(profiles[6], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax[0].legend(loc="upper right", fontsize=15)
ax[1].legend(loc="upper right", fontsize=15)
ax[2].legend(loc="upper right", fontsize=15)

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + "_all_gas_only.png", bbox_inches="tight")

# surface plots
f = plt.figure(figsize=utils.figsize[::-1] * 2)
gs = gridspec.GridSpec(3, 4)
ax1 = plt.subplot(gs[:-1, 0])
ax2 = plt.subplot(gs[-1, 0])
ax3 = plt.subplot(gs[:-1, 1])
ax4 = plt.subplot(gs[-1, 1])
ax5 = plt.subplot(gs[:-1, 2])
ax6 = plt.subplot(gs[-1, 2])
ax7 = plt.subplot(gs[:-1, 3])
ax8 = plt.subplot(gs[-1, 3])

ax1.set_ylabel(
    r"$\log_{10}\left(\Sigma\ \mathrm{[%s]}\right)$" % y_label_2d, fontsize=30
)
ax2.set_ylabel(r"$\Sigma_\mathrm{non recycled}/\Sigma_\mathrm{recycled}$", fontsize=28)

for a in (ax1, ax3, ax5, ax7):
    a.tick_params(labelbottom="off")
    a.set_xlim((0, 200))
#    a.set_yscale('log')
for a in (ax2, ax4, ax6, ax8):
    #    a.set_ylim((1e-1, 1e1))
    a.set_ylim((1e-1, 1e1))
    a.set_xlim((0, 200))
    a.set_yscale("log")
    a.plot([0, 200], [1, 1], color="k")
    a.set_xlabel(r"$r\ \mathrm{[kpc]}$", fontsize=44)

ax1.set_ylim((3, 11))
ax3.set_ylim((-1, 7))
ax5.set_ylim((-3, 5))
ax7.set_ylim((-3, 5))

ax1.set_title("All gas")
ax3.set_title("HI")
ax5.set_title("MgII")
ax7.set_title("OVI")

ax1.plot(
    bins[:-1], np.percentile(profiles_surface[0], 50, axis=1), color="b", label="total"
)
ax1.plot(
    bins[:-1],
    np.percentile(profiles_surface[1], 50, axis=1),
    color="b",
    ls="-.",
    lw=2,
    label="total T>10^5.2 K",
)
ax1.plot(
    bins[:-1],
    np.percentile(profiles_surface[2], 50, axis=1),
    color="b",
    ls="dashed",
    label="total T<=10^5.2 K",
)
ax1.plot(
    bins[:-1],
    np.percentile(profiles_surface[3], 50, axis=1),
    color="r",
    label="recycled",
)
ax1.plot(
    bins[:-1],
    np.percentile(profiles_surface[6], 50, axis=1),
    color="g",
    label="non recycled",
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles_surface[0], 25, axis=1),
    np.percentile(profiles_surface[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles_surface[3], 25, axis=1),
    np.percentile(profiles_surface[3], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax1.fill_between(
    bins[:-1],
    np.percentile(profiles_surface[6], 25, axis=1),
    np.percentile(profiles_surface[6], 75, axis=1),
    alpha=0.25,
    color="g",
)
ax1.legend(loc="upper right", fontsize=20)

ax3.plot(
    bins[:-1],
    np.percentile(HI_profiles_surface[0], 50, axis=1),
    color="b",
    label="total",
)
ax3.plot(
    bins[:-1],
    np.percentile(HI_profiles_surface[1], 50, axis=1),
    color="r",
    label="recycled",
)
ax3.plot(
    bins[:-1],
    np.percentile(HI_profiles_surface[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles_surface[0], 25, axis=1),
    np.percentile(HI_profiles_surface[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles_surface[1], 25, axis=1),
    np.percentile(HI_profiles_surface[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax3.fill_between(
    bins[:-1],
    np.percentile(HI_profiles_surface[2], 25, axis=1),
    np.percentile(HI_profiles_surface[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax5.plot(
    bins[:-1],
    np.percentile(MgII_profiles_surface[0], 50, axis=1),
    color="b",
    label="total",
)
ax5.plot(
    bins[:-1],
    np.percentile(MgII_profiles_surface[1], 50, axis=1),
    color="r",
    label="recycled",
)
ax5.plot(
    bins[:-1],
    np.percentile(MgII_profiles_surface[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles_surface[0], 25, axis=1),
    np.percentile(MgII_profiles_surface[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles_surface[1], 25, axis=1),
    np.percentile(MgII_profiles_surface[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax5.fill_between(
    bins[:-1],
    np.percentile(MgII_profiles_surface[2], 25, axis=1),
    np.percentile(MgII_profiles_surface[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax7.plot(
    bins[:-1],
    np.percentile(OVI_profiles_surface[0], 50, axis=1),
    color="b",
    label="total",
)
ax7.plot(
    bins[:-1],
    np.percentile(OVI_profiles_surface[1], 50, axis=1),
    color="r",
    label="recycled",
)
ax7.plot(
    bins[:-1],
    np.percentile(OVI_profiles_surface[2], 50, axis=1),
    color="g",
    label="non recycled",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles_surface[0], 25, axis=1),
    np.percentile(OVI_profiles_surface[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles_surface[1], 25, axis=1),
    np.percentile(OVI_profiles_surface[1], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax7.fill_between(
    bins[:-1],
    np.percentile(OVI_profiles_surface[2], 25, axis=1),
    np.percentile(OVI_profiles_surface[2], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax2.plot(
    bins[:-1],
    10 ** np.percentile(profiles_surface[6], 50, axis=1)
    / 10 ** np.percentile(profiles_surface[3], 50, axis=1),
)
ax4.plot(
    bins[:-1],
    10 ** np.percentile(HI_profiles_surface[2], 50, axis=1)
    / 10 ** np.percentile(HI_profiles_surface[1], 50, axis=1),
)
ax6.plot(
    bins[:-1],
    10 ** np.percentile(MgII_profiles_surface[2], 50, axis=1)
    / 10 ** np.percentile(MgII_profiles_surface[1], 50, axis=1),
)
ax8.plot(
    bins[:-1],
    10 ** np.percentile(OVI_profiles_surface[2], 50, axis=1)
    / 10 ** np.percentile(OVI_profiles_surface[1], 50, axis=1),
)

for a in (ax1, ax3, ax5, ax7):
    plt.setp(a.yaxis.get_majorticklabels(), fontsize=20)
for a in (ax2, ax4, ax6, ax8):
    plt.setp(a.yaxis.get_majorticklabels(), fontsize=20)
    plt.setp(a.xaxis.get_majorticklabels(), fontsize=20)

f.tight_layout()

plt.savefig(filename.split("/")[-1][:-3] + "_surface.png", bbox_inches="tight")

f, ax = plt.subplots(3, figsize=utils.figsize)
for a in ax:
    a.set_ylabel(
        r"$log_{10}\left(\Sigma\ \mathrm{[%s]}\right)$" % y_label_2d, fontsize=20
    )
    #    a.set_yscale('log')
    #    a.set_ylim((1e4, 1e8))
    a.set_ylim((4, 8))
    a.set_xlim((0, 200))
ax[0].set_title("All gas")
ax[2].set_xlabel(r"$r\ \mathrm{[kpc]}$", fontsize=20)
ax[0].tick_params(labelbottom="off")
ax[1].tick_params(labelbottom="off")


ax[0].plot(
    bins[:-1], np.percentile(profiles_surface[0], 50, axis=1), color="b", label="total"
)
ax[0].plot(
    bins[:-1],
    np.percentile(profiles_surface[1], 50, axis=1),
    color="b",
    ls="-.",
    lw=2,
    label=r"total $T>10^{5.2}$ K",
)
ax[0].plot(
    bins[:-1],
    np.percentile(profiles_surface[2], 50, axis=1),
    color="b",
    ls="dashed",
    label=r"total $T\leq10^{5.2}$ K",
)
ax[0].fill_between(
    bins[:-1],
    np.percentile(profiles_surface[0], 25, axis=1),
    np.percentile(profiles_surface[0], 75, axis=1),
    alpha=0.25,
    color="b",
)
ax[1].plot(
    bins[:-1],
    np.percentile(profiles_surface[3], 50, axis=1),
    color="r",
    label="recycled",
)
ax[1].plot(
    bins[:-1],
    np.percentile(profiles_surface[4], 50, axis=1),
    color="r",
    ls="-.",
    lw=2,
    label=r"recycled $T>10^{5.2}$ K",
)
ax[1].plot(
    bins[:-1],
    np.percentile(profiles_surface[5], 50, axis=1),
    color="r",
    ls="dashed",
    label=r"recycled $T\leq10^{5.2}$ K",
)
ax[1].fill_between(
    bins[:-1],
    np.percentile(profiles_surface[3], 25, axis=1),
    np.percentile(profiles_surface[3], 75, axis=1),
    alpha=0.25,
    color="r",
)
ax[2].plot(
    bins[:-1],
    np.percentile(profiles_surface[6], 50, axis=1),
    color="g",
    label="non recycled",
)
ax[2].plot(
    bins[:-1],
    np.percentile(profiles_surface[7], 50, axis=1),
    color="g",
    ls="-.",
    lw=2,
    label=r"non recycled $T>10^{5.2}$ K",
)
ax[2].plot(
    bins[:-1],
    np.percentile(profiles_surface[8], 50, axis=1),
    color="g",
    ls="dashed",
    label=r"non recycled $T\leq10^{5.2}$ K",
)
ax[2].fill_between(
    bins[:-1],
    np.percentile(profiles_surface[6], 25, axis=1),
    np.percentile(profiles_surface[6], 75, axis=1),
    alpha=0.25,
    color="g",
)

ax[0].legend(loc="upper right", fontsize=15)
ax[1].legend(loc="upper right", fontsize=15)
ax[2].legend(loc="upper right", fontsize=15)

f.tight_layout()

plt.savefig(
    filename.split("/")[-1][:-3] + "_surface_all_gas_only.png", bbox_inches="tight"
)
