import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils
import matplotlib.cm as cm
import matplotlib.mlab as mlab

filename = __file__

su, hu, gu = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_disc-Uebler_070_470.dat')
sb, hb, gb = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_ball_070_470.dat')

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

count_h, edges_h, count = stats.binned_statistic(halo_z_1st_infall, range(len(halo_z_1st_infall)), statistic='count', bins=np.logspace(-4, -1, 30))
normalisation = .1/float(max(count_h))
count_h = np.array(count_h).astype(float) * normalisation 
count_d, edges_d, count = stats.binned_statistic(disk_z_1st_infall, range(len(halo_z_1st_infall)), statistic='count', bins=np.logspace(-4, -1, 30))
normalisation = .1/float(max(count_d))
count_d = np.array(count_d).astype(float) * normalisation 

disk_infall_avg, edges_h, count = stats.binned_statistic(halo_z_1st_infall, disk_z_1st_infall, statistic='mean', bins=np.logspace(-4, -1, 30))
# count_h = list(count_h)
# count_d = list(count_d)
# count_h.insert(0, count_h[0])
# count_d.insert(0, count_d[0])

n = 50
x = np.logspace(-4, -1, n)
y = np.logspace(-4, -1, n)
X, Y = np.meshgrid(x[:-1], y[:-1])
Z = np.histogram2d(disk_z_1st_infall, halo_z_1st_infall, bins=np.logspace(-4, -1, n))[0]

fig, ax = plt.subplots(1)
plt.tight_layout()
ax.set_aspect('equal', adjustable='box')
ax.grid(True)
ax.scatter(halo_z_1st_infall, disk_z_1st_infall)
ax.set_xlabel('z infall ball')
ax.set_ylabel('z infall disk Uebler')
ax.set_xlim((1e-4, 1e-1))
ax.set_ylim((1e-4, 1e-1))
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot([0, 1e10], [0, 1e10], c='r')
ax.bar(edges_h[:-1], count_h, width=np.diff(edges_h), log=True, alpha=.35, color='g')
ax.barh(edges_d[:-1], count_d, height=np.diff(edges_d), log=True, alpha=.35, color='r')
cont = ax.contour(X, Y, Z, np.logspace(0, 3.2, 10))
ax.step(edges_h[:-1], disk_infall_avg, where='mid')
print edges_d[np.argmax(count_d)], edges_h[np.argmax(count_h)]

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

