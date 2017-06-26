import pygad as pg
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import utils

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

fig, ax = plt.subplots(1)
plt.tight_layout()
ax.grid(True)
ax.scatter(halo_z_1st_infall, disk_z_1st_infall)
ax.set_xlabel('z infall ball .1 R200')
ax.set_ylabel('z infall disk Uebler .1 R200')
ax.set_xlim((1e-4, 1e-1))
ax.set_ylim((1e-4, 1e-1))
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot([0, 1e10], [0, 1e10], c='r')

plt.savefig(filename.split("/")[-1][:-3] + ".png", bbox_inches='tight')

