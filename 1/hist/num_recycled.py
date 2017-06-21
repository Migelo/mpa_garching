import pygad as pg
import matplotlib.pyplot as plt
import numpy as np

filename = __file__

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_disc-Uebler_070_470.dat')

x = range(max(s.gas['num_recycled'])+1)
y = [s.gas['mass'][s.gas['num_recycled']==i].sum() for i in range(max(s.gas['num_recycled'])+1)]

fig, ax = plt.subplots(1)
ax.set_xlabel("$N_{recycled}$")
ax.set_ylabel("$M_{\odot}$")
ax.bar(x, y, width=1)

plt.savefig(filename.split("/")[-1][:-3] + ".pdf")

