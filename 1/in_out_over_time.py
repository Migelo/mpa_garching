import numpy as np
import matplotlib.pyplot as plt
import pygad as pg

filename = __file__

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/out/snap_M0977_4x_470', gas_trace='/ptmp/mpa/naab/REFINED/M0977/SF_X/4x-2phase/gastrace_M0977_4x_disc-Uebler_070_470.dat')

timerange = np.linspace(0, s.cosmic_time(), 26*4)
infall, out = np.zeros(len(timerange) - 1), np.zeros(len(timerange) - 1)

metals_at_infall = s.gas['metals_at_infall'][s.gas['num_recycled'] > -1]
metals_at_ejection = s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1]
infall_time = s.gas['infall_time'][s.gas['num_recycled'] > -1]
ejection_time = s.gas['ejection_time'][s.gas['num_recycled'] > -1]

for i, particle in enumerate(metals_at_infall):
    particle = np.c_[particle, infall_time[i], metals_at_ejection[i], ejection_time[i]]
    for metals_in, time_in, metals_out, time_out in particle:
        done = [False, False]
        for j in range(len(infall)):
            if timerange[j] < time_in <= timerange[j + 1]:
                infall[j] += metals_in
                done[0] = True
                if done[0] and done[1]:
                    break
            if timerange[j] < time_out <= timerange[j + 1]:
                out[j] += metals_out
                done[1] = True
                if done[0] and done[1]:
                    break

timerange = timerange.tolist()
timerange.insert(0, timerange[0])
timerange.append(timerange[-1])

infall = infall.tolist()
infall.insert(0, infall[0])
infall.insert(0, 0)
infall.append(0)

out = out.tolist()
out.insert(0, out[0])
out.insert(0, 0)
out.append(0)

fig, ax = plt.subplots(1)
ax.grid(True)
ax.set_xlabel("Time [Gyr]")
ax.set_ylabel("Mass [$M_{\odot}$]")

ax.step(timerange, infall, label='Infall')
ax.step(timerange, out, label='Ejection')
plt.legend(loc='best')

plt.savefig(filename.split("/")[-1][:-3] + ".pdf", bbox_inches='tight')

