"""Incomplete script!"""
import pygad as pg
import matplotlib.pyplot as plt
import numpy as np

snapshots=[item[:44] + "out/snap_" + item[23:28] + '_' + item[34:36] +'_' + item[-14:-4].split('_')[2] for item in gastraces]
infall, ejecta, models = [], [], []

for i, model in enumerate(gastraces):
    s, h, g = pg.prepare_zoom(snapshots[i], gas_trace=model)
    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > 0]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > 0]]
    print np.concatenate(metals_infall).sum(), np.concatenate(metals_ejection).sum(), model[53:61]


# fix, ax = plt.subplots(10)
# ax.scatter(infall, ejectra)

# for i, txt in enumerate(models):
#     ax.annotate(txt, (infall[i], ejectra[i]))

