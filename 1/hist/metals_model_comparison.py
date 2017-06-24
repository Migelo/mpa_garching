import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('different_models.dat', dtype=[('myint','f8'),('myfloat','f8'),
('mystring','S5')])

infall, ejecta, models = [], [], []

for item in data:
    infall.append(item[0])
    ejecta.append(item[1])
    models.append(item[2])

fix, ax = plt.subplots(1)
ax.scatter(infall, ejecta)
ax.set_xlim(3e8, max(ejecta)*1.1)
ax.set_ylim(3e8, max(ejecta)*1.1)
ax.plot([0, 1e11], [0, 1e11])
ax.grid(True)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('commulative metal infall $M_{\odot}$')
ax.set_ylabel('commulative metal ejection $M_{\odot}$')

for i, txt in enumerate(models):
    ax.annotate(txt, (infall[i], ejecta[i]))

filename = __file__
plt.savefig(filename.split("/")[-1][:-3] + ".pdf", bbox_inches='tight')

