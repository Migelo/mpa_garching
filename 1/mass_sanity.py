import numpy as np
import matplotlib.pyplot as plt
import pygad as pg
from scipy import stats
import glob
from multiprocessing import Pool
import utils

s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M1196/SF_X/4x-2phase/out/snap_M1196_4x_470', gas_trace='/u/mihac/data/M1196/4x-2phase/gastrace_disc', star_form=None)
s1, h1, g1 = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/M1196/SF_X/4x-2phase/out/snap_M1196_4x_070', gas_trace=None, star_form=None)

R200_frac, Tcrit, rhocrit = [.15, '2e4 K', '1e-2 u/cm**3']
R200, M200 = pg.analysis.virial_info(s1)
s1_mass = s1.gas[pg.BallMask(R200_frac*R200) & \
                       pg.ExprMask('(temp < "%s") & (rho > "%s")' % (Tcrit,rhocrit)) ]
s1_mass = s1_mass['mass'].sum()
s_mass = s.gas[pg.BallMask(R200_frac*R200) & \
                       pg.ExprMask('(temp < "%s") & (rho > "%s")' % (Tcrit,rhocrit)) ]
s_mass = s_mass['mass'].sum()

print (g.stars['mass'].sum() - g1.stars['mass'].sum()) / (s.gas['mass_at_infall'].sum() - s.gas['mass_at_ejection'].sum() + s1_mass - s_mass) 

