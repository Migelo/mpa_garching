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

def calc_ew(los, s):
    line_name = 'Lyman_alpha'
    line = pg.analysis.absorption_spectra.lines[line_name]
    v_limits = pg.UnitArr([-500,500],'km/s')
    xaxis, yaxis = 0, 1
    tau, dens, temp, v_edges, _ = pg.analysis.mock_absorption_spectrum_of(
                s, los, line_name, v_limits,
                xaxis=xaxis, yaxis=yaxis,
                method='particles')

    # v_edges are the edges of the bins in velocity space (rest frame velocity
    # that is) convert to observed wavelengths:
    z_edges = pg.analysis.velocities_to_redshifts(v_edges, z0=s.redshift)
    l_edges = pg.UnitScalar(line['l']) * (1.0 + z_edges)
    
    EW = pg.analysis.absorption_spectra.EW(tau, l_edges).in_units_of('Angstrom')
    return EW

def plot(args):
    halo = args[0]
    definition = args[1]
    modification = ''
    if '-' in halo:
        modification = halo[5:]
        halo = halo[:5]
    if definition.split('_')[-1] != definition:
        definition_filename = definition.split('_')[-1]
        definition = definition.split('_')[0]
    else:
        definition_filename = ''
    print args
    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_???' % (halo, modification, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    print halo, definition
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase%s/out/snap_%s_4x_%s' % (halo, modification, halo, max),
        gas_trace=None, star_form=None)

    line_name = 'Lyman_alpha'
    line = pg.analysis.absorption_spectra.lines[line_name]
    
    # set the rotation angle
    rotation_angles = np.random.rand(3) * 2 * np.pi
    diag = np.diag(np.ones(3))
    v_limits = pg.UnitArr([-500, 500],'km/s')
    
    n_rot = 21
    n_lines = 33
    # ews = np.zeros((n_rot*n_lines, 2))
    ewss = []
    rs = []
    for j in range(n_rot):
        ews = []
        r = []
        for i, angle in enumerate(np.random.rand(3) * 2 * np.pi):
            pg.transformation.rot_from_axis_angle(diag[i], angle).apply(s, total=True)
        for k in range(n_lines):
            impact_parameter = (np.random.rand() * 200 + 50) * np.random.choice([-1, 1])
            x = np.random.rand() * impact_parameter * np.random.choice([-1, 1])
            y = np.sqrt(impact_parameter**2 - x**2) * np.random.choice([-1, 1])
    #         ews[j*n_rot + k, 0] = impact_parameter
    #         ews[j*n_rot + k, 1] = calc_ew(pg.UnitArr([x, y], 'kpc'), s)
            ews.append(calc_ew(pg.UnitArr([x, y], 'kpc'), s))
            r.append(impact_parameter)
        ewss.append(ews)
        rs.append(r)
        
    f, ax = plt.subplots(7, 3, figsize=(5*8, 8*4))
    f.suptitle(halo, fontsize=44)
    for j in range(n_rot):
        ax[j//3, j%3].scatter(np.abs(rs[j]), ewss[j])
        ax[j//3, j%3].set_yscale('log')
        ax[j//3, j%3].set_ylim(10**-1, 6)
        ax[j//3, j%3].set_xlim(0, 260)
        ax[j//3, j%3].set_ylabel('EW [angstrom]', fontsize=15)
        ax[j//3, j%3].grid(True)
        ax[j//3, j%3].set_xlabel(r'impact parameter $\rho$ [kpc]', fontsize=15)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '.png', bbox_inches='tight')

p = Pool(4)
combinations = utils.combinations
#combinations = [('M0858', 'ism')]
p.map(plot, combinations)
