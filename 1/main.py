import utils
import pygad as pg

def temperature(args):
    mask = s.gas['ejection_time'][:,0] > '8 Gyr'
    ej_hot = s.gas['mass_at_ejection'][mask][ s.gas['T_at_ejection'][:,0][mask]>1e5 , 0 ].sum()
    ej_cold = s.gas['mass_at_ejection'][mask][ s.gas['T_at_ejection'][:,0][mask]<=1e5 , 0 ].sum()
    print ej_hot / (ej_cold+ej_hot)
    ej_hot = s.gas['mass_at_ejection'][ s.gas['T_at_ejection']>1e5 ].sum()
    ej_cold = s.gas['mass_at_ejection'][ s.gas['T_at_ejection']<=1e5 ].sum()
    print ej_hot / (ej_cold+ej_hot)

    T_ejection = s.gas['T_at_ejection'][s.gas['num_recycled'] > -1]
    T_infall = s.gas['T_at_infall'][s.gas['num_recycled'] > -1]

    for i, temp in enumerate(T_ejection):
        if len(temp) < len(T_infall[i]):
            T_infall[i] = T_infall[i][:-1]
        elif len(temp) > len(T_infall[i]):
            print 'hmm'

    fig, ax = plt.subplots(3)
    plt.tight_layout()
    ax[0].set_xlabel("$T_{infall}$")
    ax[0].set_ylabel("$T_{ejection}$")
    ax[0].set_xlim((1e2, 1e9))
    ax[0].set_ylim((1e2, 1e9))
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].scatter(T_infall, T_ejection, alpha=.1, edgecolor='none')
    ax[0].plot([0, 1e80], [0, 1e80], color='r')

    ax[1].set_ylabel("count")
    ax[1].set_xlabel("$T$ [K]")
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].hist(T_infall.flatten(), bins=np.logspace(2, 9, 100/2), alpha=.5, label='infall')
    ax[1].hist(T_ejection.flatten(), bins=np.logspace(2, 9, 100/2), alpha=.5, label='ejection')
    lgd1 = ax[1].legend(loc='best')

    ax[2].set_title('Initial')
    ax[2].set_ylabel("count")
    ax[2].set_xlabel("$T$ [K]")
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].hist(T_infall[:,0], bins=np.logspace(2, 9, 100/2), alpha=.5, label='initial infall')
    ax[2].hist(T_ejection[:,0], bins=np.logspace(2, 9, 100/2), alpha=.5, label='initial ejection')
    lgd2 = ax[2].legend(loc='best')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')


def profiles(s, h, g):
    r = s.gas['r'][s.gas['num_recycled'] > -1] 
    nrc = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    z = s.gas['metallicity'][s.gas['num_recycled'] > -1]

    nrc_avg, edges_nrc, count_nrc = stats.binned_statistic(r, nrc, statistic='mean', bins=np.linspace(0, 3000, 3000))
    z_avg, edges_z, count_z = stats.binned_statistic(r, z, statistic='mean', bins=np.linspace(0, 3000, 3000))
    
    nrc_avg = utils.prepare_step(nrc_avg)
    z_avg = utils.prepare_step(z_avg)

    fig, (ax1, ax2) = plt.subplots(2)
    plt.tight_layout()
    ax1.grid(True)
    ax1.set_xlim((0, 100))
    ax1.set_ylim((0, np.max(nrc[r < ax1.get_xlim()[1]])))
    ax1.set_xlabel('r [kpc]')
    ax1.set_ylabel('number of recycles')
    ax1.scatter(r, nrc)
    ax1.step(edges_nrc, nrc_avg, c='r')

    ax2.grid(True)
    ax2.set_xlim((0, 100))
    ax2.set_ylim((5e-4, 1e-1))
    ax2.set_yscale(('log'))
    ax2.set_xlabel('r [kpc]')
    ax2.set_ylabel('z')
    ax2.scatter(r, z)
    ax2.step(edges_z, z_avg, c='r')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')


def num_recycled(s, h, g):
    halo = args[0]
    definition = args[1]

    path = '/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_???' % (halo, halo)
    max = int(sorted(glob.glob(path))[-1][-3:])
    s, h, g = pg.prepare_zoom('/ptmp/mpa/naab/REFINED/%s/SF_X/4x-2phase/out/snap_%s_4x_%s' % (halo, halo, max), gas_trace='/u/mihac/data/%s/4x-2phase/gastrace_%s' % (halo, definition), star_form=None)

    x = range(s.gas['num_recycled'].max() + 1)
    y = [s.gas['mass'][s.gas['num_recycled'] == i].sum() for i in range(s.gas['num_recycled'].max() + 1)]
    y = np.array(y)

    fig, ax = plt.subplots(1)
    plt.rcParams.update({'font.size': 16})
    plt.tight_layout()

    ax.set_xlabel("$N_{recycled}$")
    ax.set_ylabel("$M_{\odot}$")
    ax.bar(x, y, width=1)

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')


def nrc_vs_z(s, h, g):
    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    infall_time = [item[item > 0] for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
    ejection_time = [item[item > 0] for item in s.gas['ejection_time'][s.gas['num_recycled'] > -1]]
    number_of_recycles = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    most_recent_z = []

    for i, item in enumerate(metals_infall):
        metals = np.vstack([np.column_stack([item, infall_time[i]]), np.column_stack([metals_ejection[i], ejection_time[i]])])
        most_recent_z.append(metals[np.argmax(metals[:, 1])][0])

    average_nor, edges_nor, count_nor = stats.binned_statistic(most_recent_z, number_of_recycles, statistic='mean', bins=np.linspace(0, .1, 100))

    fix, ax = plt.subplots(1)
    ax.grid(True)
    ax.set_xlabel("z")
    ax.set_ylabel("number of recycles")
    ax.set_xlim((0, .08))
    ax.scatter(most_recent_z, number_of_recycles, label='number of recycles')
    ax.plot(edges_nor[: -1], average_nor, c='r', label='average')
    plt.legend(loc='best')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')


def nrc_vs_r(s, h, g):
    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    cycle_r_max = [item[item > 0] for item in s.gas['cycle_r_max'][s.gas['num_recycled'] > -1]]
    number_of_recycles = s.gas['num_recycled'][s.gas['num_recycled'] > -1]
    extended_nor = []

    for i, particle in enumerate(cycle_r_max):
        for line in particle:
            extended_nor.append(number_of_recycles[i])

    average_nor, edges_nor, count_nor = stats.binned_statistic(np.concatenate(cycle_r_max), extended_nor, statistic='mean', bins=np.linspace(0, 200, 201))
    average_z, edges, count = stats.binned_statistic(np.concatenate(cycle_r_max), np.concatenate(metals_ejection), statistic='mean', bins=np.linspace(0, 200, 201))

    average_z = utils.prepare_step(average_z)
    average_nor = utils.prepare_step(average_nor)

    fig, ax = plt.subplots(1)
    ax.grid(True)
    ax.set_xlabel(('cycle r max [kpc]'))
    ax.set_ylabel('z', color='b')
    ax.step(edges, average_z)
    ax.tick_params('y', colors='b')

    ax2 = ax.twinx()
    ax2.step(edges_nor, average_nor, color='r')
    ax2.set_ylabel('average number of recycles', color='r')
    ax2.tick_params('y', colors='r')
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')


def metals(s, h, g):
    metals_ejection = [item[item > 0] for item in s.gas['metals_at_ejection'][s.gas['num_recycled'] > -1] / s.gas['mass_at_ejection'][s.gas['num_recycled'] > -1]]
    metals_infall = [item[item > 0] for item in s.gas['metals_at_infall'][s.gas['num_recycled'] > -1] / s.gas['mass_at_infall'][s.gas['num_recycled'] > -1]]
    metals_infall_full = copy(metals_infall)

    for i, temp in enumerate(metals_ejection):
        if len(temp) < len(metals_infall[i]):
            metals_infall[i] = metals_infall[i][:-1]
        elif len(temp) > len(metals_infall[i]):
            print 'hmm'

    metals_ejection = [item for sublist in metals_ejection for item in sublist]
    metals_infall = [item for sublist in metals_infall for item in sublist]
    metals_infall_full = [item for sublist in metals_infall_full for item in sublist]

    fig, ax = plt.subplots(3)
    plt.tight_layout()
    ax[0].set_xlabel("z at infall")
    ax[0].set_ylabel("z at ejection")
    ax[0].set_xlim((0, .1))
    ax[0].set_ylim((0, .08))
    ax[0].scatter(metals_infall, metals_ejection)
    ax[0].plot([0, 1], [0, 1], color='r')

    ax[1].set_ylabel("count")
    ax[1].set_xlabel("z")
    ax[1].set_xlim((0, .1))
    ax[1].hist(metals_infall_full, bins=np.linspace(0, .1, 101), alpha=.5, label='infall')
    ax[1].hist(metals_ejection, bins=np.linspace(0, .1, 101), alpha=.5, label='ejection')
    leg1 = plt.legend(loc='best')

    ax[2].set_ylabel("count")
    ax[2].set_xlabel("z")
    ax[2].set_xscale('log')
    ax[2].set_yscale('log')
    ax[2].hist(metals_infall_full, bins=np.logspace(-4, -1, 31), alpha=.5, label='infall')
    ax[2].hist(metals_ejection, bins=np.logspace(-4, -1, 31), alpha=.5, label='ejection')
    leg2 = plt.legend(loc='best')

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')


def in_out_over_time(s, h, g):
    timerange = np.linspace(0, s.cosmic_time(), 13 * 8)
    rec = s.gas[s.gas['num_recycled'] > -1]
    metals_infall = rec['metals_at_infall']
    metals_infall_initial = metals_infall[:, 0]
    metals_infall_reac = metals_infall[:, 1:]
    metals_ejection = rec['metals_at_ejection'] 
    metals_ejection_initial = metals_ejection[:, 0]
    metals_ejection_reac = metals_ejection[:, 1:]

    mass_infall = rec['mass_at_infall']
    mass_infall_initial = mass_infall[:, 0]
    mass_infall_reac = mass_infall[:, 1:]
    mass_ejection = rec['mass_at_ejection'] 
    mass_ejection_initial = mass_ejection[:, 0]
    mass_ejection_reac = mass_ejection[:, 1:]
    
    infall_time = rec['infall_time']
    infall_time_initial = infall_time[:, 0]
    infall_time_reac = infall_time[:, 1:]
    ejection_time = rec['ejection_time']
    ejection_time_initial = ejection_time[:, 0]
    ejection_time_reac = ejection_time[:, 1:]

    metals_infall, edges, count = stats.binned_statistic(np.concatenate(infall_time), np.concatenate(metals_infall), statistic='sum', bins=timerange)
    metals_ejection, edges, count = stats.binned_statistic(np.concatenate(ejection_time), np.concatenate(metals_ejection), statistic='sum', bins=timerange)
    metals_infall_initial, edges_init, count = stats.binned_statistic(infall_time_initial, metals_infall_initial, statistic='sum', bins=timerange)
    metals_ejection_initial, edges_init, count = stats.binned_statistic(ejection_time_initial, metals_ejection_initial, statistic='sum', bins=timerange)
    metals_infall_reac, edges_reac, count = stats.binned_statistic(np.concatenate(infall_time_reac), np.concatenate(metals_infall_reac), statistic='sum', bins=timerange)
    metals_ejection_reac, edges_reac, count = stats.binned_statistic(np.concatenate(ejection_time_reac), np.concatenate(metals_ejection_reac), statistic='sum', bins=timerange)
    mass_infall, edges, count = stats.binned_statistic(np.concatenate(infall_time), np.concatenate(mass_infall), statistic='sum', bins=timerange)
    mass_ejection, edges, count = stats.binned_statistic(np.concatenate(ejection_time), np.concatenate(mass_ejection), statistic='sum', bins=timerange)
    mass_infall_initial, edges_init, count = stats.binned_statistic(infall_time_initial, mass_infall_initial, statistic='sum', bins=timerange)
    mass_ejection_initial, edges_init, count = stats.binned_statistic(ejection_time_initial, mass_ejection_initial, statistic='sum', bins=timerange)
    mass_infall_reac, edges_reac, count = stats.binned_statistic(np.concatenate(infall_time_reac), np.concatenate(mass_infall_reac), statistic='sum', bins=timerange)
    mass_ejection_reac, edges_reac, count = stats.binned_statistic(np.concatenate(ejection_time_reac), np.concatenate(mass_ejection_reac), statistic='sum', bins=timerange)

    dt = timerange[1]*1e9

    metals_infall /=  dt
    metals_ejection /= dt
    metals_infall_initial /= dt
    metals_ejection_initial /= dt
    metals_infall_reac /= dt
    metals_ejection_reac /= dt
    mass_infall /= dt
    mass_ejection /= dt
    mass_infall_initial /= dt
    mass_ejection_initial /= dt
    mass_infall_reac /= dt
    mass_ejection_reac /= dt


    metals_infall = utils.prepare_step(metals_infall)
    metals_ejection = utils.prepare_step(metals_ejection)
    metals_infall_initial = utils.prepare_step(metals_infall_initial)
    metals_ejection_initial = utils.prepare_step(metals_ejection_initial)
    metals_infall_reac = utils.prepare_step(metals_infall_reac)
    metals_ejection_reac = utils.prepare_step(metals_ejection_reac)
    mass_infall = utils.prepare_step(mass_infall)
    mass_ejection = utils.prepare_step(mass_ejection)
    mass_infall_initial = utils.prepare_step(mass_infall_initial)
    mass_ejection_initial = utils.prepare_step(mass_ejection_initial)
    mass_infall_reac = utils.prepare_step(mass_infall_reac)
    mass_ejection_reac = utils.prepare_step(mass_ejection_reac)


    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(25, 20))
    plt.tight_layout()
    for ax in (ax1, ax2, ax3, ax4):
        ax.grid(True)
        ax.set_xlabel("Time [yr]")
        ax.set_xlim((0, s.cosmic_time()))

    ax1.set_ylabel("mass rate [$M_{\odot}/yr$]")
    ax2.set_ylabel("mass rate [$M_{\odot}/yr$]")
    ax3.set_ylabel("metal rate [$M_{\odot}/yr$]")
    ax4.set_ylabel("metal rate [$M_{\odot}/yr$]")

    ax1.set_title("Mass")
    ax1.step(edges, mass_infall, label='Total infall')
    ax1.step(edges_init, mass_infall_initial, label='First infall')
    ax1.step(edges_reac, mass_infall_reac, label='Subsequent infall')
    lgd1 = ax1.legend(loc='upper left')
    ax2.step(edges, mass_ejection, label='Total ejection')
    ax2.step(edges_init, mass_ejection_initial, label='First ejection')
    ax2.step(edges_reac, mass_ejection_reac, label='Subsequent ejection')
    lgd2 = ax2.legend(loc='upper left')

    ax3.set_title("Metals")
    ax3.step(edges, metals_infall, label='Total infall')
    ax3.step(edges_init, metals_infall_initial, label='First infall')
    ax3.step(edges_reac, metals_infall_reac, label='Subsequent ejection')
    lgd3 = ax3.legend(loc='upper left')
    ax4.step(edges, metals_ejection, label='Total ejection')
    ax4.step(edges_init, metals_ejection_initial, label='First ejection')
    ax4.step(edges_reac, metals_ejection_reac, label='Subsequent ejection')
    lgd4 = ax4.legend(loc='upper left')

    utils.mkdir('./%s' % (halo))
    plt.savefig('./%s/' % (halo) + ("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')

def time_in(s, array):
    array = np.array(array[array > 0])
    return pg.UnitArr(np.r_[array[1:], float(s.cosmic_time())] - array, 'Gyr')[:-1]

def infall_time(s, h, g):
    data = [time_in(s, item) for item in s.gas['infall_time'][s.gas['num_recycled'] > -1]]
    data = [item for sublist in data for item in sublist]

    fig, ax = plt.subplots(1)
    ax.set_xlabel("infall time [Gyr]")
    ax.set_ylabel("count")
    ax.hist(data, bins=np.linspace(0, 1.5, 50))

    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

def density_profiles(s, h, g):
    f, ax = plt.subplots(3, 2, figsize=(20,12))
    
    for axes in ax.flatten():
        axes.grid(True)

    pg.plotting.profile(s.gas, '200 kpc', 'mass', ax=ax[0,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'mass', ax=ax[0,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'mass', ax=ax[0,0], label='recycled')
    ax[0,0].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'MgII', ax=ax[1,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'MgII', ax=ax[1,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'MgII', ax=ax[1,0], label='recycled')
    ax[1,0].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'CIV', ax=ax[0,1], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'CIV', ax=ax[0,1], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'CIV', ax=ax[0,1], label='recycled')
    ax[0,1].set_ylim((1e-4, 1e0))
    ax[0,1].legend(loc='best')

    pg.plotting.profile(s.gas, '200 kpc', 'HI', ax=ax[1,1], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'HI', ax=ax[1,1], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'HI', ax=ax[1,1], label='recycled')
    ax[1,1].legend(loc='best')
    
    pg.plotting.profile(s.gas, '200 kpc', 'OVI', ax=ax[2,0], label='total')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] == -1], '200 kpc', 'OVI', ax=ax[2,0], label='no_recycling')
    pg.plotting.profile(s.gas[s.gas['num_recycled'] > -1], '200 kpc', 'OVI', ax=ax[2,0], label='recycled')
    ax[2,0].set_ylim((1e-3, 1e-1))
    ax[2,0].legend(loc='best')
    
    f.tight_layout()
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + type + ".png", bbox_inches='tight')

def cycle_r_max(s, h, g):
    cycle_r_max = s.gas['cycle_r_max'] #[s.gas['num_recycled'] > -1]
    particle_mass = s.gas['mass_at_ejection'] #[s.gas['num_recycled'] > -1]
    mask = np.max(s.gas['cycle_r_max'], axis=1) > '100 kpc'
    cycle_r_max_lim = cycle_r_max[mask]
    ID_mask = pg.IDMask(s.gas['ID'][mask])


    mass_in_bins, edges2, count2 = stats.binned_statistic(cycle_r_max.flatten(),
                                     particle_mass.flatten(), statistic='sum', bins=np.logspace(-1.1, 4, 90))
    mass_in_bins /= 1e10
    
    fig, ax = plt.subplots(4, figsize=(16,12))
    ax[0].set_xlabel("$r_{max}$ [kpc]")
    ax[0].set_ylabel("mass [$10^{10}\ M_{\odot}$]")
    ax[0].set_xscale('log')
    ax[0].set_xlim((1e-2, 1e4))
    ax[0].set_yscale('log')
    ax[0].grid(True)
    ax[0].bar(edges2[:-1], mass_in_bins, width=np.diff(edges2))

    ax[1].set_xlabel("$r_{max}$ [kpc]")
    ax[1].set_ylabel("count")
    ax[1].grid(True)
    ax[1].set_yscale('log')
    ax[1].hist(cycle_r_max.flatten(), bins=np.linspace(0, 300, 301))
    
    ax[2].set_title("same as above but logscale")
    ax[2].set_xlabel("$r_{max}$ [kpc]")
    ax[2].set_ylabel("count")
    ax[2].grid(True)
    ax[2].set_yscale('log')
    ax[2].set_xscale('log')
    ax[2].hist(cycle_r_max.flatten(), bins=np.logspace(-1.1, 4, 90), log=True)

    ax[3].set_title("same as above but just particles with max(cycle_r_max) > 100)")
    ax[3].set_xlabel("$r_{max}$ [kpc]")
    ax[3].set_ylabel("count")
    ax[3].grid(True)
    ax[3].set_yscale('log')
    ax[3].set_xscale('log')
    ax[3].hist(cycle_r_max_lim.flatten(), bins=np.logspace(-1.1, 4, 90), log=True)

    fig.tight_layout()
    
    plt.savefig(filename.split("/")[-1][:-3] + '_' + halo + '_' + definition + ".png", bbox_inches='tight')
    plt.clf()
    plt.close()
    #for snap in sorted(glob.glob(path))[70::20]:
    #    s, h, g = pg.prepare_zoom(snap, gas_trace=None, star_form=None)
    #    fig, ax = plt.subplots(1, 3)
    #    fig.suptitle(s.cosmic_time())
    #    pltargs = dict(clim=(10**3.5, 10**6.5))
    #    pg.plotting.image(s.gas, extent='400 kpc', ax=ax[0], **pltargs)
    #    pg.plotting.image(s.gas[ID_mask], extent='400 kpc', ax=ax[1], **pltargs)
    #    pg.plotting.image(s.gas[~ID_mask], extent='400 kpc', ax=ax[2], **pltargs)
    #    plt.tight_layout()
    #    print '%s/%s/%s_%s.png' % (halo, definition, definition, snap.split('/')[-1])
    #    plt.savefig('%s/%s/%s_%s.png' % (halo, definition, definition, snap.split('/')[-1]))
    #    plt.clf()
    #    plt.close()

