import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.time import Time as aTime
from pytransit import RoadRunnerModel

from Transit_Fitting import do_fitting


def func(model):
    # def f(_time, k, a, mid, inclination, e, w, ldc1, ldc2, tm):
    #     tm.set_data(time=_time)
    #     return tm.evaluate_ps(t0=mid, k=k, ldc=np.array([ldc1, ldc2]), p=period,
    #                           a=a, i=inclination * np.pi, e=e, w=w * np.pi)
    #
    # def draw_my_models(local_model, color):
    #     tm = RoadRunnerModel(local_model)
    #     axs.plot(time_bjd_tdb - ZERO, f(time_bjd_tdb, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, tm), 'y.',
    #              markersize=3, zorder=4,
    #              linewidth=0.5, label='model, Limb darkening: ' + local_model + ', k = ' + '{:.3f}'.format(k_model) +
    #                                   ', a = ' + '{:.2f}'.format(a) + ', i = ' +
    #                                   '{:.3f}'.format(i_model) + 'Pi, '
    #                                   + 'e = ' + '{:.2f}'.format(e) + ', w = ' + '{:.3f}'.format(w) + 'Pi',
    #              linestyle='solid', color=color)  # рисуем модель

    Scope = 'The George Mason University Observatory - 0.8'
    period = 2.103193
    fil = 'R'
    # model = 'power-2'
    # model = 'quadratic'
    # model = 'uniform'
    # model = 'all'
    path = r'D:\docs\sci\tess\good\for_cursuch\plavchan'
    name = 'TIC229510866-01_20190923_GMU_R_measurements'
    file_path = path + '/' + name + '.xls'
    info_file = ascii.read(file_path, delimiter='\t')

    rel_flux_t_1 = info_file['rel_flux_T1']
    rel_flux_err_t_1 = info_file['rel_flux_err_T1']

    Date = aTime(info_file['JD_UTC'], format='jd')
    time_bjd_tdb = info_file['BJD_TDB']
    Depth = 1 / 10 ** (7.246375 / 2500)

    # create Title
    Title = 'TIC229510866.01'
    Title += '\n'
    Title += Scope + ', ' + Date[0].datetime.strftime('%Y%m%d')
    Title += ', filter =' + fil
    Exp = info_file['EXPTIME'][0]
    Title += ', Extime=' + '{:.1f}'.format(Exp)
    Title += '\n'
    Title += 'Observation start=' + Date[0].datetime.strftime('%Y %m %d %H:%M:%S')
    Title += ', '
    Title += 'end=' + Date[-1].datetime.strftime('%Y %m %d %H:%M:%S')
    Title += ', ' + str(np.round((Date[-1].jd - Date[0].jd) * 24 * 60, 0)) + ' min.'
    Title += ', ' + str(len(Date)) + ' frames'

    const = (1 - np.mean(rel_flux_t_1[-30:]))
    print('const =', const)

    ZERO = int(time_bjd_tdb[0])
    print('zero =', ZERO)

    # modeling
    f, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3, perr = \
        do_fitting(time_bjd_tdb[0], time_bjd_tdb[-1], model, time_bjd_tdb,
                   rel_flux_t_1 + const, period, (1 - Depth) * 1000, fil)

    # draw light curves
    fig, axs = plt.subplots(1, 1, figsize=(8, 4), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
    fig.suptitle(Title, fontsize=8)

    # k_model = 0.05174269848625799
    # a = 3.6382803690514085
    # mid_time = 2458749.60140155
    # i_model = 0.502134121987164
    # e = 0.13693814892993625
    # w = 0.5438555025458212
    # ldc1 = 0.7127979625328335
    # ldc2 = 2.5975168344883337
    #
    # draw_my_models('power-2', 'k')
    #
    # k_model = 0.05968874417686405
    # a = 3.0008410601507713
    # mid_time = 2458749.602221113
    # i_model = 0.5705203987817031
    # e = 1.198594291633736e-07
    # w = 0.4000000000000001
    # ldc1 = 0.9999999999999999
    # ldc2 = 0.29057987674903496
    #
    # draw_my_models('quadratic', 'r')
    #
    # k_model = 0.05445961339098987
    # a = 3.601393618341107
    # mid_time = 2458749.6048065666
    # i_model = 0.5000000015116634
    # e = 0.10241006966808901
    # w = 0.4900235531678831
    # ldc1 = 1
    # ldc2 = 0.5
    #
    # draw_my_models('uniform', 'lime')

    axs.plot(time_bjd_tdb - ZERO, f, 'y.', markersize=3, zorder=4,
             linewidth=0.5, label='model, Limb darkening: ' + model + ', k = ' + '{:.3f}'.format(k_model) +
                                  ', a = ' + '{:.2f}'.format(a) + ', i = ' +
                                  '{:.3f}'.format(i_model) + 'Pi, '
                                  + 'e = ' + '{:.2f}'.format(e) + ', w = ' + '{:.3f}'.format(w) + 'Pi',
             linestyle='solid')  # рисуем модель

    axs.errorbar(time_bjd_tdb - ZERO, rel_flux_t_1 + const,
                 rel_flux_err_t_1, fmt='b.', zorder=3,
                 label='TIC 229510866.01' + r', 1$\sigma$ errorbar')
    # ' depth = ' + '{:.3f}'.format(transit_depth_median[num]) + '±' +
    # '{:.3f}'.format(std_transit_depth[num]))  # рисуем данные
    axs.legend(loc=2, fontsize=6)
    axs.set_ylabel('Normalized flux', fontsize=6)

    locs = axs.get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        # print(str(round(x.tdb.jd - ZERO, 2)))
        x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    axs.set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    axs.set_xlabel('BJD_TDB - ' + str(ZERO), fontsize=6)

    axs.tick_params(axis='both', labelsize=6, direction='in')
    axs.grid()

    plt.savefig(path + '/' + name + '_%s.pdf' % model)
    fig.show()
    out = 'k = {}\na = {}\nmid_time_bjd = {}\ninclination = {}\ne = {}\nw = {}\nldc = [{}, {}]\nt1 = {}, t2 = {}\nt0 = {' \
          '}, t3 = {}\nperr = {}'.format(k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3, perr)

    o = open(path + '/' + 'out_%s.txt' % model, 'w')
    o.write(out)


warnings.simplefilter("ignore")
func('power-2')
func('quadratic')
func('uniform')

# func('all')
