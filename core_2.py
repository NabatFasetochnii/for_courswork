import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.time import Time as aTime
from pytransit import RoadRunnerModel

from Transit_Fitting import do_fitting


def func(model):
    def fuck(_time, k, a, mid, inclination, e, w, ldc1, ldc2, tm):
        tm.set_data(time=_time)
        return tm.evaluate_ps(t0=mid, k=k, ldc=np.array([ldc1, ldc2]), p=period,
                              a=a, i=inclination * np.pi, e=e, w=w * np.pi)

    def draw_my_models(local_model, color):
        tm = RoadRunnerModel(local_model)
        axs.plot(time_bjd_tdb - ZERO, fuck(time_bjd_tdb, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, tm), 'y.',
                 markersize=3, zorder=4,
                 linewidth=0.5, label='model, Limb darkening: ' + local_model + ', k = ' + '{:.3f}'.format(k_model) +
                                      ', a = ' + '{:.2f}'.format(a) + ', i = ' +
                                      '{:.3f}'.format(i_model) + 'Pi, '
                                      + 'e = ' + '{:.2f}'.format(e) + ', w = ' + '{:.3f}'.format(w) + 'Pi',
                 linestyle='solid', color=color)  # рисуем модель

    Scope = 'Коуровская астрономическая обсерватория - 0.4'
    period = 2.103193
    fil = 'V'
    # model = 'power-2'
    # model = 'quadratic'
    # model = 'uniform'
    # model = 'all'
    path = r'D:\docs\sci\tess\good\for_cursuch\me'
    name = 'TIC229510866.01_20201211_Kourovka_0.4_V_lightcurve'
    file_path = path + '/' + name + '.dat'
    info_file = ascii.read(file_path, delimiter='\t')

    rel_flux_t_1 = info_file['Rel_Flux_Target']
    rel_flux_err_t_1 = info_file['Rel_Flux_Err_Target']

    Date = aTime(info_file['JD'], format='jd')
    time_bjd_tdb = info_file['BJD_TDB']
    Depth = 1 / 10 ** (7.246375 / 2500)

    # create Title
    Title = 'TIC229510866.01'
    Title += '\n'
    Title += Scope + ', ' + Date[0].datetime.strftime('%Y%m%d')
    Title += ', filter =' + fil
    Exp = 50
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

    # draw light curves
    fig, axs = plt.subplots(1, 1, figsize=(8, 4), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
    fig.suptitle(Title, fontsize=8)

    if model != 'all':
        # modeling
        f, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3, perr = \
            do_fitting(time_bjd_tdb[0], time_bjd_tdb[-1], model, time_bjd_tdb,
                       rel_flux_t_1, period, (1 - Depth) * 1000, fil)

        axs.plot(time_bjd_tdb - ZERO, f, 'y.', markersize=3, zorder=4,
                 linewidth=0.5, label='model, Limb darkening: ' + model + ', k = ' + '{:.3f}'.format(k_model) +
                                      ', a = ' + '{:.2f}'.format(a) + ', i = ' +
                                      '{:.3f}'.format(i_model) + 'Pi, '
                                      + 'e = ' + '{:.2f}'.format(e) + ', w = ' + '{:.3f}'.format(w) + 'Pi',
                 linestyle='solid')  # рисуем модель
    else:
        k_model = 0.07389518126082138
        a = 3.4018730472398375
        mid_time = 2459195.472651039
        i_model = 0.5686121161939944
        e = 0.22719820847355754
        w = 0.6969896782139561
        ldc1 = 0.2660520578224939
        ldc2 = 0.015410786651013053

        draw_my_models('power-2', 'k')

        k_model = 0.06456762018805841
        a = 3.7267069254338403
        mid_time = 2459195.4734767354
        i_model = 0.4822899765082747
        e = 0.1766343217140779
        w = 0.5450184875231457
        ldc1 = 0.7319079875392563
        ldc2 = 0.3299215584825351

        draw_my_models('quadratic', 'r')

        k_model = 0.0688420913579282
        a = 3.7368329673526715
        mid_time = 2459195.4776497125
        i_model = 0.5000023451450417
        e = 0.190390459382966
        w = 0.5120376737176133
        ldc1 = 0.7
        ldc2 = 0.05

        draw_my_models('uniform', 'lime')

    axs.errorbar(time_bjd_tdb - ZERO, rel_flux_t_1,
                 rel_flux_err_t_1, fmt='b.', zorder=3,
                 label='TIC 229510866.01' + r', 1$\sigma$ errorbar')
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
    if model != 'all':
        out = 'k = {}\na = {}\nmid_time_bjd = {}\ninclination = {}\n' \
              'e = {}\nw = {}\nldc = [{}, {}]\nt1 = {}, t2 = {}\nt0 = {},' \
              ' t3 = {}\nperr = {}'.format(k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3, perr)

        o = open(path + '/' + 'out_%s.txt' % model, 'w')
        o.write(out)


warnings.simplefilter("ignore")
# func('power-2')
# func('quadratic')
# func('uniform')
func('all')
