from pytransit import RoadRunnerModel
import numpy as np
from scipy import optimize


def do_fitting(start_transit, end_transit, model, time, flux, period, depth, fil):
    print('depth = ', depth)
    _mid = (end_transit + start_transit) / 2.
    _k = 0.08
    # choose the limb darkening coefficients depending on the filter
    if fil == 'I':
        _ldc1 = 0.3
        _ldc2 = -0.1
    elif fil == 'R':
        _ldc1 = 0.5
        _ldc2 = 0.01
    elif fil == 'V':
        _ldc1 = 0.7
        _ldc2 = 0.0
    elif fil == 'B':
        _ldc1 = 1
        _ldc2 = -0.1
    else:
        _ldc1 = 0.1
        _ldc2 = -0.05

    if depth > 50:
        _a = 15
    elif depth > 40:
        _a = 30
    elif depth > 30:
        _a = 25
    elif depth > 25:
        _a = 18
    elif depth > 20:
        _a = 16
    elif depth > 15:
        _a = 12
        # _k = 0.1
    elif depth > 10:
        _a = 8
    elif depth > 5:
        _a = 5
    else:
        _a = 4
    if end_transit - start_transit < 0.2:
        _e = 0.2
    else:
        _e = 0.0
    _w = 0.0
    _i = 0.5
    # _a = 15
    tm = RoadRunnerModel(model)  # 'power-2'
    bad_counter = 0
    max_bad = 30
    buf = _k, _a, _w, _i, _e

    # def f(_time, k, a, i, e, w, ldc1, ldc2):
    #     tm.set_data(time=_time)
    #     return tm.evaluate_ps(t0=_mid, k=k, ldc=np.array([ldc1, ldc2]), p=period,
    #                           a=a, i=i * np.pi, e=e, w=w * np.pi)

    def f(_time, k, a, mid, inclination, e, w, ldc1, ldc2):
        tm.set_data(time=_time)
        return tm.evaluate_ps(t0=mid, k=k, ldc=np.array([ldc1, ldc2]), p=period,
                              a=a, i=inclination * np.pi, e=e, w=w * np.pi)
    while True:
        # popt, pcov = optimize.curve_fit(f, time, flux, p0=[_k, _a, _i, _e, _w, _ldc1, _ldc2],
        #                                 bounds=[[0.07, 3, 0.0, 0.0, 0.01, 0, 0],
        #                                         [0.9, 7, 0.6, 0.6, 0.2, 1, 1]])
        popt, pcov = optimize.curve_fit(f, time, flux, p0=[_k, _a, _mid, _i, _e, _w, _ldc1, _ldc2],
                                        bounds=[[0.0, 3, start_transit, 0.0, 0.0, 0.0, 0, -1],
                                                [0.5, 20, end_transit, 0.7, 1, 0.5, 1, 0.1]])
        new_flux = f(time, *popt)
        t_min_index = np.where(abs(new_flux - np.amin(new_flux)) < 2e-3)[0]
        if len(t_min_index) < 2:
            t_min_index = [0, -1]
        else:
            t_min_index = [t_min_index[0], t_min_index[-1]]
        t_max_index = np.where(new_flux == np.amax(new_flux))[0]
        for i in range(1, len(t_max_index)):
            if abs(t_max_index[i] - t_max_index[i - 1]) > 1:
                t_max_index = [t_max_index[i - 1], t_max_index[i]]
                break
        if len(t_max_index) > 2 or len(t_max_index) == 1:
            if time[t_max_index[0]] > _mid:
                t_max_index = [0, t_max_index[0]]
            else:
                t_max_index = [t_max_index[-1], -1]

        if time[t_max_index[0]] > _mid:  # Если левый максимум находится где-то справа от центра,
            # то это фигня какая-то и мы говорим, чтобы максимум был первой точкой слева
            t_max_index[0] = 0

        if time[t_max_index[1]] < _mid:  # Аналогично для правого максимума
            t_max_index[1] = -1

        if t_max_index[0] == t_min_index[0] and t_max_index[-1] == t_min_index[-1] or \
                t_max_index[0] == t_min_index[1]:
            if bad_counter == max_bad:
                print('bad_counter =', max_bad)
                break
            elif bad_counter < 5:
                _a = _a + 0.5
            elif bad_counter < 10:
                _a = buf[1]
                _k = _k + 0.001
            elif bad_counter < 20:
                _k = buf[0]
                _w = _w + 0.01
            elif bad_counter < 25:
                _w = buf[2]
                _e = _e + 0.01
            elif bad_counter < 30:
                _e = buf[4]
                _i = _i + 0.02
            print('bad_counter =', bad_counter)
            bad_counter = bad_counter + 1
        else:
            break
    print('pcov: ')
    print(pcov)
    print('k = ', popt[0])
    print('a = ', popt[1])
    print('mid_time (jd) = ', popt[2])
    print('i Pi = ', popt[3])
    print('e = ', popt[4])
    print('w Pi = ', popt[5])
    print('ldc1 = ', popt[6])
    print('ldc2 = ', popt[7])
    # print('i Pi = ', popt[2])
    # print('e = ', popt[3])
    # print('w Pi = ', popt[4])
    # print('ldc1 = ', popt[5])
    # print('ldc2 = ', popt[6])

    # print(new_flux)
    # t_min_index = np.argmin(new_flux)
    # t_max_index = np.argmax(new_flux)

    # if len(t_max_index) == 1:
    #     if time[t_max_index[0]] > _mid:
    #         t_max_index = [0, t_max_index[0]]
    #     else:
    #         t_max_index = [t_max_index[0], -1]
    print('t_min_index = ', t_min_index)
    print('t_max_index = ', t_max_index)
    # plt.figure(figsize=(7, 7))
    # plt.plot(time, new_flux)
    # plt.plot(time, flux)
    # # plt.axvline(time[t_min_index[0]])
    # # plt.axvline(time[t_min_index[1]])
    # # plt.axvline(time[t_max_index[0]])
    # # plt.axvline(time[t_max_index[1]])
    # plt.grid()
    # plt.show()
    # # plt.savefig(r'D:\docs\sci\tess\good\sended\TIC2305928001\3970\d.pdf')
    # print('t_max_index = ', t_max_index)
    # return new_flux, popt[0], popt[1], _mid, popt[2], popt[3], popt[4], popt[5], popt[6], *t_min_index, *t_max_index
    return new_flux, *popt, *t_min_index, *t_max_index

# f, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3
