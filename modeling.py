import numpy as np
from astropy.io import ascii
from astropy.table import Table
from pytransit import RoadRunnerModel


def find_chi2(model, path, k, a, mid_time_bjd, inclination, e, w, ldc, perr):
    period = 2.103193

    def fuck(_time, k, a, mid, inclination, e, w, ldc1, ldc2, tm):
        tm.set_data(time=_time)
        return tm.evaluate_ps(t0=mid, k=k, ldc=np.array([ldc1, ldc2]), p=period,
                              a=a, i=inclination * np.pi, e=e, w=w * np.pi)

    file_name = path.split('/')[-1]
    file_name = file_name.split('.')[:-1]
    info_file = ascii.read(path, delimiter='\t')

    rel_flux_t_1 = info_file['rel_flux_T1']
    rel_flux_err_t_1 = info_file['rel_flux_err_T1']
    time_bjd_tdb = info_file['BJD_TDB']
    const = (1 - np.mean(rel_flux_t_1[-30:]))
    print('const =', const)
    rel_flux_t_1 = rel_flux_t_1 + const

    tm = RoadRunnerModel(model)
    model_flux = fuck(time_bjd_tdb, k, a, mid_time_bjd, inclination, e, w, ldc[0], ldc[1], tm)

    Y = sum(((rel_flux_t_1 - model_flux) / rel_flux_err_t_1) ** 2)
    print('Y =', Y)
    chi2 = Y / (len(rel_flux_t_1) - 8)
    print('chi2 =', chi2)

    info = Table(
        [[model], [k], [perr[0]], [a], [perr[1]], [mid_time_bjd], [perr[2]], [inclination], [perr[3]], [e], [perr[4]],
         [w], [perr[5]], [ldc[0]], [perr[6]], [ldc[1]], [perr[7]], [Y], [chi2]],
        names=(
            'model', 'k', 'err_k', 'a', 'err_a', 'mid_time_bjd', 'err_mid', 'inclination', 'err_i', 'e',
            'err_e',
            'w', 'err_w', 'ldc1', 'err_ldc1', 'ldc2', 'err_ldc2', 'Y', 'chi2'),
        meta={'name': 'info'})

    ascii.write(info, file_name[0] + file_name[1] + file_name[1] + '_' + model + '.txt', delimiter='\t',
                fast_writer=False)


path = r'D:\docs\курсач3\for_cursuch\bieryla\TIC229510886.01_UT2019.1014_KeplerCam_z_measurements.xls'
find_chi2('power-2', path,
          k=0.0866614780821258,
          a=3.6553997658544817,
          mid_time_bjd=2458764.354557003,
          inclination=0.534385904288784,
          e=0.13242893868072736,
          w=0.6125616852176173,
          ldc=[0.9427851107990236, 0.8284505596556796],
          perr=[6.00991384e-02, 1.50828492e+02, 6.06192167e-03, 2.69830954e+00,
                4.15834516e+01, 4.29486113e+01, 5.38389765e+00, 7.53288138e+00])
find_chi2('quadratic', path,
          k=0.08744962402402959,
          a=3.6529308647345315,
          mid_time_bjd=2458764.354580735,
          inclination=0.5312827742623156,
          e=0.1581028296762522,
          w=0.4249430238768369,
          ldc=[0.7116513556964006, 0.16120624304586534],
          perr=[5.17352756e-02, 1.66988384e+02, 6.08041369e-03, 2.56903921e+00,
                4.79936021e+01, 2.86355103e+01, 3.48495730e+00, 3.73189292e+00])
find_chi2('uniform', path,
          k=0.0803647073272396,
          a=3.5998110956479605,
          mid_time_bjd=2458764.3352059447,
          inclination=0.5000000002186413,
          e=0.09989261792686133,
          w=0.4800067801660832,
          ldc=[0.7, 0.05],
          perr=[2.55650655e-03, 2.94084775e+01, 4.69117776e-03, 4.77251602e+04,
                8.37768664e+00, 1.68663131e+01, 0.00000000e+00, 0.00000000e+00])
