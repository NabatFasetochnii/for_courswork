# #disable warnings
import copy
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.time import Time as aTime
from astropy.timeseries import TimeSeries, aggregate_downsample
from scipy import stats

from Condition_Report import condition_report
from Transit_Fitting import do_fitting

# #from astropy.convolution import convolve, Gaussian1DKernel

warnings.simplefilter("ignore")


def draw_my_annotate(ax, t):
    ax.annotate('{:.3f}'.format(t), xy=(t, 1.011), xytext=(t, 0.98), va='center', ha='center',
                arrowprops={'width': 0.1, 'headwidth': 0, 'linestyle': '-', 'color': 'g'},
                fontsize=4)


def what_transit(time, T0, T1):
    deltaTime = [time[i] - time[i - 1] for i in range(1, len(time + 1))]
    gap = max(deltaTime) > stats.mode(deltaTime)[0]
    ingress = time[0] < T0
    egress = time[-1] > T1
    if ingress and egress and not gap:
        return ' a full transit'
    elif ingress and egress and gap:
        return ' a full with gapped'
    elif ingress and not egress and not gap:
        return ' an ingress '
    elif ingress and not egress and gap:
        return ' a gapped ingress'
    elif not ingress and egress and not gap:
        return ' an egress'
    elif not ingress and egress and gap:
        return ' a gapped egress'


def do_post_proc(Path2Data):
    # Path2Data = path
    Scope = 'Kourovka_0.4'

    #
    ########################################################################
    # InfoFile = 'TTF_Info.txt'
    CooFile = 'Coo.txt'
    TimeFile = 'Time.txt'
    FluxFiles = ["Flux1.txt", "Flux2.txt", "Flux3.txt"]
    SNFiles = ["SN1.txt", "SN2.txt", "SN3.txt"]
    MagFiles = ["Mag1.txt", "Mag2.txt", "Mag3.txt"]
    MerrFiles = ["Merr1.txt", "Merr2.txt", "Merr3.txt"]
    SkyFile = "Sky.txt"
    MaxFile = 'Max.txt'

    # FluxFiles.append("Flux4.txt")
    # SNFiles.append("SN4.txt")
    # MagFiles.append("Mag4.txt")
    # MerrFiles.append("Merr4.txt")

    ########################################################################
    #  #read data
    # Info = ascii.read(Path2Data + '/' + InfoFile, delimiter='\t', format='commented_header')
    Cat = [ascii.read(Path2Data + '/' + CooFile, delimiter='\t', format='commented_header')
           for item in range(len(FluxFiles))]
    Time = [ascii.read(Path2Data + '/' + TimeFile, delimiter='\t') for item in range(len(FluxFiles))]
    Flux = [np.genfromtxt(Path2Data + '/' + FluxFile) for FluxFile in FluxFiles]
    SN = [np.genfromtxt(Path2Data + '/' + SNFile) for SNFile in SNFiles]
    Mag = [np.genfromtxt(Path2Data + '/' + MagFile) for MagFile in MagFiles]
    Error = [np.genfromtxt(Path2Data + '/' + MerrFile) for MerrFile in MerrFiles]
    Sky = [np.genfromtxt(Path2Data + '/' + SkyFile) for item in range(len(FluxFiles))]
    Max = [np.genfromtxt(Path2Data + '/' + MaxFile) for item in range(len(FluxFiles))]

    NEB_Flux = copy.deepcopy(Flux)

    #  #transet start-end
    # T0 = aTime(Info['jd_start'][0] + 2450000, format='jd')
    # T1 = aTime(Info['jd_end'][0] + 2450000, format='jd')

    ########################################################################
    #  #set prefics
    Date = aTime(Time[0]['JD'][0], format='jd')
    Pref = "TIC229510866.01" '_' + \
           Date.datetime.strftime('%Y%m%d') + '_' + Scope + '_' + Time[0]['FILTER'][0]

    #  #condition report
    condition_report(Time[0], Max[0], Path2Data + '/' + Pref + '_conditions.pdf', Scope)

    ########################################################################
    #  #plot raw data
    # plt.figure(figsize=(7, 7))
    # plt.plot(Time['JD'], Mag)
    # plt.ylabel('Mag', fontsize=6)
    # plt.title('Raw data', fontsize=8)
    # plt.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
    # locs, labels = plt.xticks()
    # t = aTime(locs, format='jd')
    # x_ticks_labels = []
    # for x in t:
    #     x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
    # plt.xticks(locs, x_ticks_labels, rotation='vertical', fontsize=6)
    # plt.xlabel('Date-Time (UTC), ' + Time['DATE-OBS'][0].split('T')[0], fontsize=6)
    # #  #plt.xlim(locs[1], locs[-2])
    # plt.tick_params(axis='both', labelsize=6, direction='in')
    # plt.grid()
    # plt.show()

    ########################################################################
    #  #delete bad stars
    #  #print(np.isnan(Flux)[0])
    Indexes = [np.where(np.isnan(item)[0]) for item in Flux]
    #  #print(Index)
    #  #print(Cat)
    for item in range(len(Indexes)):
        Flux[item] = np.delete(Flux[item], Indexes[item], axis=1)
        SN[item] = np.delete(SN[item], Indexes[item], axis=1)
        Mag[item] = np.delete(Mag[item], Indexes[item], axis=1)
        Error[item] = np.delete(Error[item], Indexes[item], axis=1)
        Sky[item] = np.delete(Sky[item], Indexes[item], axis=1)
        Max[item] = np.delete(Max[item], Indexes[item], axis=1)
        Cat[item].remove_rows(Indexes[item])

    ########################################################################
    #  Search bad frames
    bad_threshold = 0.25
    # print('Search bad frames, threshold =', bad_threshold)
    diffs = [mag[:, 0] - np.nanmedian(mag[:, 0]) for mag in Mag]
    # plt.figure(figsize=(7, 7))
    # plt.plot(Time['JD'], diff)
    # plt.hlines(Bad_threshold, locs[1], locs[-2])
    # plt.ylabel('Extinction, mag', fontsize=6)
    # plt.title('Bad frames', fontsize=8)
    # plt.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
    # plt.xticks(locs, x_ticks_labels, rotation='vertical', fontsize=6)
    # plt.xlabel('Date-Time (UTC), ' + Time['DATE-OBS'][0].split('T')[0], fontsize=6)
    # #  plt.xlim(locs[1], locs[-2])
    # plt.tick_params(axis='both', labelsize=6, direction='in')
    # plt.grid()
    # plt.show()

    #  delete bad frames
    for item in range(len(diffs)):
        Index = np.where(diffs[item] < bad_threshold)[0]
        Time[item] = Time[item][Index]
        Flux[item] = Flux[item][Index]
        SN[item] = SN[item][Index]
        Mag[item] = Mag[item][Index]
        Error[item] = Error[item][Index]
        Sky[item] = Sky[item][Index]
        Max[item] = Max[item][Index]
        NEB_Flux[item] = NEB_Flux[item][Index]

    # plt.figure(figsize=(7, 7))
    # plt.plot(Time['JD'], Mag)
    # plt.plot(Time['JD'], Mag[:, 0], 'r.', label='target')
    # plt.ylabel('Mag', fontsize=6)
    # plt.title('Bad frames deleted', fontsize=8)
    # plt.axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
    # plt.xticks(locs, x_ticks_labels, rotation='vertical', fontsize=6)
    # plt.xlabel('Date-Time (UTC), ' + Time['DATE-OBS'][0].split('T')[0], fontsize=6)
    # #  #plt.xlim(locs[1], locs[-2])
    # plt.tick_params(axis='both', labelsize=6, direction='in')
    # plt.legend()
    # plt.grid()
    # plt.show()

    #  #######################################################################
    #  #####delete bad frames
    #  ###Index = np.where(Airmass[:, 2]<4000)[0]
    #  #####print(Index)
    #  ###Mag = Mag[Index]
    #  ###Error = Error[Index]
    #  ###Time =  Time[Index]
    #  ###Flux = Flux[Index]
    #  ###Airmass = Airmass[Index]
    #  ###plt.plot(Mag)
    #  ###plt.plot(Mag[:, 0], 'r.')
    #  ###plt.show()
    #  #
    #  #
    #  ###
    #  #######search bad frames
    #  #####for ii in range(0, Mag.shape[1]):
    #  #####    Bad_threshold = 0.10
    #  #####    diff = Mag[:, ii] - np.nanmedian(Mag[:, ii], 0)#[:, np.newaxis]
    #  #####    diff = diff - np.nanmedian(diff, 0)
    #  #######    plt.plot(diff)
    #  #######    plt.hlines(Bad_threshold, 0, 800)
    #  #######    plt.show()
    #  #####
    #  #####    ##delete bad frames
    #  #####    Index = np.where(diff<Bad_threshold)[0]
    #  #######    print(Index)
    #  #####    Mag = Mag[Index]
    #  #####    Error = Error[Index]
    #  #####    Time =  Time[Index]
    #  #####    Flux = Flux[Index]
    #  #####    Airmass = Airmass[Index]
    #  #######    plt.plot(Mag)
    #  #######    plt.plot(Mag[:, 0], 'r.')
    #  #######    plt.show()
    #  ###
    #  #####Mag = np.delete(Mag, 4, axis=1)

    ########################################################################

    for item in range(len(FluxFiles)):
        flag = 0
        Merror = np.nanmedian(Error[item], 0)
        while flag != 1:
            #  #delete strong trend
            TFlux = Flux[item][:, 1:Flux[item].shape[1]]
            diff = TFlux / np.nanmedian(TFlux, 0)  # [:, np.newaxis]
            diff = diff / np.nanmedian(diff, 0)
            trend = sigma_clipped_stats(diff, sigma=3,
                                        cenfunc=np.nanmedian, stdfunc=np.nanstd, axis=1)[1]
            # plt.plot(diff)
            # plt.plot(trend, 'r--')
            # plt.xlabel('#Frame')
            # plt.ylabel('Relative flux')
            # plt.title('Main trend')
            # plt.show()

            Flux[item] = Flux[item] / trend[:, np.newaxis]
            NEB_Flux[item] = NEB_Flux[item] / trend[:, np.newaxis]
            # plt.plot(Flux[:, 1:])
            # plt.plot(Flux[:, :1], 'r.')
            # plt.xlabel('#Frame')
            # plt.ylabel('Flux')
            # plt.grid()
            # plt.show()

            #  #delete bad stars
            Std = np.nanstd(Flux[item][:, 1:], 0) / (np.nanmedian(Flux[item][:, 1:], 0) * Merror[1:])
            if np.max(Std) > 2:
                Index = np.argmax(Std) + 1
                # print('Delete star #', Index, ' with STD =', np.max(Std))
                Flux[item] = np.delete(Flux[item], Index, axis=1)
                SN[item] = np.delete(SN[item], Index, axis=1)
                Mag[item] = np.delete(Mag[item], Index, axis=1)
                Error[item] = np.delete(Error[item], Index, axis=1)
                Sky[item] = np.delete(Sky[item], Index, axis=1)
                Max[item] = np.delete(Max[item], Index, axis=1)
                Merror = np.delete(Merror, Index)
                Cat[item].remove_row(Index)
            else:
                flag = 1

    ########################################################################
    #   check result
    #   plot std vs mag
    M = [20. - 2.5 * np.log10(item) for item in Flux]
    S = [np.nanstd(m, axis=0) for m in M]
    M = [np.nanmedian(m, axis=0) for m in M]
    # M = [m + (Info['V'][0] - m[0]) for m in M]
    M = [m + (10.582 - m[0]) for m in M]
    # plt.plot(M, S, '.')
    # plt.plot(M[0], S[0], 'r*', label='target')
    # plt.ylabel('std(mag)', fontsize=6)
    # plt.xlabel('~Vmag', fontsize=6)
    # plt.tick_params(axis='both', labelsize=6, direction='in')
    # plt.legend()
    # plt.title('STD vs Mag')
    # plt.grid()
    # plt.show()

    #  plot stars
    s_flux = [item[:, 1:] / np.nanmedian(item[:, 1:], axis=0) for item in Flux]
    s_sn = [item[:, 1:] for item in SN]
    id_s = [cat['ID'][1:] for cat in Cat]
    _STD = [np.nanstd(s_f, 0) for s_f in s_flux]
    sort = [np.argsort(_std) for _std in _STD]
    s_flux = [s_f[:, sort[num]] for num, s_f in enumerate(s_flux)]
    s_sn = [s_sn[:, sort[num]] for num, s_sn in enumerate(s_sn)]
    id_s = [id_s[num][sort] for num, sort in enumerate(sort)]

    ########################################################################
    #  start report
    Title = 'TIC229510866.01'
    Title += '\n'
    Title += Scope + ', ' + Time[0]['DATE-OBS'][0].split('T')[0]
    Fil = Time[0]['FILTER'][0]
    Title += ', filter=' + Fil
    Exp = Time[0]['EXPTIME'][0]
    Title += ', Extime=' + '{:.1f}'.format(Exp)
    Title += '\n'
    Title += 'Observation start=' + Time[0]['DATE-OBS'][0].split('.')[0]
    Title += ', '
    Title += 'end=' + Time[0]['DATE-OBS'][-1].split('.')[0]
    Title += ', ' + str(np.round((Time[0]['JD'][-1] - Time[0]['JD'][0]) * 24 * 60, 0)) + ' min.'
    # Title += '\nTransit start=' + Info['start time'][0].replace(' ', 'T')
    # Title += ', '
    # Title += 'end=' + Info['end time'][0].replace(' ', 'T')
    Title += ', ' + str(len(Time[0]['JD'])) + ' frames'
    Report_Name = Path2Data + '/' + Pref + '_lightcurve.pdf'

    #  zero level and dip
    Depth = 1 / 10 ** (7.246375 / 2500)  # ?????????????? ???????????????? ???? ?????????????????????? ???????????????????? TODO WARNING
    # Depth = 1 / 10 ** (Info['depth(mmag)'][0] / 2500)  # ?????????????? ???????????????? ???? ?????????????????????? ????????????????????
    Target_Report_Tbl = Time
    model = 'power-2'
    model_params = []
    for item in range(len(Flux)):
        Index = range(0, 10)
        # Index = np.where((Time[item]['JD'] < T0.jd) | (Time[item]['JD'] > T1.jd))[0]
        Zero_Flux = sigma_clipped_stats(Flux[item][Index, :],
                                        sigma=3, cenfunc=np.nanmedian, stdfunc=np.nanstd, axis=0)[1]
        Target_Report_Tbl[item].remove_columns(['EXPTIME', 'FILTER', 'FLAG'])
        Target_Report_Tbl[item]['X_Shift'] = np.round(Target_Report_Tbl[item]['X_Shift'] -
                                                      np.mean(Target_Report_Tbl[item]['X_Shift']), 2)
        Target_Report_Tbl[item]['Y_Shift'] = np.round(Target_Report_Tbl[item]['Y_Shift'] -
                                                      np.mean(Target_Report_Tbl[item]['Y_Shift']), 2)
        Target_Report_Tbl[item]['Rel_Flux_Target'] = np.round(Flux[item][:, 0] / Zero_Flux[0], 5)
        Target_Report_Tbl[item]['Rel_Flux_Err_Target'] = np.round(1 / SN[item][:, 0], 5)

        # fitting
        model_params.append(do_fitting(Time[0]['JD'][0], Time[0]['JD'][-1], model, Target_Report_Tbl[item]['BJD_TDB'],
                                       Target_Report_Tbl[item]['Rel_Flux_Target'],
                                       2.103193, (1 - Depth) * 1000, Fil))
        # f, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3 = \
        #     do_fitting(T0, T1, model, Target_Report_Tbl['JD'],
        #                Target_Report_Tbl['Rel_Flux_Target'], Info['period(days)'][0], (1 - Depth) * 1000)

    # bin
    min_bin_size = 3
    bin_sizes = [int(len(Time[item]['DATE-OBS']) / 10) for item in range(len(Flux))]
    time_series = []
    ts_binned = []
    ts_std = []
    transit_depth = []
    std_transit_depth = []
    transit_depth_median = []
    for item_num, bin_size in enumerate(bin_sizes):

        t1 = model_params[item_num][9]
        t2 = model_params[item_num][10]
        if bin_size > 10:
            bin_size = 10
        if bin_size >= min_bin_size:
            time_series.append(TimeSeries(time=Time[item_num]['DATE-OBS']))
            time_series[-1]['flux'] = Target_Report_Tbl[item_num]['Rel_Flux_Target']
            ts_binned.append(aggregate_downsample(time_series=time_series[-1], time_bin_size=bin_size * u.min,
                                                  aggregate_func=np.nanmean))
            ts_std.append(aggregate_downsample(time_series=time_series[-1], time_bin_size=bin_size * u.min,
                                               aggregate_func=np.nanstd))
            transit_depth_index = [item for item in range(len(ts_binned[-1]))
                                   if Target_Report_Tbl[item_num]['JD'][t2]  # t2
                                   > ts_binned[-1].time_bin_start.jd[item]
                                   > Target_Report_Tbl[item_num]['JD'][t1]]  # t1
            transit_depth.append(ts_binned[-1]['flux'][transit_depth_index[0]:transit_depth_index[-1]])
        else:
            transit_depth.append(Target_Report_Tbl[item_num]['Rel_Flux_Target'][t1:t2])
        depth_sigma_clipped = sigma_clipped_stats(transit_depth[-1], sigma=3,
                                                  cenfunc=np.nanmedian, stdfunc=np.nanstd, axis=0)
        std_transit_depth.append(depth_sigma_clipped[2])
        if bin_size >= min_bin_size:
            std_transit_depth[-1] = np.sqrt(std_transit_depth[-1] ** 2 + np.nansum(np.square(ts_std[-1]['flux'])))
        transit_depth_median.append(depth_sigma_clipped[1])

    #  draw light curves
    fig, axs = plt.subplots(5, 1, figsize=(6, 12), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
    fig.suptitle(Title, fontsize=8)
    ZERO = int(Target_Report_Tbl[0]['BJD_TDB'][0])
    print('zero =', ZERO)
    for num in range(len(Flux)):
        # print(Target_Report_Tbl[num]['BJD_TDB'] - ZERO)
        f, k_model, a, mid_time, i_model, e, w, ldc1, ldc2, t1, t2, t0, t3 = model_params[num]
        pos = [0.125, 0.8 - num * 0.156, 0.8, 0.12]
        axs[num].set_position(pos)
        axs[num].plot(Target_Report_Tbl[num]['BJD_TDB'] - ZERO, f, 'y.', markersize=3, zorder=4,
                      linewidth=0.5, label='model, Limb darkening: ' + model + ', k = ' + '{:.2f}'.format(k_model) +
                                           ', a = ' + '{:.1f}'.format(a) + ', i = ' +
                                           '{:.2f}'.format(i_model) + 'Pi, '
                                           + 'e = ' + '{:.1f}'.format(e) + ', w = ' + '{:.2f}'.format(w) + 'Pi',
                      linestyle='solid')  # ???????????? ????????????
        axs[num].errorbar(Target_Report_Tbl[num]['BJD_TDB'] - ZERO,
                          Target_Report_Tbl[num]['Rel_Flux_Target'],
                          Target_Report_Tbl[num]['Rel_Flux_Err_Target'], fmt='b.', zorder=3,
                          label='TIC 229510866.01' + r', 1$\sigma$ errorbar,' +
                                                  ' depth = ' + '{:.3f}'.format(transit_depth_median[num]) + '??' +
                                                  '{:.3f}'.format(std_transit_depth[num]))  # ???????????? ????????????
        if bin_sizes[num] >= min_bin_size:
            axs[num].errorbar(x=ts_binned[num].time_bin_start.tdb.jd - ZERO,
                              y=ts_binned[num]['flux'],
                              yerr=ts_std[num]['flux'],
                              fmt='g.', markersize=9, zorder=4,
                              label='x' + str(bin_sizes[num]) + ' binned' + r', 1$\sigma$ errorbar')  # ???????????? ????????

        axs[num].legend(loc=2, fontsize=6)
        axs[num].set_ylabel('Normalized flux', fontsize=6)
        # axs[num].axvspan(T0.tdb.jd - ZERO, T1.tdb.jd - ZERO, facecolor='k', alpha=0.2)  # ?????????????????????????????????? ????????????????
        # ???? ?????????????????????? ????????????????????
        locs = axs[0].get_xticks()
        t = aTime(locs, format='jd')
        x_ticks_labels = []
        for x in t:
            # print(str(round(x.tdb.jd - ZERO, 2)))
            x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
        axs[num].set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
        axs[num].set_xlabel('BJD_TDB - ' + str(ZERO), fontsize=6)

        axs[num].axhline(Depth, linewidth=0.5, color='r')  # ?????????????? ???????????????? ???? ?????????????????????? ????????????????????
        axs[num].axhline(transit_depth_median[num], linewidth=0.5, color='g')  # ?????????????? ???????????????? ???? ????????????????????
        draw_my_annotate(axs[num], Target_Report_Tbl[num]['BJD_TDB'][t0] - ZERO)  # ???????????? ???????????????? ???? ????????????????????
        draw_my_annotate(axs[num], Target_Report_Tbl[num]['BJD_TDB'][t1] - ZERO)  # ??????1
        draw_my_annotate(axs[num], Target_Report_Tbl[num]['BJD_TDB'][t2] - ZERO)  # ??????2
        draw_my_annotate(axs[num], Target_Report_Tbl[num]['BJD_TDB'][t3] - ZERO)  # ?????????? ???????????????? ???? ????????????????????

        axs[num].tick_params(axis='both', labelsize=6, direction='in')
        axs[num].set_title('Target, aperture radius - ' +
                           '{:.1f}'.format(Time[num]['APER' + str(1 + num)][0]) + 'pix', loc='left', fontsize=6)
        axs[num].grid()

    #  plot stars
    pos = [0.125, 0.23, 0.8, 0.22]
    axs[len(Flux)].set_position(pos)
    _k = 0
    axs[len(Flux)].plot(Time[0]['BJD_TDB'] - ZERO, Flux[0][:, 0] / np.median(Flux[0][:, 0]), 'r*', label='Target')
    for _i in range(0, len(sort[0])):
        if _k < 7:
            axs[len(Flux)].plot(Time[0]['BJD_TDB'] - ZERO, s_flux[0][:, _i] - 0.05 * (_i + 1),
                                '.', label='Ref#' + str(id_s[0][_i]))
        _k = _k + 1
    axs[len(Flux)].legend(loc=2, fontsize=6)
    axs[len(Flux)].set_ylabel('Normalized flux', fontsize=6)
    # axs[len(Flux)].axvspan(T0.tdb.jd - ZERO, T1.tdb.jd - ZERO, facecolor='k', alpha=0.2)
    locs = axs[len(Flux)].get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        x_ticks_labels.append(str('{:.3f}'.format(x.tdb.jd)))
    axs[len(Flux)].set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    axs[len(Flux)].set_xlabel('BJD_TDB - ' + str(ZERO), fontsize=6)
    axs[len(Flux)].tick_params(axis='both', labelsize=6, direction='in')
    axs[len(Flux)].set_title('Reference stars, aperture radius - ' +
                             '{:.1f}'.format(Time[0]['APER1'][0]) + 'pix', loc='left', fontsize=6)
    axs[len(Flux)].grid()

    #  plot errors
    pos = [0.125, 0.03, 0.8, 0.16]
    axs[len(Flux) + 1].set_position(pos)
    axs[len(Flux) + 1].plot(M[0], S[0], '.')
    axs[len(Flux) + 1].plot(M[0][0], S[0][0], 'r*', label='target')
    axs[len(Flux) + 1].set_ylabel('std(mag)', fontsize=6)
    axs[len(Flux) + 1].set_xlabel('~Vmag', fontsize=6)
    axs[len(Flux) + 1].tick_params(axis='both', labelsize=6, direction='in')
    axs[len(Flux) + 1].legend(loc=2, fontsize=6)
    axs[len(Flux) + 1].set_title('Errors vs magnitudes, aperture radius - ' +
                                 '{:.1f}'.format(Time[0]['APER1'][0]) + 'pix', loc='left', fontsize=6)
    axs[len(Flux) + 1].grid()

    plt.savefig(Report_Name)
    # plt.show()

    #  NEB
    fig1, ax = plt.subplots(1, 1, figsize=(6, 7), dpi=125)
    fig1.suptitle(Title, fontsize=8)
    #  plot stars
    pos = [0.125, 0.07, 0.8, 0.8]
    ax.set_position(pos)
    _k = 0
    NEB_Flux[0] = NEB_Flux[0] / np.median(NEB_Flux[0], axis=0)  # TODO ???????????????????? ??????????????????
    #                                                              ?????????????? ???????????? ?? ???????? ???????????? (?????????????)

    # NEB_Flux = [NEB_Flux[:, q] for q in range(1, NEB_Flux) if np.nanstd(NEB_Flux[:, q]) < 2]
    ax.plot(Time[0]['BJD_TDB'] - ZERO, NEB_Flux[0][:, 0], 'r*',
            label='Target, RMS=' + str(np.round(np.nanstd(NEB_Flux[0][:, 0]), 3)))
    for _i in range(1, NEB_Flux[0].shape[0]):
        try:
            # print('std', np.nanstd(NEB_Flux[0][:, _i]))
            if np.nanstd(NEB_Flux[0][:, _i]) > 3:
                continue
        except:
            pass
        if _k < 10:
            ax.plot(Time[0]['BJD_TDB'] - ZERO, NEB_Flux[0][:, _i] - 0.05 * _i,
                    '.', label='Ref#' + str(_i) + ',' + 'Dist=' +
                               str(np.round(Cat[0]['Dist'][_i] * 3600., 1)) + '", ' +
                               'RMS=' + str(np.round(np.nanstd(NEB_Flux[0][:, _i]), 3)))

            _k = _k + 1
    ax.legend(loc=2, fontsize=6)
    ax.set_ylabel('Normalized flux', fontsize=6)
    # ax.axvspan(T0.tdb.jd - ZERO, T1.tdb.jd - ZERO, facecolor='k', alpha=0.2)
    locs = ax.get_xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        x_ticks_labels.append(str('{:.3f}'.format(x.tdb.jd)))
    ax.set_xticklabels(x_ticks_labels, fontsize=5)  # rotation='vertical',
    ax.set_xlabel('BJD_TDB - ' + str(ZERO), fontsize=6)
    ax.tick_params(axis='both', labelsize=6, direction='in')
    ax.set_title(str(_k) + ' nearest stars, sorted by distance from target' +
                 ', aperture radius - ' + '{:.1f}'.format(Time[0]['APER1'][0]), loc='left', fontsize=6)
    ax.grid()
    Report_Name = Path2Data + '/' + Pref + '_NEBs_checking.pdf'
    plt.savefig(Report_Name)
    # plt.show()

    # report
    _k = 0
    for _i in range(0, len(sort)):
        if _k < 3:
            F = 'Rel_Flux_Ref#' + str(id_s[0][_i])
            E = 'Rel_Flux_Err_Ref#' + str(id_s[0][_i])
            Target_Report_Tbl[0][F] = np.round(s_flux[0][:, _i], 5)
            Target_Report_Tbl[0][E] = np.round(1 / s_sn[0][:, _i], 5)
        _k = _k + 1

    ascii.write(Target_Report_Tbl[0], Path2Data + '/' + Pref + '_lightcurve.dat',
                overwrite=True, delimiter='\t', fast_writer=False,
                format='commented_header')


# start = time.time()
# do_post_proc(r'D:\docs\sci\tess\good\for_cursuch\TIC22951086601\3988\TIC229510866.01_20201118_Kourovka_0.4_V')
# print('time =', time.time() - start)
