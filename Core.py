import os
# disable warnings
import time
import warnings

import astropy.io.fits as fits
from astropy.io import ascii
# import astropy modules
import astropy.wcs as wcs
import numpy as np

# import my modules
from Draw_Map import Draw_Map
from Get_Cat import Get_UCAC, Get_GAIA
from Headers import Get_Info
from PSF2FWHM import PSF2FWHM
from Phot import Phot
# from TTF import TTF_Query
from get_PSF import get_PSF
import Post_Proc
from unzip import unzip

warnings.simplefilter("ignore")


def core(Path2Data, RAper=[4.0, 6.0, 8.0]):
    ######################################################################
    # set paths. catalog and target name (or none for name from header)

    # Path2Coo = 'D:/Temp/new_from_master/GPX/TF8A1859/coo_1.txt'
    Path2Coo = 'GAIA'
    # Target = 'TIC+147950620.01'
    Target = None
    Scope = 'Kourovka_0.4'
    # priority = 5  # set 5 for old data

    ######################################################################
    # set parameters for photometry
    saturation = 45000
    FWHM_e = 2.0
    gain = 1.3
    r_noise = 10.
    v_lim = 15.
    max_mag = 9.
    ######################################################################
    # unzip
    unzip(Path2Data)
    ######################################################################
    # read directory and create list of fits-files
    file_list = []
    dir_content = os.listdir(Path2Data)
    for ii in range(0, len(dir_content)):
        if dir_content[ii].count('.fit') or dir_content[ii].count('.fits'):
            file_list.append(Path2Data + '/' + dir_content[ii])

    # read first frame for object name and other information
    file_name = file_list[0]
    print('First frame: ' + file_name.split('/')[-1])
    hdulist = fits.open(file_name)
    Header = hdulist[0].header
    # Data = hdulist[0].data.copy()
    hdulist.verify('fix')
    hdulist.close()
    Info = Get_Info(Header)

    # get transit info from TTF
    DT = Info[6].datetime.strftime('%m-%d-%Y')
    if Target == None:
        Target = Header['OBJNAME']
        Target = Target[0:3] + Target[3:].replace('_', '.')
    print('Date: ', DT, ' Target: ', Target)
    # TTF = TTF_Query(DT, 6, 0.5, 16, Target)
    # print('TTF: ', TTF)
    # try:
    #     if not str(TTF['r_planet'][0]).isdigit():
    #         TTF['r_planet'][0] = '-1'
    # except:
    #     pass
    # print('TTF, R = ', TTF['r_planet'][0])
    # print('Object: ', Target, '=', TTF['Name'][0])
    # print('Start-end: ', TTF['start time'][0], '-', TTF['end time'][0])
    # print('Depth: ', TTF['depth(mmag)'][0])
    # print('Priority:', TTF['priority'][0])
    # print('Comments: ', TTF['comments'][0])

    # set paths and suffics for output files
    Suff = Header['FILTER']  # filter
    Pref = Target + '_' + Info[6].datetime.strftime('%Y%m%d') + '_' + Scope + '_' + Suff
    # Path2Save = os.path.split(Path2Data)[0] + '/' + Suff + '{0:.1f}'.format(RAper[0]) + '{0:.1f}'.format(RAper[1]) + \
    #                                           '{0:.1f}'.format(RAper[2])  # + '{0:.1f}'.format(RAper[3])
    Path2Save = os.path.split(Path2Data)[0] + '/' + Pref
    if not os.path.exists(Path2Save):
        os.makedirs(Path2Save)

    # save transit info
    # ascii.write(TTF, Path2Save + '/' + 'TTF_Info.txt', fast_writer=False,
    #             overwrite=True, delimiter='\t', format='commented_header')

    # read coo file or create from UCAC or GAIA

    if Path2Coo == 'UCAC':
        Catalog = Get_UCAC(Info[0], Info[1], (Info[2] - 20) * Info[3] / 2., v_lim)
    elif Path2Coo == 'GAIA':
        Catalog = Get_GAIA(Info[0], Info[1], (Info[2] - 20) * Info[3] / 2., v_lim)
    else:
        Catalog = np.genfromtxt(Path2Coo, skip_header=1)

    # save catalog
    ascii.write(Catalog, Path2Save + '/Coo.txt', fast_writer=False,
                overwrite=True, delimiter='\t', format='commented_header')

    # open log file
    df = open(Path2Save + '/Time.txt', 'w')
    df.write('DATE-OBS\tJD\tHJD\tBJD_TDB\t' +
             'EXPTIME\tFILTER\tAPER1\t' +
             'APER2\tAPER3\tAIRMASS\t' +
             'EXTINC\tFWHM\tELL\t' +
             'X_Shift\tY_Shift\t' +
             'Sky\tFLAG\n')

    Flux = []
    SN = []
    Mag = []
    MErr = []
    Sky = []
    Max = []

    counter = len(file_list)
    for ii in range(0, len(file_list)):
        print('----------<>-----------')
        counter = counter - 1
        print('Frames: ' + str(counter))
        print(file_list[ii].split('/')[-1])

        hdulist = fits.open(file_list[ii])
        Header = hdulist[0].header
        Data = hdulist[0].data
        hdulist.verify('fix')
        hdulist.close()
        Info = Get_Info(Header)
        T = Info[6].datetime.isoformat(timespec='milliseconds')

        # make WCS object
        w = wcs.WCS(Header)
        XY = w.all_world2pix(Catalog['Ra'], Catalog['Dec'], 0)
        XY = np.vstack((XY[0], XY[1])).T
        #     print(Header)
        #     print(w.wcs.cd)

        # # find fwhm
        # if is_fwhm_to_r_aper:
        #     FWHM = find_fwhm(Data, XY=XY[0])
        #     r_aper_from_fwhm.append(FWHM_Scale*FWHM)
        #     print('{:.1f}'.format(FWHM_Scale) + 'FWHM = ', r_aper_from_fwhm[-1])

        # make subarray of PSF stars
        Index = np.where((Catalog['R'] > max_mag) & (Catalog['R'] < max_mag + 5))[0]
        PSF_stars = XY[Index, :]
        print(str(len(PSF_stars)) + ' stars use for PSF model')
        # making of PSF model
        PSF_model = get_PSF(Data, PSF_stars, FWHM_e)
        FWHM1, FWHM2, Phi, PSky = PSF2FWHM(PSF_model)
        FWHM = (FWHM1 + FWHM2) / 2.
        print('FWHM1={0:.2f}'.format(FWHM1), '\tFWHM2={0:.2f}'.format(FWHM2),
              '\tPhi={0:.2f}'.format(Phi * 180. / np.pi),
              '\tSky={0:.1f}'.format(PSky))

        if ii == 0:
            # draw picture and save
            Pic_Name = Path2Save + '/' + Pref + '_field.pdf'
            Draw_Map(Data, Header, Catalog, Pic_Name, RAper[0], Target, Info, XY)

        df.write(T + '\t' + '{:.7f}'.format(Info[6].jd) + '\t')  # DATE-OBS JD
        df.write('{:.7f}'.format(Info[7]) + '\t' + '{:.7f}'.format(Info[8]) + '\t')  # HJD BJD
        df.write('{:.1f}'.format(Header['EXPTIME']) + '\t')  # EXPTIME
        df.write(Header['FILTER'] + '\t')  # FILTER
        for r in RAper:
            df.write('{:.2f}'.format(r) + '\t')  # Aperture
        df.write('{:.2f}'.format(Header['AIRMASS']) + '\t')  # Airmass
        df.write('{:.2f}'.format(Header['EXTINCT']) + '\t')  # EXTINC
        df.write('{:.2f}'.format(Header['SEXFWHM']) + '\t')  # FWHM
        df.write('{:.2f}'.format(Header['SEXELL']) + '\t')  # Ell
        df.write(str(Info[4]) + '\t' + str(Info[5]) + '\t')  # X_Shift Y_Shift

        df.write('{:.2f}'.format(PSky) + '\t')  # Sky

        # photometry Flux, SN, Mag, MErr, Sky, Max
        Phot_Result = Phot(Data, XY, PSF_model, RAper, FWHM, gain, r_noise, saturation)

        Flux.append(Phot_Result[0])
        SN.append(Phot_Result[1])
        Mag.append([item + 2.5 * np.log10(Header['EXPTIME']) for item in Phot_Result[2]])
        MErr.append(Phot_Result[3])
        Sky.append(Phot_Result[4])
        Max.append(Phot_Result[5])

        if np.isnan(Phot_Result[2][0][0]):
            df.write('ERR' + '\n')
        elif Phot_Result[5][0] > saturation:
            df.write('SAT' + '\n')
        else:
            df.write('OK' + '\n')

    df.close()

    Flux = np.array(Flux)
    SN = np.array(SN)
    Mag = np.array(Mag)
    MErr = np.array(MErr)
    Sky = np.array(Sky)
    Max = np.array(Max)

    np.savetxt(Path2Save + '/Flux1.txt', Flux[:, 0], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/Flux2.txt', Flux[:, 1], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/Flux3.txt', Flux[:, 2], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/SN1.txt', SN[:, 0], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/SN2.txt', SN[:, 1], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/SN3.txt', SN[:, 2], fmt='%.1f', delimiter='\t')
    np.savetxt(Path2Save + '/Mag1.txt', Mag[:, 0], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Mag2.txt', Mag[:, 1], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Mag3.txt', Mag[:, 2], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Merr1.txt', MErr[:, 0], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Merr2.txt', MErr[:, 1], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Merr3.txt', MErr[:, 2], fmt='%.4f', delimiter='\t')
    np.savetxt(Path2Save + '/Sky.txt', Sky, fmt='%.2f', delimiter='\t')
    np.savetxt(Path2Save + '/Max.txt', Max, fmt='%.1f', delimiter='\t')

    # np.savetxt(Path2Save + '/Flux4.txt', Flux[:, 3], fmt='%.1f', delimiter='\t')
    # np.savetxt(Path2Save + '/SN4.txt', SN[:, 3], fmt='%.1f', delimiter='\t')
    # np.savetxt(Path2Save + '/Mag4.txt', Mag[:, 3], fmt='%.4f', delimiter='\t')
    # np.savetxt(Path2Save + '/Merr4.txt', MErr[:, 4], fmt='%.4f', delimiter='\t')

    Post_Proc.do_post_proc(Path2Save)


start = time.time()
core(Path2Data=r'D:\docs\sci\tess\good\TIC22951086601\4012\EAST_V')
print('time =', time.time() - start)
