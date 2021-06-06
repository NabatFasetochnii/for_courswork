import numpy as np
from math import isnan

from Master_Processing.PSF2FWHM import PSF2FWHM
from Master_Processing.get_PSF import centroid, D2_moffat_fitter


def find_fwhm(data, XY):
    FWHM_e = 2.0
    R1 = np.ceil(FWHM_e)  # estimation of aperture radii
    R2 = np.ceil(FWHM_e * 3.)  # sky annulus inner radii
    R3 = np.ceil(FWHM_e * 6.)  # sky annulus outer radii

    x_coo = XY[0]
    y_coo = XY[1]
    #         print(x_coo, y_coo)
    ROI = np.copy(data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])  # copy small area
    offset = centroid(R1, R2, R3, ROI)  # search centroid, Gauss sigma and mean sky

    if isnan(offset[0]) == False and isnan(offset[1]) == False:
        x_coo = x_coo + offset[0]
        y_coo = y_coo + offset[1]
        MSKY = offset[2]
        ROI = np.copy(data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])
        param = D2_moffat_fitter(ROI, MSKY, x_coo, y_coo, R3)  # fit 2D moffat psf
        if param != None:
            FWHM1, FWHM2, Phi, PSky = PSF2FWHM(param)
            FWHM = (FWHM1+FWHM2)/2.
            return FWHM
    else:
        pass
