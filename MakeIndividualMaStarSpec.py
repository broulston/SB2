import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

main_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/SB2_Composites/"
data_dir ="data/"

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3.fits')
mastarSpec = fits.open(main_dir+'mastar-goodspec-v2_4_3-v1_0_2.fits')

MANGAID = mastar[1].data.field('MANGAID')
NVISITS = mastar[1].data.field('NVISITS')
RA = mastar[1].data.field('RA')
DEC = mastar[1].data.field('DEC')

outputDIR = main_dir+"MaStar_spec/"
header = "wavelength, flux, err"

for ii, ID in enumerate(MANGAID):
    all_visits = np.where(mastar[2].data.field('MANGAID') == ID)[0]
    for jj in range(all_visits.size):
        this_visit = all_visits[jj]
        wavelength = mastarSpec[1].data.field('WAVE')[this_visit]
        flux = mastarSpec[1].data.field('FLUX')[this_visit]
        ivar = mastarSpec[1].data.field('IVAR')[this_visit]
        var = 1.0/ivar
        error = np.sqrt(var)
        error[~np.isfinite(error)] = 999.99
        MJD = mastar[2].data.field('MJD')[this_visit]
        ra_string = '{:0>9.5f}'.format(RA[ii])
        dec_string = '{:0=+10.5f}'.format(DEC[ii])
        MANGAID_string = ID
        mjd_string = '{:0>5}'.format(str(np.int(MJD)))
        filename_string = ra_string+dec_string+"_"+mjd_string+"_"+MANGAID_string
        spec = np.stack((wavelength, flux, error), axis=1)
        np.savetxt(outputDIR+filename_string+".txt", spec, delimiter=",", header=header, fmt="%10.5f")
        #print(str(np.round((ii+1)*100./MANGAID.size,2))+"%")