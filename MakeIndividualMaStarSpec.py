import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from progressbar import ProgressBar

main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
data_dir = "data/"

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3.fits')
mastarSpec = fits.open(main_dir+'mastar-goodspec-v2_4_3-v1_0_2.fits')

MANGAID = mastar[1].data.field('MANGAID')
NVISITS = mastar[1].data.field('NVISITS')
RA = mastar[1].data.field('OBJRA')
DEC = mastar[1].data.field('OBJDEC')

R_EST = mastar[1].data.field('R_EST') * u.pc

orginal_outputDIR = main_dir+"MaStar_spec/"
converted_outputDIR = main_dir+"MaStar_specLum/"
header = "wavelength, flux, err"

all_Lums = []
all_filenames = []
progress = ProgressBar(MANGAID.size, fmt=ProgressBar.FULL)
for ii, ID in enumerate(MANGAID):
    progress.current += 1
    progress()
    all_visits = np.where(mastar[2].data.field('MANGAID') == ID)[0]
    for jj in range(all_visits.size):
        this_visit = all_visits[jj]
        wavelength = mastarSpec[1].data.field('WAVE')[this_visit] * u.angstrom
        flux = (mastarSpec[1].data.field('FLUX')[this_visit]
                * u.erg / u.s / u.cm**2 / u.angstrom)
        ivar = mastarSpec[1].data.field('IVAR')[this_visit]
        with np.errstate(divide='ignore'):
            var = 1.0/ivar
        error = np.sqrt(var)
        error[~np.isfinite(error)] = 999.99
        error = error * u.erg / u.s / u.cm**2 / u.angstrom
        MJD = mastar[2].data.field('MJD')[this_visit]
        ra_string = '{:0>9.5f}'.format(RA[ii])
        dec_string = '{:0=+10.5f}'.format(DEC[ii])
        MANGAID_string = ID
        mjd_string = '{:0>5}'.format(str(np.int(MJD)))
        filename_string = (ra_string + dec_string + "_"
                           + mjd_string + "_" + MANGAID_string)
        spec = np.stack((wavelength, flux, error), axis=1)
        np.savetxt(orginal_outputDIR+filename_string+".txt", spec,
                   delimiter=",", header=header, fmt="%5.4f, %0.9e, %0.9e")
        Lum = 4.0*np.pi*1e-17 * (flux * wavelength * R_EST[ii].to(u.cm)**2)
        Lum = Lum.to(u.W)
        Lum_error = (4.0*np.pi*1e-17
                     * (error * wavelength * R_EST[ii].to(u.cm)**2))
        Lum_error = Lum_error.to(u.W)
        all_Lums.append(Lum.sum()/const.L_sun)
        # print(Lum.sum()/const.L_sun)
        spec = np.stack((wavelength, Lum, Lum_error), axis=1)
        all_filenames.append(filename_string)
        np.savetxt(converted_outputDIR+filename_string+".txt", spec,
                   delimiter=",", header=header, fmt="%5.4f, %0.9e, %0.9e")

progress.done()

all_filenames_array = np.array(all_filenames)

np.savetxt(data_dir+"list_of_all_MaStar_spec.txt", all_filenames_array,
           delimiter=",", fmt="%s")