import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from progressbar import ProgressBar

main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
data_dir = "data/"

mastar = fits.open(data_dir + 'mastarall-gaia-v2_4_3_SDSSDR12.fits')
mastarSpec = fits.open(main_dir + 'mastar-goodspec-v2_4_3-v1_0_2.fits')

PHOT_G_MEAN_MAG = mastar[1].data.field('PHOT_G_MEAN_MAG')
BP_RP = mastar[1].data.field('BP_RP')
gmag_SDSSDR12 = mastar[1].data.field('gmag_SDSSDR12')
imag_SDSSDR12 = mastar[1].data.field('imag_SDSSDR12')
R_EST = mastar[1].data.field('R_EST')
M_G = PHOT_G_MEAN_MAG + 5.0 - 5.0 * np.log10(R_EST)

cut_index = np.where((gmag_SDSSDR12 > 15.5)
                    & (imag_SDSSDR12 > 15.5) & (M_G > 5.0))[0]


MANGAID = mastar[1].data.field('MANGAID')[cut_index]
NVISITS = mastar[1].data.field('NVISITS')[cut_index]
RA = mastar[1].data.field('OBJRA')[cut_index]
DEC = mastar[1].data.field('OBJDEC')[cut_index]

R_EST = mastar[1].data.field('R_EST')[cut_index] * u.pc

orginal_outputDIR = main_dir + "MaStar_spec_cuts/"
converted_outputDIR = main_dir + "MaStar_specLum_cuts/"
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
        wavelength = (mastarSpec[1].data.field('WAVE')[this_visit]
                      * u.angstrom)
        flux = (mastarSpec[1].data.field('FLUX')[this_visit]
                * u.erg / u.s / u.cm**2 / u.angstrom)
        # wavelength[np.where(flux == 0.0)[0]] = np.nan
        ivar = mastarSpec[1].data.field('IVAR')[this_visit]
        with np.errstate(divide='ignore'):
            var = 1.0 / ivar
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
        # plt.plot(wavelength, flux)
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
        all_filenames.append(filename_string+".txt")
        np.savetxt(converted_outputDIR+filename_string+".txt", spec,
                   delimiter=",", header=header, fmt="%5.4f, %0.9e, %0.9e")

np.savetxt(data_dir+"list_of_all_MaStar_specCUT.txt", all_filenames,
           delimiter=",", fmt="%s")
progress.done()
# plt.show()
# plt.clf()
# plt.close()
