import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

main_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/SB2_Composites/"
data_dir ="data/"

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3.fits')
mastarSpec = fits.open(main_dir+"mastar-goodspec-v2_4_3-v1_0_2.fits")

MANGAID = mastar[1].data.field('MANGAID')
NVISITS = mastar[1].data.field('NVISITS')
RA = mastar[1].data.field('RA')
DEC = mastar[1].data.field('DEC')

PARALLAX = mastar[1].data.field('PARALLAX')
PARALLAX_ERROR = mastar[1].data.field('PARALLAX_ERROR')
PARALLAX_SNR = np.abs(PARALLAX / PARALLAX_ERROR)
R_EST = mastar[1].data.field('R_EST') * u.pc

PyHammerResult = np.genfromtxt(data_dir+"PyHammerResults_MaStar.csv", dtype='str', delimiter=",")
PyHammerResultCondensed = PyHammerResult[:,0:4:3]

for ii in range(PyHammerResultCondensed[:,0].size):
    PyHammerResultCondensed[ii,0] = PyHammerResultCondensed[ii,0].replace("/", " ").split()[-1][:-4]

outputDIR = main_dir+"MaStar_specLum/"
header = "wavelength, Lum [W], err"
all_Lums = []
all_filenames = []
for ii, ID in enumerate(MANGAID):
    all_visits = np.where(mastar[2].data.field('MANGAID') == ID)[0]
    for jj in range(all_visits.size):
        this_visit = all_visits[jj]
        MJD = mastar[2].data.field('MJD')[this_visit]
        ra_string = '{:0>9.5f}'.format(RA[ii])
        dec_string = '{:0=+10.5f}'.format(DEC[ii])
        MANGAID_string = ID
        mjd_string = '{:0>5}'.format(str(np.int(MJD)))
        filename_string = ra_string+dec_string+"_"+mjd_string+"_"+MANGAID_string

        if np.isin(filename_string, PyHammerResultCondensed[:,0]):
            wavelength = mastarSpec[1].data.field('WAVE')[this_visit] * u.angstrom 
            flux = mastarSpec[1].data.field('FLUX')[this_visit] * u.erg / u.s /u.cm**2 / u.angstrom 
            ivar = mastarSpec[1].data.field('IVAR')[this_visit] 
            var = 1.0/ivar
            error = np.sqrt(var)
            error[~np.isfinite(error)] = 999.99
            error = error * u.erg / u.s /u.cm**2 / u.angstrom 
            Lum = 4.0*np.pi*1e-17 * (flux * wavelength * R_EST[ii].to(u.cm)**2 )
            Lum = Lum.to(u.W)
            Lum_error = 4.0*np.pi*1e-17 * (error * wavelength * R_EST[ii].to(u.cm)**2 )
            Lum_error = Lum_error.to(u.W)
            all_Lums.append(Lum.sum()/const.L_sun)
            #print(Lum.sum()/const.L_sun)
            spec = np.stack((wavelength, Lum, Lum_error), axis=1)
            all_filenames.append(filename_string) 
            np.savetxt(outputDIR+filename_string+".txt", spec, delimiter=",", header=header, fmt="%5.4f, %0.9e, %0.9e")

        print(str(np.round((ii+1)*100./MANGAID.size,2))+"%")

np.savetxt(data_dir+"all_MaStar_spec.txt", all_filenames, delimiter=",", fmt="%s")