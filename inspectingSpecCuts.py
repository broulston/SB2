import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table

main_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/"
data_dir = "data/"

mastar = fits.open(data_dir + 'mastarall-gaia-v2_4_3_SDSSDR12_specTypes.fits')
mastarSpec = fits.open(main_dir + 'mastar-goodspec-v2_4_3-v1_0_2.fits')

mastar1_Table = Table.read(mastar[1])
mastar2_Table = Table.read(mastar[2])

PHOT_G_MEAN_MAG = mastar[1].data.field('PHOT_G_MEAN_MAG')
BP_RP = mastar[1].data.field('BP_RP')
gmag_SDSSDR12 = mastar[1].data.field('gmag_SDSSDR12')
imag_SDSSDR12 = mastar[1].data.field('imag_SDSSDR12')
R_EST = mastar[1].data.field('R_EST')
M_G = PHOT_G_MEAN_MAG + 5.0 - 5.0 * np.log10(R_EST)

MANGAID = mastar[1].data.field('MANGAID')
NVISITS = mastar[1].data.field('NVISITS')
RA = mastar[1].data.field('OBJRA')
DEC = mastar[1].data.field('OBJDEC')

R_EST = mastar[1].data.field('R_EST') * u.pc

specTypes = mastar[1].data.field('Guessed Spectral Type')
unique_specTypes = np.unique(specTypes)

specType_SB2_combos = np.array(np.meshgrid(unique_specTypes,
                               unique_specTypes)).T.reshape(-1, 2)

orginal_outputDIR = main_dir + "MaStar_spec/"
converted_outputDIR = main_dir + "MaStar_specLum/"
header = "wavelength, flux, err"

all_filenames = np.genfromtxt(data_dir + "list_of_all_MaStar_specCUT.txt",
                              dtype="U")
for ii in range(all_filenames.size):
    spec = np.loadtxt(orginal_outputDIR+all_filenames[ii], delimiter=",")
    wavelength = spec[:, 0]
    flux = spec[:, 1]
    plt.plot(wavelength, flux+ii)
    plt.show()
    plt.clf()
    plt.close()
