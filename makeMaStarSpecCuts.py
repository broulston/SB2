import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
from progressbar import ProgressBar
from subprocess import *
import os
import shutil


main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
data_dir = "/Users/benjaminroulston/Dropbox/GitHub/SB2/data/"

mastar = fits.open(data_dir + 'mastarall-gaia-v2_4_3_specTypes_PS1_SDSSDR12_ALLFILE.fits')

RA = mastar[1].data.field("OBJRA")
DEC = mastar[1].data.field("OBJDEC")

# **********************************************************************************
# **********************************************************************************
c_icrs = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
c_gal = c_icrs.galactic 

col1 = fits.Column(name='b', format='D',
                   array=c_gal.b.value)
col2 = fits.Column(name='l', format='D',
                   array=c_gal.l.value)
coldefs = fits.ColDefs([col1, col2])
hdu = fits.BinTableHDU.from_columns(coldefs)

new_columns = mastar[1].columns + hdu.columns
new_hdu = fits.BinTableHDU.from_columns(new_columns, name="ALLFILE")

new_HDUList = fits.HDUList(hdus=[mastar[0], new_hdu])
new_HDUList.writeto(data_dir+"mastarall-gaia-v2_4_3_specTypes_PS1_SDSSDR12_ALLFILE.fits", overwrite=True)
# **********************************************************************************
# **********************************************************************************
mastar = fits.open(data_dir + 'mastarall-gaia-v2_4_3_specTypes_PS1_SDSSDR12_ALLFILE.fits')

RA = mastar[1].data.field("OBJRA")
DEC = mastar[1].data.field("OBJDEC")

PHOT_G_MEAN_MAG = mastar[1].data.field('PHOT_G_MEAN_MAG')
BP_RP = mastar[1].data.field('BP_RP')
R_EST = mastar[1].data.field('R_EST')
M_G = PHOT_G_MEAN_MAG + 5.0 - 5.0 * np.log10(R_EST)

gmag_SDSSDR12 = mastar[1].data.field('gmag_SDSSDR12')
imag_SDSSDR12 = mastar[1].data.field('imag_SDSSDR12')

b = mastar[1].data.field('b') * u.deg

# cut_index = np.where((gmag_SDSSDR12 > 15.5)
#                      & (imag_SDSSDR12 > 15.5)
#                      & (c_gal.b > 10.0*u.deg)
#                      & (M_G > 5.0)
#                      )[0]

cut_index = np.where((b > 10.0*u.deg)
                     & (M_G > 5.0)
                     )[0]


cut_filenames = []
shutil.rmtree(data_dir+"MaStar_spec_cuts/")
shutil.rmtree(data_dir+"MaStar_specLum_cuts/")
if not os.path.exists(data_dir+"MaStar_spec_cuts/"):
        os.mkdir(data_dir+"MaStar_spec_cuts/")

if not os.path.exists(data_dir+"MaStar_specLum_cuts/"):
        os.mkdir(data_dir+"MaStar_specLum_cuts/")

progress = ProgressBar(cut_index.size, fmt=ProgressBar.FULL)
for ii in range(cut_index.size):
    progress.current += 1
    progress()
    this_star = cut_index[ii]
    MJD = mastar[1].data.field('MJD')[this_star]
    ra_string = '{:0>9.5f}'.format(mastar[1].data.field('OBJRA')[this_star])
    dec_string = '{:0=+10.5f}'.format(mastar[1].data.field('OBJDEC')[this_star])
    MANGAID_string = mastar[1].data.field('MANGAID')[this_star]
    mjd_string = '{:0>5}'.format(str(np.int(MJD)))
    filename_string = (ra_string + dec_string + "_"
                       + mjd_string + "_" + MANGAID_string+".txt")
    cut_filenames.append(filename_string)
    copied1 = check_output(["cp",main_dir+"MaStar_spec/"
                           +filename_string,data_dir
                           +"MaStar_spec_cuts/"])
    copied2 = check_output(["cp",main_dir+"MaStar_specLum/"
                           +filename_string,data_dir
                           +"MaStar_specLum_cuts/"])

cut_filenames_array = np.array(cut_filenames)

np.savetxt(data_dir+"list_of_all_MaStar_specCUT.txt", cut_filenames_array,
           delimiter=",", fmt="%s")
progress.done()

# In [61]: np.unique(mastar[1].data.field('SpecType')[cut_index])
# Out[61]:
# chararray(['A2', 'A3', 'A6', 'A7', 'A9', 'B6', 'F1', 'F2', 'F3', 'F4',
#            'F5', 'F6', 'F7', 'F8', 'F9', 'G0', 'G1', 'G2', 'G3', 'G4',
#            'G5', 'G6', 'G7', 'G8', 'G9', 'K0', 'K1', 'K2', 'K3', 'K4',
#            'K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
#            'WD4', 'WD5', 'WD6'], dtype='<U3')

# In [62]: np.unique(mastar[1].data.field('SpecType')[cut_index]).shape
# Out[62]: (43,)

cols = mastar[1].columns
data = mastar[1].data
new_hdu = fits.BinTableHDU(data[cut_index], name="ALLFILE")

col1 = fits.Column(name='Filename', format='41A',
                   array=cut_filenames_array)

coldefs = fits.ColDefs([col1])
hdu = fits.BinTableHDU.from_columns(coldefs)

new_columns = new_hdu.columns + hdu.columns
final_hdu = fits.BinTableHDU.from_columns(new_columns, name="ALLFILE")

new_HDUList = fits.HDUList(hdus=[mastar[0], final_hdu])
new_HDUList.writeto(data_dir+"mastarall-gaia-v2_4_3_specTypes_PS1_"
                    +"SDSSDR12_ALLFILEcutsFilename.fits", overwrite=True)
