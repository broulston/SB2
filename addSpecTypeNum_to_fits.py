import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable


def splitSpecType(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return head, tail


main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
data_dir = "data/"

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3_SDSSDR12_specTypes.fits')
mastarSpec = fits.open(main_dir+'mastar-goodspec-v2_4_3-v1_0_2.fits')

MANGAID = mastar[1].data.field('MANGAID')
NVISITS = mastar[1].data.field('NVISITS')
RA = mastar[1].data.field('OBJRA')
DEC = mastar[1].data.field('OBJDEC')

PARALLAX = mastar[1].data.field('PARALLAX')
PARALLAX_ERROR = mastar[1].data.field('PARALLAX_ERROR')
PARALLAX_SNR = np.abs(PARALLAX / PARALLAX_ERROR)
R_EST = mastar[1].data.field('R_EST') * u.pc

BP_RP = mastar[1].data.field('BP_RP')
m_G = mastar[1].data.field('PHOT_G_MEAN_MAG')
M_G = m_G + 5.0 - 5.0*np.log10(R_EST.value)

guessed_SpecType_alpha = mastar[1].data.field('Guessed Spectral Type')


spectypeNUM = np.empty(guessed_SpecType_alpha.size, dtype=int)
# Spectral type, 0 for O to 7 for L, 8 = C, 9 = WD
specType_toNUM_NUM = np.arange(10)
specType_toNUM_alph = np.array(['O', 'B', 'A', 'F', 'G',
                                'K', 'M', 'L', 'C', 'WD'])
for ii in range(spectypeNUM.size):
    try:
        current_specType = guessed_SpecType_alpha[ii]
        current_mainspecType, current_subType = splitSpecType(current_specType)
        current_subType = np.int(current_subType)
        spectypeNUM[ii] = (np.where(
                                    specType_toNUM_alph == current_mainspecType
                                    )[0][0]
                           * 10 + current_subType)
    except ValueError:
        spectypeNUM[ii] = -1


col1 = fits.Column(name='SpecTypeNum', format='I',
                   array=spectypeNUM)
coldefs = fits.ColDefs([col1])
hdu = fits.BinTableHDU.from_columns(coldefs)

new_columns = mastar[1].columns + hdu.columns
new_hdu = fits.BinTableHDU.from_columns(new_columns, name="Unique Stars")
#new_hdu.writeto('mastarall-gaia-v2_4_3_selfMatched_specTypes.fits')

hdu2 = fits.BinTableHDU.from_columns(mastar[2].columns, name="All Vists")

new_HDUList = fits.HDUList(hdus=[mastar[0], new_hdu, hdu2])
new_HDUList.writeto(data_dir+"mastarall-gaia-v2_4_3_SDSSDR12_specTypesNum.fits", overwrite=True)
