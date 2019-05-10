import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from progressbar import ProgressBar

def splitSpecType(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return head, tail

main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
data_dir = "data/"

orginal_outputDIR = main_dir+"MaStar_spec/"
converted_outputDIR = main_dir+"MaStar_specLum/"

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3.fits')
mastarSpec = fits.open(main_dir+'mastar-goodspec-v2_4_3-v1_0_2.fits')

pyhammer_match = np.genfromtxt(data_dir+"PyHammerResults_05-03-2019.csv",
                               delimiter=",", skip_header=1, dtype="str")

ordered_specType = []

num_spec = mastar[2].data.size
progress = ProgressBar(num_spec, fmt=ProgressBar.FULL)
for ii in range(num_spec):
    progress.current += 1
    progress()
    MANGAID_string = mastar[2].data.field('MANGAID')[ii]
    MJD = mastar[2].data.field('MJD')[ii]
    ra_string = '{:0>9.5f}'.format(mastar[2].data.field('OBJRA')[ii])
    dec_string = '{:0=+10.5f}'.format(mastar[2].data.field('OBJDEC')[ii])
    mjd_string = '{:0>5}'.format(str(np.int(MJD)))
    filename_string = (ra_string + dec_string + "_"
                       + mjd_string + "_" + MANGAID_string)
    try:
        pyhammer_index = (np.where(pyhammer_match[:,0] ==
                                   orginal_outputDIR
                                   +filename_string
                                   +".txt")[0][0])
        ordered_specType.append(pyhammer_match[pyhammer_index, 3])
    except IndexError:
        ordered_specType.append("NO")
        
progress.done()

ordered_specType_array = np.array(ordered_specType)


spectypeNUM = np.empty(ordered_specType_array.size, dtype=int)
# Spectral type, 0 for O to 7 for L, 8 = C, 9 = WD
specType_toNUM_NUM = np.arange(10)
specType_toNUM_alph = np.array(['O', 'B', 'A', 'F', 'G',
                                'K', 'M', 'L', 'C', 'WD'])
for ii in range(spectypeNUM.size):
    try:
        current_specType = ordered_specType_array[ii]
        current_mainspecType, current_subType = splitSpecType(current_specType)
        current_subType = np.int(current_subType)
        spectypeNUM[ii] = (np.where(
                                    specType_toNUM_alph == current_mainspecType
                                    )[0][0]
                           * 10 + current_subType)
    except ValueError:
        spectypeNUM[ii] = -1


col1 = fits.Column(name='SpecType', format='3A',
                   array=ordered_specType_array)
col2 = fits.Column(name='SpecTypeNum', format='I',
                   array=spectypeNUM)
coldefs = fits.ColDefs([col1, col2])
hdu = fits.BinTableHDU.from_columns(coldefs)

new_columns = mastar[2].columns + hdu.columns
new_hdu = fits.BinTableHDU.from_columns(new_columns, name="GOODVISITS")
#new_hdu.writeto('mastarall-gaia-v2_4_3_selfMatched_specTypes.fits')

hdu1 = fits.BinTableHDU.from_columns(mastar[1].columns, name="GOODSTARS")

new_HDUList = fits.HDUList(hdus=[mastar[0], hdu1, new_hdu])
new_HDUList.writeto(data_dir+"mastarall-gaia-v2_4_3__specTypes.fits", overwrite=True)
