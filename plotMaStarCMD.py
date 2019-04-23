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

mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3_selfMatched.fits')
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

mastar_file_filenames = []
for ii, ID in enumerate(MANGAID):
    MJD = mastar[1].data.field('MJD')[ii]
    ra_string = '{:0>9.5f}'.format(RA[ii])
    dec_string = '{:0=+10.5f}'.format(DEC[ii])
    MANGAID_string = ID
    mjd_string = '{:0>5}'.format(str(np.int(MJD)))
    mastar_file_filenames.append(ra_string + dec_string + "_"
        + mjd_string + "_" + MANGAID_string)

mastar_file_filenames = np.array(mastar_file_filenames)

outputDIR = main_dir+"MaStar_specLum/"
all_filenames = np.genfromtxt(data_dir+"all_MaStar_spec.txt", dtype='str')

PyHammerResult = np.genfromtxt(data_dir+"PyHammerResults_MaStar.csv",
                               dtype='str', delimiter=",")
PyHammerResultCondensed = PyHammerResult[:, 0:4:3]

for ii in range(PyHammerResultCondensed[:, 0].size):
    PyHammerResultCondensed[ii, 0] = PyHammerResultCondensed[ii, 0].replace(
        "/", " ").split()[-1][:-4]

mathced_types_to_mastar = np.empty(M_G.size, dtype="U3")
for ii in range(mastar_file_filenames.size):
    try:
        index = np.where(PyHammerResultCondensed[:, 0]
                        == mastar_file_filenames[ii])[0][0]
        mathced_types_to_mastar[ii] = PyHammerResultCondensed[index, 1]
    except IndexError:
        print(mastar_file_filenames[ii])
        mathced_types_to_mastar[ii] = np.nan

spectypeNUM = np.empty(mathced_types_to_mastar.size, dtype=int)
# Spectral type, 0 for O to 7 for L, 8 = C, 9 = WD
specType_toNUM_NUM = np.arange(10)
specType_toNUM_alph = np.array(['O', 'B', 'A', 'F', 'G',
                                'K', 'M', 'L', 'C', 'WD'])
for ii in range(spectypeNUM.size):
    try:
        current_specType = mathced_types_to_mastar[ii]
        current_mainspecType, current_subType = splitSpecType(current_specType)
        current_subType = np.int(current_subType)
        spectypeNUM[ii] = (np.where(
                                    specType_toNUM_alph == current_mainspecType
                                    )[0][0]
                            * 10 + current_subType)
    except ValueError:
        spectypeNUM[ii] = -1

usable_spec_index = np.where(spectypeNUM != -1)[0]

cm = plt.cm.get_cmap('viridis')
sc = plt.scatter(BP_RP[usable_spec_index], M_G[usable_spec_index], s=5,
                 c=spectypeNUM[usable_spec_index], cmap=cm)
ax = plt.gca()
divider1 = make_axes_locatable(ax)
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(sc, cax=cax1)
cbar1.ax.get_yaxis().labelpad = 15
cbar1.ax.set_ylabel('SpecType', rotation=270)
cbar1.set_ticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
cbar1.set_ticklabels(['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'C', 'WD'])
# plt.scatter(np.log10(all_Per_ls[this_object_index]),
# np.log10(all_Amp_ls[this_object_index]),
# s=150.0, marker="X", color=single_point_color, edgecolors='red')
# ax = plt.gca()
ax.invert_yaxis()
ax.set_xlabel('BP - RP')
ax.set_ylabel('M$_{G}$')
plt.savefig("MaStar_CMD_spectypes.eps", dpi=600)
# plt.show()
plt.clf()
plt.close()

col1 = fits.Column(name='SpecType', format='3A',
                   array=mathced_types_to_mastar)
coldefs = fits.ColDefs([col1])
hdu = fits.BinTableHDU.from_columns(coldefs)

new_columns = mastar[1].columns + hdu.columns
new_hdu = fits.BinTableHDU.from_columns(new_columns)
new_hdu.writeto('mastarall-gaia-v2_4_3_selfMatched_specTypes.fits')
