from pyhamimports import *
from spectrum import Spectrum
from eyecheck import Eyecheck
from gui_utils import *
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import seaborn as sns

spec = Spectrum()

specPath = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/HARD_COPY_ORGINAL_DATA/SDSS_spec/ALL_VARSTAR_SPEC/ASCII/"
spec_fileName = "000.01306-002.43007_7850-56956-0321.txt"

ftype = None
message, ftype = spec.readFile(specPath+spec_fileName, ftype)

snVal = spec.calcSN()
spec.normalizeFlux()
measuredLines = spec.measureLines()
spec._lines = spec.measureLines()

lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
linesLabels = np.array(list(spec._lines.keys()))[np.argsort(list(spec._lines.keys()))]


CstarTemp_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/C_star_temps/"
WDTemp_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/WD_temps/"

Cstar_temp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/list_of_C_star_temps_11-17-2018.txt",dtype="U100")
WD_temp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/list_of_meaned_temps.txt",dtype="U100")
for ll in range(Cstar_temp_list.size):
    message, ftype = spec.readFile(CstarTemp_dir+Cstar_temp_list[ll], ftype)
    #snVal = spec.calcSN()
    #spec.normalizeFlux()
    measuredLines = spec.measureLines()
    spec._lines = spec.measureLines()

for kk in range(WD_temp_list.size):
    #spec = Spectrum()
    message, ftype = spec.readFile(WDTemp_dir+WD_temp_list[kk], ftype)
    #snVal = spec.calcSN()
    #spec.normalizeFlux()
    measuredLines = spec.measureLines()
    spec._lines = spec.measureLines()
    plt.plot(spec.wavelength,spec.flux)
    plt.show()

for kk in range(WD_temp_list.size):
    #spec = Spectrum()
    spec = fits.open(WDTemp_dir+WD_temp_list[kk])
    plt.plot(10.0**spec[1].data.field('loglam'),spec[1].data.field('flux'))
    plt.show()
#**************************************************************************
#**************************************************************************
pklPath = os.path.join(spec.thisDir, 'resources', 'tempLines.pickle')
with open(pklPath, 'rb') as pklFile:
    tempLines = pickle.load(pklFile)
    # Read in indices measured from templates
    # tempLines is a list of arrays with format: [spts, subs, fehs, lums, lines]
    # lines is a list of 2D arrays with indices and variances for each line
    # index for each spectrum that goes into a template
    # spts = 0,1,2,3,4,5,6,7 == O, B, A, F, G, K, M, L

tempLines_specTypesLabels = np.array(['O','B','A','F','G','K','M','L',])
tempLines_specTypes = np.array(tempLines[0])
tempLines_specSubTypes = np.array(tempLines[1])
tempLines_met = np.array(tempLines[2])

tempLines_arr = np.array(tempLines[4])

orginal_tempSpecLoc = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/PyHammer/PyHammer-master/resources/templates/"
CstarTemp_indivualSpec_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/spec/"
WDTemp_indivualSpec_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/WD_magCut_spec/"

otherTemp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/PyHammer/PyHammer-master/resources/allTempList.txt",dtype="U100")

new_tempLines_0 = tempLines[0]
new_tempLines_1 = tempLines[1]
new_tempLines_2 = tempLines[2]
new_tempLines_3 = tempLines[3]
new_tempLines_4 = []

ftype = None
for ii in range(tempLines_specTypes.size):
    line_name = ""
    line_name = line_name + tempLines_specTypesLabels[new_tempLines_0][ii]
    line_name = line_name + str(new_tempLines_1[ii])



new_tempLines_0 = np.append(np.append(tempLines[0], np.ones(Cstar_temp_list.size)*8), np.ones(WD_temp_list.size)*9)
new_tempLines_1 = np.append(np.append(tempLines[1], np.arange(Cstar_temp_list.size)+1), np.arange(WD_temp_list.size)+1)
new_tempLines_2 = np.append(np.append(tempLines[2], np.zeros(Cstar_temp_list.size)), np.zeros(WD_temp_list.size))
new_tempLines_3 = np.append(np.append(tempLines[3], np.ones(Cstar_temp_list.size)*5), np.ones(WD_temp_list.size)*5)
new_tempLines_4 = tempLines[4]

ftype = None
for ll in range(Cstar_temp_list.size):
    temp_list = []
    spec_used_in_this_temp = np.genfromtxt(CstarTemp_dir+Cstar_temp_list[ll][:-5]+"_specUsedInTemp.txt", dtype='U100')
    for jj in range(spec_used_in_this_temp.size):
        message, ftype = spec.readFile(CstarTemp_indivualSpec_dir+spec_used_in_this_temp[jj], ftype)
        measuredLines = spec.measureLines()
        spec._lines = spec.measureLines()
        lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
        linesLabels = np.array(list(spec._lines.keys()))[np.argsort(list(spec._lines.keys()))]
        temp_list.append(lines)
    new_tempLines_4.append(temp_list)

ftype = None
for ll in range(WD_temp_list.size):
    temp_list = []
    spec_used_in_this_temp = np.genfromtxt(WDTemp_dir+WD_temp_list[ll][:-5]+"_specUsedInTemp.txt", dtype='U100')
    for jj in range(spec_used_in_this_temp.size):
        message, ftype = spec.readFile(WDTemp_indivualSpec_dir+spec_used_in_this_temp[jj], ftype)
        measuredLines = spec.measureLines()
        spec._lines = spec.measureLines()
        lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
        linesLabels = np.array(list(spec._lines.keys()))[np.argsort(list(spec._lines.keys()))]
        temp_list.append(lines)
    new_tempLines_4.append(temp_list)

new_tempLines = [new_tempLines_0, new_tempLines_1, new_tempLines_2, new_tempLines_3, new_tempLines_4]

pklPath = os.path.join(spec.thisDir, 'resources', 'tempLines_01-17-2018_withClines.pickle')
with open(pklPath, 'wb') as pklFile:
    pickle.dump(new_tempLines, pklFile)
#**************************************************************************
#**************************************************************************

for ll in range(Cstar_temp_list.size):
    spec_used_in_this_temp = np.genfromtxt(CstarTemp_dir+Cstar_temp_list[ll][:-5]+"_specUsedInTemp.txt", dtype='U100')
    for jj in range(spec_used_in_this_temp.size):
        message, ftype = spec.readFile(CstarTemp_indivualSpec_dir+spec_used_in_this_temp[jj], ftype)
        spec._lines = spec.measureLines()
        spec.guessSpecType()
        print(jj, spec.guess)

for ll in range(WD_temp_list.size):
    spec_used_in_this_temp = np.genfromtxt(WDTemp_dir+WD_temp_list[ll][:-5]+"_specUsedInTemp.txt", dtype='U100')
    for jj in range(spec_used_in_this_temp.size):
        message, ftype = spec.readFile(WDTemp_indivualSpec_dir+spec_used_in_this_temp[jj], ftype)
        spec._lines = spec.measureLines()
        spec.guessSpecType()
        print(jj, spec.guess)

all_C_star_spec = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/all_C_star_spec.txt", dtype="U100")
guess_type = np.empty(all_C_star_spec.size, dtype=int)
for ll in range(all_C_star_spec.size):
    message, ftype = spec.readFile(CstarTemp_indivualSpec_dir+all_C_star_spec[ll], ftype)
    spec._lines = spec.measureLines()
    spec.guessSpecType()
    guess_type[ll] = spec.guess['specType']
    #print(ll, spec.guess)

C_percent_correct = np.where(guess_type == 8)[0].shape[0] / np.where(guess_type >= 0)[0].shape[0] #0.8909214092140921

all_WD_star_spec = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/all_WD_star_spec.txt", dtype="U100")
guess_type = np.empty(all_WD_star_spec.size, dtype=int)
WDTemp_indivualSpecALL_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/WD_ALL_spec/"
for ll in range(all_WD_star_spec.size):
    message, ftype = spec.readFile(WDTemp_indivualSpecALL_dir+all_WD_star_spec[ll], ftype)
    spec._lines = spec.measureLines()
    spec.guessSpecType()
    guess_type[ll] = spec.guess['specType']
    #print(ll, spec.guess)

WD_percent_correct = np.where(guess_type == 9)[0].shape[0] / np.where(guess_type >= 0)[0].shape[0]

for ii in range(10):
    print(np.where(guess_type == ii)[0].shape[0] / np.where(guess_type >= 0)[0].shape[0])

#**************************************************************************
#**************************************************************************
plt_dir = "/Users/benjaminroulston/Desktop/"
num_of_lines = 34
num_of_specTypes = tempLines_specTypesLabels.size
combined_line_data = np.ones((num_of_lines, num_of_specTypes, 2000), dtype=np.float64) * np.nan
for ii in range(num_of_lines):
    for jj in range(num_of_specTypes):
        this_specType_index = np.where(tempLines_specTypes == jj)[0]
        num_of_subTypes = this_specType_index.size
        for kk in range(num_of_subTypes):
            this_subType_data = np.array(tempLines_arr[this_specType_index][kk])
            this_subType_num_of_spec = this_subType_data.shape[0]
            combined_line_data[ii,jj,:this_subType_num_of_spec] = this_subType_data[:,ii,0]

# fig = plt.figure(figsize=(12,12), constrained_layout=True)
# gs = GridSpec(7, 5, figure=fig, hspace=0.5, wspace=0.5)#, height_ratios=[1, 1, 1, 1, 1, 1, 1], width_ratios=[1, 1, 1, 1, 1])
row = 0
column = 0
line_xlims = [[-0.2, 1.3], [-0.2, 1.2], [-1.25, 1.5], [0.875, 1.05], [0.7, 1.3], [0.6, 1.6], [-2,2], [-2, 2], [0.7, 1.8], [0.6, 1.75], [0.4, 1.5],
              [-1,2], [-1,2], [0.6, 1.6], [-1,2], [0,2], [-1,2], [-1,2], [0,2], [0,2], [0.5, 1.1], [0.75, 1.7], [0,1.5], [0.2,2.5],
               [0.3,1.3], [0.8, 1.8], [0.5,1.25], [0.9,1.3], [0.2,1.2],[-1, 11.3], [-0.5, 8],[-0.5, 3.5], [0,7], [-2,2]]
for ii in range(num_of_lines):
    if column > 4:
        column = 0
        row += 1
    #ax = fig.add_subplot(gs[row, column])
    for jj in range(8):
        data = combined_line_data[ii,jj,:][~np.isnan(combined_line_data[ii,jj,:])]
        #sns.kdeplot(data, label=tempLines_specTypesLabels[jj], ax=ax)
        sns.kdeplot(data, label=tempLines_specTypesLabels[jj])
    for kk in range(WD_temp_list.size):
        message, ftype = spec.readFile(WDTemp_dir+WD_temp_list[kk], ftype)
        snVal = spec.calcSN()
        spec.normalizeFlux()
        measuredLines = spec.measureLines()
        spec._lines = spec.measureLines()
        plt.axvline(x=spec.lines[linesLabels[ii]][0], color='blue', ls='dashed', alpha=0.5)
    for ll in range(Cstar_temp_list.size):
        message, ftype = spec.readFile(CstarTemp_dir+Cstar_temp_list[ll], ftype)
        snVal = spec.calcSN()
        spec.normalizeFlux()
        measuredLines = spec.measureLines()
        spec._lines = spec.measureLines()
        plt.axvline(x=spec.lines[linesLabels[ii]][0], color='red', ls='dashed', alpha=0.5)
    #ax.legend(loc='best')
    #ax.set_xlabel('Index')
    #ax.set_ylabel('PDF')
    #ax.set_title(linesLabels[ii])
    plt.xlim(line_xlims[ii])
    plt.legend(loc='best')
    plt.xlabel('Index')
    plt.ylabel('PDF')
    plt.title(linesLabels[ii])
    column +=1
    plt.savefig(plt_dir+"/plots/"+linesLabels[ii]+".eps",dpi=600,bbox_inches='tight')
    plt.clf()
    plt.close()

# h,l=ax.get_legend_handles_labels()
# ax = fig.add_subplot(gs[row, column])
# ax.legend(h,l) 

# plt.savefig(plt_dir+"/plots/PyHammer_IndexPlots.eps",dpi=600,bbox_inches='tight')
# #plt.show()
# plt.clf()
# plt.close()



