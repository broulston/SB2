from pyhamimports import *
from spectrum import Spectrum
from eyecheck import Eyecheck
from gui_utils import *
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import seaborn as sns

CstarTemp_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/C_star_temps/"
WDTemp_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/WD_temps/"
otherTemp_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/PyHammer/PyHammer-master/resources/templates/"

Cstar_temp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/C_temps/list_of_C_star_temps_11-17-2018.txt",dtype="U100")
WD_temp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/WD_PyHammer_temps/list_of_meaned_temps.txt",dtype="U100")
otherTemp_list = np.genfromtxt("/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/Spectral_fitting/PyHammer/PyHammer-master/resources/allTempList.txt",dtype="U100")

spec = Spectrum()

pklPath = os.path.join(spec.thisDir, 'resources', 'tempLines.pickle')
#pklPath = 'resources/tempLines.pickle'
#pklPath = 'resources/tempLines_12-04-2018.pickle'

with open(pklPath, 'rb') as pklFile:
    tempLines = pickle.load(pklFile)

tempLines_0 = tempLines[0]
tempLines_1 = tempLines[1]
tempLines_2 = tempLines[2]
tempLines_3 = tempLines[3]
tempLines_4 = tempLines[4]
# 0, 1, 2, 3, 4, 5, 6, 7 = O, B, A, F, G, K, M, L

specTypes = np.array(['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'C', 'WD'])

new_tempLines_0 = np.empty(otherTemp_list.size, dtype=int)
new_tempLines_1 = np.empty(otherTemp_list.size, dtype=int)
new_tempLines_2 = np.empty(otherTemp_list.size, dtype=np.float64)
new_tempLines_3 = np.ones(otherTemp_list.size, dtype=int) * 5
new_tempLines_4 = []

for ii in range(otherTemp_list.size):
    new_tempLines_0[ii] = np.where(specTypes == otherTemp_list[ii][0])[0][0]
    new_tempLines_1[ii] = otherTemp_list[ii][1]
    if len(otherTemp_list[ii].replace("_"," ").split()) == 1:
        new_tempLines_2[ii] = 0.
    else:
        new_tempLines_2[ii] = np.float64(otherTemp_list[ii].replace("_"," ").split()[1])

new_tempLines_0 = np.append(np.append(new_tempLines_0, np.ones(Cstar_temp_list.size)*8), np.ones(WD_temp_list.size)*9)
new_tempLines_1 = np.append(np.append(new_tempLines_1, np.arange(Cstar_temp_list.size)+1), np.arange(WD_temp_list.size)+1)
new_tempLines_2 = np.append(np.append(new_tempLines_2, np.zeros(Cstar_temp_list.size)), np.zeros(WD_temp_list.size))
new_tempLines_3 = np.append(np.append(new_tempLines_3, np.ones(Cstar_temp_list.size)*5), np.ones(WD_temp_list.size)*5)

all_spec_list = np.append(np.append(otherTemp_list, Cstar_temp_list), WD_temp_list)
ftype = None
for ii in range(all_spec_list.size):
    if all_spec_list[ii][0] == "C":
        message, ftype = spec.readFile(CstarTemp_dir+all_spec_list[ii], ftype)
        spec._lines = spec.measureLines()
        lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
        new_tempLines_4.append(lines)
    elif all_spec_list[ii][0] == "W":
        message, ftype = spec.readFile(WDTemp_dir+all_spec_list[ii], ftype)
        spec._lines = spec.measureLines()
        lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
        new_tempLines_4.append(lines)
    else:
        message, ftype = spec.readFile(otherTemp_dir+all_spec_list[ii], ftype)
        spec._lines = spec.measureLines()
        lines = np.array(list(spec._lines.values()))[np.argsort(list(spec._lines.keys()))]
        new_tempLines_4.append(lines)

new_tempLines = [new_tempLines_0, new_tempLines_1, new_tempLines_2, new_tempLines_3, new_tempLines_4]

pklPath = os.path.join(spec.thisDir, 'resources', 'tempLines_01-19-2018_C+WD_lines.pickle')
with open(pklPath, 'wb') as pklFile:
    pickle.dump(new_tempLines, pklFile)





