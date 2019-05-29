import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

import makeSB2temp as SB2

individualSpec_dir = "SB2_IndividualSpec/"
compositeSpec_dir = "SB2_CompositeSpec/"

spec_types = np.array(["A",
                       "F",
                       "G",
                       "K",
                       "M",
                       "WD"])



spec_list1 = "M_list.txt"
spec_list2 = "G_list.txt"

specFilenames1 = np.loadtxt(individualSpec_dir+spec_list1, dtype="U")
specFilenames2 = np.loadtxt(individualSpec_dir+spec_list2, dtype="U")

for ii, filename1 in enumerate(specFilenames1):
    for jj, filename2 in enumerate(specFilenames2):
        filename1 = filename1
        filename2 = filename2
        spectype1 = filename1.replace("/"," ").split()[1][:-4].replace("_"," ").split()[0].upper()
        spectype2 = filename2.replace("/"," ").split()[1][:-4].replace("_"," ").split()[0].upper()
        this_SB2 = SB2.SB2(individualSpec_dir+filename1, spectype1,
                           individualSpec_dir+filename2, spectype2)
        combinedFilename = (filename1.replace("/"," ").split()[1][:-4].replace("_"," ").split()[0].upper() +
                            "+"+
                            filename2.replace("/"," ").split()[1][:-4].replace("_"," ").split()[0].upper())

        this_SB2.plotCompositeSB2(filename=compositeSpec_dir+"plots/"+combinedFilename+".eps")
        this_SB2.saveCompositeSB2(filename=compositeSpec_dir+combinedFilename+".txt")