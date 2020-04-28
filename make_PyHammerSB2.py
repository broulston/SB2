# %load_ext autoreload
# %autoreload 2

import numpy as np
from SB2 import SB2


def run_SB2_combos(spec_list1, spec_list2,
                   individualSpec_dir, compositeSpec_dir):
    def split_filename(filename):
        filename = filename.replace("/", " ").split()[1][:-4]
        return filename.replace("_", " ").split()[0].upper()

    specFilenames1 = np.loadtxt(individualSpec_dir + spec_list1, dtype="U")
    specFilenames2 = np.loadtxt(individualSpec_dir + spec_list2, dtype="U")

    for ii, filename1 in enumerate(specFilenames1):
        for jj, filename2 in enumerate(specFilenames2):
            filename1 = filename1
            filename2 = filename2
            spectype1 = filename1.replace(
                "/", " ").split()[1][:-5].replace("_", " ").split()[0].upper()
            spectype2 = filename2.replace(
                "/", " ").split()[1][:-5].replace("_", " ").split()[0].upper()
            this_SB2 = SB2(individualSpec_dir + filename1, spectype1,
                           individualSpec_dir + filename2, spectype2)
            print(spectype1, spectype2, this_SB2.LratioPercent)
            if this_SB2.isCOMBO:
                combinedFilename = (split_filename(filename1) +
                                    "+" +
                                    split_filename(filename2))

                this_SB2.plotCompositeSB2(filename=compositeSpec_dir +
                                          "plots/" +
                                          combinedFilename +
                                          ".eps")
                this_SB2.saveCompositeSB2(filename=compositeSpec_dir +
                                          combinedFilename)


individualSpec_dir = "SB2_IndividualSpec/"
compositeSpec_dir = "SB2_CompositeSpec_04-27-2020/"

spec_types = np.array(["O", "B", "A", "F", "G", "K", "M", "C", "WD"])

# run_SB2_combos("A_list.txt", "F_list.txt",
#                individualSpec_dir, compositeSpec_dir)

for spectype1 in spec_types[:-1]:
    print(f"Making {spectype1} star SB2s:")
    secondary_types = spec_types[np.where(spec_types == spectype1)[0][0] + 1:]
    for spectype2 in secondary_types:
        run_SB2_combos(f"{spectype1}_list.txt", f"{spectype2}_list.txt",
                       individualSpec_dir, compositeSpec_dir)
