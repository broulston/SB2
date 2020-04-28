from SB2 import SB2


def test_makeSB2():
    individualSpec_dir = "SB2_IndividualSpec/"
    spectype1 = "C2"
    spectype2 = "WD3"

    filename1 = "C/" + spectype1 + ".txt"
    filename2 = "WD/" + spectype2 + ".txt"

    this_SB2 = SB2(individualSpec_dir + filename1, spectype1,
                   individualSpec_dir + filename2, spectype2)

    if this_SB2.isCOMBO:
        combinedFilename = spectype1 + spectype2
        this_SB2.plotCompositeSB2(filename=combinedFilename +
                                  ".eps")
        this_SB2.saveCompositeSB2(filename=combinedFilename)
