import makeSB2temp as SB2


def test_makeSB2():
    individualSpec_dir = "SB2_IndividualSpec/"
    spectype1 = "C2"
    spectype2 = "WD3"

    filename1 = spectype1 + ".fits"
    filename2 = spectype2 + ".fits"

    this_SB2 = SB2.SB2(individualSpec_dir + filename1, spectype1,
                       individualSpec_dir + filename2, spectype2)

    if this_SB2.isCOMBO:
        combinedFilename = spectype1 + spectype2
        this_SB2.plotCompositeSB2(filename=combinedFilename +
                                  ".eps")
        this_SB2.saveCompositeSB2(filename=combinedFilename)
