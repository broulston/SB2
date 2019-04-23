import warnings
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import bisect

import makeSB2temp

SB2 = makeSB2temp.SB2("179.21537+018.78220_57502_3-112133641.txt", "G1","164.12292-002.58138_57434_3-110152864.txt", "WD1")
SB2.plotCompositeSB2() 
SB2.saveCompositeSB2()



#"206.48379+001.07330_57057_4-21387.txt", "M0"