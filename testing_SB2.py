import warnings
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import bisect

import makeSB2temp
#SB2 = makeSB2temp.SB2("A2","F4", 3.0)

SB2 = makeSB2temp.SB2("M1","WD2", 3.0)
SB2.plotCompositeSB2() 
SB2.saveCompositeSB2()