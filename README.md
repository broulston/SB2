# SB2

[![GitHub release](https://img.shields.io/github/release-pre/broulston/SB2.svg)](https://github.com/broulston/SB2/releases/latest)
[![GitHub commits](https://img.shields.io/github/commits-since/broulston/SB2/v0.0.svg)](https://github.com/broulston/SB2/commits/master)
[![GitHub issues](https://img.shields.io/github/issues/broulston/SB2.svg)](https://github.com/broulston/SB2/issues)
[![license](https://img.shields.io/github/license/broulston/SB2.svg)](https://github.com/broulston/SB2/blob/master/license.txt)
[![Python Supported](https://img.shields.io/badge/Python%20Supported-3-brightgreen.svg)](conda)
[![Maintenance](https://img.shields.io/maintenance/yes/2020.svg)]()
[![Powered by Astropy](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org)

`SB2` is a Python package to create composite spectra for binary stars. 

`SB2` constructs spectra in natural units (i.e. Watts/&#x212B;) using GaiaDR2 distances. This is done using a library of luminosity normalized spectra. This library was created from a combination of the [MaStar](https://www.sdss.org/surveys/mastar/) survey from [SDSS-IV](https://www.sdss.org) and the [Pickles+1998](https://ui.adsabs.harvard.edu/abs/1998PASP..110..863P/abstract) library. The Pickles library was used for OBAF stars while MaStar and SDSS was used for the GKM, C, WD stars. 

`SB2` was created as a part of the v2.0.0 release of [PyHammer](https://github.com/BU-hammerTeam/PyHammer). More details of how the luminosity normalized spectra were created can be found in the corresponding paper [Roulston+2020]()

![Vi_plot](./LumSpec.png?raw=true)

The program simply takes input of the spectral types (with paths) and then computes if the luminosity values are within a specific range. If they are, the SB2 is created. The user can change any of these parameters to allow SB2s to always be made. 

The code can also take additional user templates provided/created templates as long as they are in the same units. 

---

## Example

```python
import makeSB2temp as SB2

this_SB2 = SB2("C2","WD3")
this_SB2.plotCompositeSB2() 
this_SB2.saveCompositeSB2()
```

---
![Vi_plot](./C2+WD3.png?raw=true)


