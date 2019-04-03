# SB2

[![GitHub release](https://img.shields.io/github/release-pre/broulston/SB2.svg)](https://github.com/broulston/SB2/releases/latest)
[![GitHub commits](https://img.shields.io/github/commits-since/broulston/SB2/v0.0.svg)](https://github.com/broulston/SB2/commits/master)
[![GitHub issues](https://img.shields.io/github/issues/broulston/SB2.svg)](https://github.com/broulston/SB2/issues)
[![license](https://img.shields.io/github/license/broulston/SB2.svg)](https://github.com/broulston/SB2/blob/master/license.txt)
[![Python Supported](https://img.shields.io/badge/Python%20Supported-3-brightgreen.svg)](conda)
[![Maintenance](https://img.shields.io/maintenance/yes/2019.svg)]()
[![Powered by Astropy](https://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org)

`SB2` is a Python package to create composite spectra for binary stars. 

`SB2` constructs spectra in natural units (i.e. Watts) using GaiaDR2 distances. This is done using the [MaStar](https://www.sdss.org/surveys/mastar/) survey from [SDSS-IV](https://www.sdss.org). Currently this codes uses the `mastarall-gaia-v2_4_3.fits` and `mastar-goodspec-v2_4_3-v1_0_2.fits` data files which can be found at the [MaStar Data Access](https://www.sdss.org/dr15/mastar/mastar-data-access/) page.

Once the MaStar spectra are converted to natural units, they are simple added together based on user specified spectral type inputs (e.g. A4, WD1, M4, etc.)

---

## Example

```python
import makeSB2temp

SB2 = makeSB2temp.SB2("M1","WD2", snrCut=3.0)
SB2.plotCompositeSB2() 
SB2.saveCompositeSB2()
```

---



