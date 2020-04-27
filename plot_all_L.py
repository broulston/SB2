import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def removeSdssStitchSpike(wavelength, flux):
    """
    All SDSS spectrum have a spike in the spectra
    between 5569 and 5588 angstroms where the two
    detectors meet. This method will remove the spike
    at that point by linearly interpolating across that gap.
    """
    # Make a copy so as to not alter the original, passed in flux
    flux = flux.copy()
    # Search for the indices of the bounding wavelengths on the spike. Use the
    # fact that the wavelength is an array in ascending order to search quickly
    # via the searchsorted method.
    lower = np.searchsorted(wavelength, 5569)
    upper = np.searchsorted(wavelength, 5588)
    # Define the flux in the stitch region to be
    # linearly interpolated values between
    # the lower and upper bounds of the region.
    flux[lower:upper] = np.interp(wavelength[lower:upper],
                                  [wavelength[lower], wavelength[upper]],
                                  [flux[lower], flux[upper]])
    return flux


def smoothFlux(flux):
    # Define our own nancumsum method since numpy's nancumsum was only
    # added recently and not everyone will have the latest numpy version
    def nancumsum(x):
        return np.ma.masked_array(x,
                                  mask=(np.isnan(x) | np.isinf(x))).cumsum().filled(np.nan)
    N = max(int(len(flux) / 600), 20)
    cumsum = nancumsum(np.insert(flux, 0, 0))
    smoothFlux = (cumsum[N:] - cumsum[:-N]) / N
    smoothFlux = np.append(flux[int(np.floor((N - 1) / 2))], smoothFlux)
    smoothFlux = np.append(smoothFlux, flux[-int(np.floor(N / 2)):])
    return smoothFlux


filenames = np.genfromtxt("SB2_IndividualSpec/ALLSPEC.txt", dtype="U")
specType_order = np.array(['O', 'B', 'A', 'F', 'G', 'K', 'M', 'C', 'WD'])
colors = ['black', 'purple', 'darkorchid', 'blue',
          'lightseagreen', 'green', 'darkorange', 'red', 'dimgray']

# n = specType_order.size
# colors = plt.cm.nipy_spectral(np.linspace(0,1,n))

xmin = 3800
xmax = 9000
ymin = 2e17
ymax = 5e27
sig_range = 3.0
xmajor_tick_space = 1000
xminor_tick_space = 100
ymajor_tick_space = 1e28
yminor_tick_space = 1e27

previous_spectype = ""
fig = plt.figure(figsize=(6, 9))
for ii, spectype_filename in enumerate(filenames):
    spectype = spectype_filename.replace("/", " ").split()[0]
    this_color = colors[np.where(specType_order == spectype)[0][0]]
    data = np.loadtxt("SB2_IndividualSpec/" + spectype_filename,
                      delimiter=",", skiprows=1)
    first_of_spectype = (spectype != previous_spectype)
    if first_of_spectype:
        plt.plot(smoothFlux(data[:, 0]), smoothFlux(removeSdssStitchSpike(
            data[:, 0], data[:, 1])), color=this_color,
            label=spectype, linewidth=1.0, zorder=1)
        first_of_spectype = False
        previous_spectype = spectype
    else:
        plt.plot(smoothFlux(data[:, 0]), smoothFlux(removeSdssStitchSpike(
            data[:, 0], data[:, 1])), color=this_color,
            linewidth=1.0, zorder=1)
        previous_spectype = spectype

plt.xlabel(r"Wavelength [$\rm{\AA}$]")
plt.ylabel(r"Luminosity [W $\rm{\AA}^{-1}}$]")
# plt.legend(loc='best', ncol=3)
leg = plt.legend(loc='best', ncol=3, frameon=False)
# get the lines and texts inside legend box
leg_lines = leg.get_lines()
leg_texts = leg.get_texts()
# bulk-set the properties of all lines and texts
plt.setp(leg_lines, linewidth=1.5)
# plt.setp(leg_texts, fontsize='x-large')
ax = plt.gca()
ax.set_yscale('log')
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.xaxis.set_major_locator(ticker.MultipleLocator(xmajor_tick_space))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(xminor_tick_space))

locmaj = ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
ax.yaxis.set_major_locator(locmaj)

locmin = ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(ticker.NullFormatter())

plt.tight_layout()
plt.savefig("SB2_IndividualSpec/LumSpec.eps", dpi=600)
plt.clf()
plt.close()
