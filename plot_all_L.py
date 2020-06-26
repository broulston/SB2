import numpy as np
from astropy.io import fits
import astropy.units as u
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
                                  mask=(np.isnan(x) |
                                        np.isinf(x))).cumsum().filled(np.nan)
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

xmajor_tick_space = 1000
xminor_tick_space = 100

previous_spectype = ""
fig = plt.figure(figsize=(8, 9))
ax1 = plt.gca()
ax2 = ax1.twinx()
for ii, spectype_filename in enumerate(filenames):
    spectype = spectype_filename.replace("/", " ").split()[0]
    this_color = colors[np.where(specType_order == spectype)[0][0]]
    data = fits.open("SB2_IndividualSpec/" + spectype_filename)
    try:
        lam = data[1].data['lam']
    except KeyError:
        lam = 10.0**data[1].data['LogLam']

    Lum = data[1].data['Lum']

    smooth_lam = smoothFlux(lam)
    # smooth_Lum_w = smoothFlux(removeSdssStitchSpike(lam, Lum))
    # smooth_Lum_erg = (smooth_Lum_w * u.W / u.angstrom).to(u.erg / u.s / u.angstrom).value

    smooth_Lum_erg = smoothFlux(removeSdssStitchSpike(lam, Lum))
    smooth_Lum_w = (smooth_Lum_erg * u.erg / u.s / u.angstrom).to(u.W / u.angstrom).value
    first_of_spectype = (spectype != previous_spectype)
    if first_of_spectype:
        if spectype=='WD':
            ax1.plot(smooth_lam, smooth_Lum_erg,
                     color=this_color, label='DA',
                     linewidth=1.0, zorder=1)
            ax2.plot(smooth_lam, smooth_Lum_w,
                     color=this_color, label='DA',
                     linewidth=1.0, zorder=1)
            first_of_spectype = False
            previous_spectype = spectype
        elif spectype=='C':
            ax1.plot(smooth_lam, smooth_Lum_erg,
                     color=this_color, label='dC',
                     linewidth=1.0, zorder=1)
            ax2.plot(smooth_lam, smooth_Lum_w,
                     color=this_color, label='dC',
                     linewidth=1.0, zorder=1)
            first_of_spectype = False
            previous_spectype = spectype
        else:
            ax1.plot(smooth_lam, smooth_Lum_erg,
                     color=this_color, label=spectype,
                     linewidth=1.0, zorder=1)
            ax2.plot(smooth_lam, smooth_Lum_w,
                     color=this_color, label=spectype,
                     linewidth=1.0, zorder=1)
            first_of_spectype = False
            previous_spectype = spectype
    else:
        ax1.plot(smooth_lam, smooth_Lum_erg,
                 color=this_color,
                 linewidth=1.0, zorder=1)
        ax2.plot(smooth_lam, smooth_Lum_w,
                 color=this_color,
                 linewidth=1.0, zorder=1)
        previous_spectype = spectype

ax1.set_xlabel("Wavelength [\AA]")
ax1.set_ylabel("Luminosity [erg s$^{-1}$ \AA$^{-1}$]")
ax2.set_ylabel("Luminosity [W \AA$^{-1}$]")
# plt.legend(loc='best', ncol=3)
leg = plt.legend(loc='best', ncol=3, frameon=False)
# get the lines and texts inside legend box
leg_lines = leg.get_lines()
leg_texts = leg.get_texts()
# bulk-set the properties of all lines and texts
plt.setp(leg_lines, linewidth=1.5)
# plt.setp(leg_texts, fontsize='x-large')
ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.set_xlim([xmin, xmax])
ax1.set_ylim([(ymin * u.W).to(u.erg / u.s).value,
              (ymax * u.W).to(u.erg / u.s).value])

ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])

ax1.xaxis.set_major_locator(ticker.MultipleLocator(xmajor_tick_space))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(xminor_tick_space))

locmaj = ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
ax1.yaxis.set_major_locator(locmaj)

locmaj = ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
ax2.yaxis.set_major_locator(locmaj)

locmin = ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100)
ax1.yaxis.set_minor_locator(locmin)

locmin = ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100)
ax2.yaxis.set_minor_locator(locmin)

ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
ax1.yaxis.set_minor_formatter(ticker.NullFormatter())

plt.tight_layout()
plt.savefig("LumSpec.pdf", dpi=600)
plt.clf()
plt.close()
