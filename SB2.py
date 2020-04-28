import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import bisect


class SB2(object):
    """SB2 spectrum class

    """

    def __init__(self, file1, spectype1, file2, spectype2):
        self.spectype1 = spectype1
        self.spectype2 = spectype2
        self.file1 = file1
        self.file2 = file2

        self._makeCompositeSB2(self.file1, self.file2)

    def _splitSpecType(self, s):
        head = s.rstrip('0123456789')
        tail = s[len(head):]
        return head, tail

    def _interpOntoGrid(self, wavelength, flux):
        """
        Description:
            A method to put the spectrum flux and variance onto the same
            wavelength grid as the templates (5 km/s equally spaced bins)
        """
        # Interpolate flux and variance onto the wavelength grid
        waveStart = 3_550
        waveEnd = 10_500
        waveNum = 65_000
        # 65,000 wavelengths gives 5km/s resolution across this region
        # dv = 2.9979e5 * (np.diff(waveGrid) / waveGrid[1:])
        self.waveGrid = np.logspace(np.log10(waveStart), np.log10(waveEnd), num=waveNum)

        interpFlux = np.interp(self.waveGrid, wavelength,
                               flux, right=np.nan, left=np.nan)

        # cut the grids off at 3650 and 10200 like the templates
        startIndex = bisect.bisect_right(self.waveGrid, 3650)
        stopIndex = bisect.bisect_right(self.waveGrid, 10200)

        wavelength = self.waveGrid[startIndex:stopIndex]
        flux = interpFlux[startIndex:stopIndex]

        return wavelength, flux

    def _makeCompositeSB2(self, spectype1, spectype2):
        self.flux = np.zeros(61617)
        self.wavelength = np.zeros(61617)
        self.error = np.zeros(61617)

        spec1 = fits.open(self.file1)
        self.wavelengthComponent1 = spec1[1].data['lam']
        self.fluxComponent1 = spec1[1].data['Lum']
        self.errorComponent1 = spec1[1].data['LumErr']

        self.specComponent1 = np.stack(
            (self.wavelengthComponent1,
             self.fluxComponent1, self.errorComponent1), axis=1)

        spec2 = fits.open(self.file2)
        self.wavelengthComponent2 = spec2[1].data['lam']
        self.fluxComponent2 = spec2[1].data['Lum']
        self.errorComponent2 = spec2[1].data['LumErr']

        self.specComponent2 = np.stack(
            (self.wavelengthComponent2, self.fluxComponent2,
                self.errorComponent2), axis=1)

        self.wavelength, spec1_interpFlux = self._interpOntoGrid(
            self.wavelengthComponent1, self.fluxComponent1)

        self.wavelength, spec1_interpError = self._interpOntoGrid(
            self.wavelengthComponent1, self.errorComponent1)

        self.wavelength, spec2_interpFlux = self._interpOntoGrid(
            self.wavelengthComponent2, self.fluxComponent2)

        self.wavelength, spec2_interpError = self._interpOntoGrid(
            self.wavelengthComponent2, self.errorComponent2)

        Lratio = spec1_interpFlux / spec2_interpFlux
        Lratio_20percent_20percent = np.where((Lratio >=
                                               0.8) & (Lratio <= 5.0))[
            0].shape[0] / Lratio.shape[0]

        self.LratioPercent = Lratio_20percent_20percent
        if Lratio_20percent_20percent >= 0.2:
            self.flux = spec1_interpFlux + spec2_interpFlux
            self.error = np.sqrt(spec1_interpError**2 + spec2_interpError**2)
            self.spec = np.stack(
                (self.wavelength, self.flux, self.error), axis=1)

            normFluxConst = self.flux[np.argmin(np.abs(self.wavelength -
                                                       8000.0))]
            self.normedFlux = self.flux / normFluxConst
            self.normederror = self.error / normFluxConst
            self.normedspec = np.stack(
                (self.wavelength, self.normedFlux, self.normederror), axis=1)
            self.isCOMBO = True
        else:
            self.isCOMBO = False
            # raise ValueError('This combo is invalid', self.spectype1,
            #                  self.spectype2)

    def plotCompositeSB2(
            self, saveplot=True,
            plotIndividual=True, filename=None, PyHammer=False):
        xmin = 3800
        xmax = 9000
        # sig_range = 3.0
        xmajor_tick_space = 1000
        xminor_tick_space = 100
        # ymajor_tick_space = 2
        # yminor_tick_space = 0.5
        flux_copy = self.flux.copy()
        if PyHammer:
            self.flux = self.normedFlux
        if plotIndividual:
            plt.plot(self.wavelengthComponent1, self.fluxComponent1,
                     label=self.spectype1, linewidth=0.50, color='#0392ff')
            plt.plot(self.wavelengthComponent2, self.fluxComponent2,
                     label=self.spectype2, linewidth=0.50,
                     color=(0.22745098039215686,
                            0.6901960784313725,
                            0.011764705882352941,
                            1.0))
            plt.plot(self.wavelength, self.flux, label='combined',
                     linewidth=0.5, color='#ff1c03')
            plt.legend(loc='best')
            plt.xlabel(r"Wavelength [$\rm{\AA}$]")  # , fontdict=font)
            plt.ylabel(r"Luminosity [W $\rm{\AA}_{-1}$]")
            ax = plt.gca()
            ax.set_xlim([xmin, xmax])
            # ax.set_yscale('log')
            # plt.set_ylim([ymin,ymax])
            ax.xaxis.set_major_locator(
                ticker.MultipleLocator(xmajor_tick_space))
            ax.xaxis.set_minor_locator(
                ticker.MultipleLocator(xminor_tick_space))

            if filename is None:
                filename = self.spectype1 + "+" + self.spectype2 + ".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
        else:
            plt.plot(self.wavelength, self.flux)
            plt.xlabel(r"Wavelength [$\AA$]")
            plt.ylabel("Luminosity [W]")
            if filename is None:
                filename = self.spectype1 + "+" + self.spectype2 + ".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
        plt.clf()
        plt.close()
        self.flux = flux_copy

    def saveCompositeSB2(self, filename=None, saveMethod='both'):
        if filename is None:
            filename = self.spectype1 + "+" + self.spectype2
        if saveMethod == 'PyHammer' or saveMethod == 'both':
            t = Table()
            t['LogLam'] = np.log10(self.wavelength)
            t['Flux'] = self.normedFlux
            t['PropErr'] = self.normederror
            t['Std'] = self.normederror
            t.write(filename + ".fits", format='fits', overwrite=True)
        if saveMethod == 'Lum' or saveMethod == 'both':
            t = Table()
            t['Lam'] = self.wavelength
            t['Lum'] = self.flux
            t['LumVar'] = self.error**2
            t['LumErr'] = self.error
            t.write(filename + ".fits", format='fits', overwrite=True)
