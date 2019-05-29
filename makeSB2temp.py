import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import bisect

main_dir = ("/Users/benjaminroulston/Dropbox/Research/TDSS/"
            + "Variable_Stars/WORKING_DIRECTORY/SB2_Composites/")
upper_data_dir = "data/"
normal_data_dir = "data/MaStar_specLum_cuts/"

WD_data_dir = "data/Levenhagen2017_data/Lspec/"

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
        self.waveGrid = 10**(5*0.43429448190325182 / 299792.458
                             * np.arange(0, 65000) + 3.55)

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

        spec1 = np.loadtxt(self.file1,
                           skiprows=1, delimiter=",")
        self.wavelengthComponent1 = spec1[:, 0]
        self.fluxComponent1 = spec1[:, 1]
        self.errorComponent1 = spec1[:, 2]

        self.specComponent1 = np.stack(
            (self.wavelengthComponent1,
             self.fluxComponent1, self.errorComponent1), axis=1)

        spec2 = np.loadtxt(self.file2,
                           skiprows=1, delimiter=",")
        self.wavelengthComponent2 = spec2[:, 0]
        self.fluxComponent2 = spec2[:, 1]
        self.errorComponent2 = spec2[:, 2]

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

        self.flux = spec1_interpFlux + spec2_interpFlux
        self.error = np.sqrt(spec1_interpError**2 + spec2_interpError**2)
        self.spec = np.stack((self.wavelength, self.flux, self.error), axis=1)

        normFluxConst = self.flux[np.where(
                                           abs(self.wavelength
                                               - 8000.0)
                                           == np.min(abs(self.wavelength 
                                                         - 8000.0))
                                           )[0]][0]
        self.normedFlux = self.flux / normFluxConst
        self.normederror = self.error / normFluxConst
        self.normedspec = np.stack((self.wavelength, self.normedFlux, self.normederror), axis=1)

    def plotCompositeSB2(
            self, saveplot=True,
            plotIndividual=True, filename=None, PyHammer=False):
        flux_copy = self.flux.copy()
        if PyHammer:
            self.flux = self.normedFlux
        if plotIndividual:
            plt.plot(self.wavelengthComponent1, self.fluxComponent1,
                     label=self.spectype1)
            plt.plot(self.wavelengthComponent2, self.fluxComponent2,
                     label=self.spectype2)
            plt.plot(self.wavelength, self.flux, label='combined')
            plt.legend(loc='best')
            plt.xlabel(r"Wavelength [$\AA$]")
            plt.ylabel("Luminosity [W]")
            plt.xlim([self.waveGrid.min(), self.waveGrid.max()])
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

    def saveCompositeSB2(self, filename=None, PyHammer=False):
        if filename is None:
            filename = self.spectype1+"+"+self.spectype2+".txt"

        header = "wavelength, Lum, Lum_error"
        if PyHammer:
            np.savetxt(filename, self.normedspec, delimiter=",",
                       header=header, fmt="%5.4f, %0.9e, %0.9e")
        else:
            np.savetxt(filename, self.spec, delimiter=",",
                       header=header, fmt="%5.4f, %0.9e, %0.9e")
