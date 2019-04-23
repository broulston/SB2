import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import bisect

main_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/SB2_Composites/"
data_dir ="data/"

class SB2(object):
    """SB2 spectrum class

    """
    def __init__(self, file1, spectype1, file2, spectype2):
        self.mastar = fits.open(data_dir
            + 'mastarall-gaia-v2_4_3_SDSSDR12_specTypes.fits')

        self.specTypes = self.mastar[1].data.field('Guessed Spectral Type') 
        self.unique_specTypes = np.unique(self.specTypes)

        if (np.where(self.unique_specTypes == spectype1)[0].size==1 
           & np.where(self.unique_specTypes == spectype2)[0].size==1):
            self.spectype1 = spectype1
            self.spectype2 = spectype2
            self.file1 = file1
            self.file2 = file2
        else:
            if (np.where(self.unique_specTypes == spectype1)[0].size==0 
               & np.where(self.unique_specTypes == spectype2)[0].size==0):
                raise ValueError("spectype1 AND spectype2 are not a valid choice!")

            if np.where(self.unique_specTypes == spectype1)[0].size==0:
                raise ValueError("spectype1 is not a valid choice!")

            if np.where(self.unique_specTypes == spectype2)[0].size==0:
                raise ValueError( "spectype2 is not a valid choice!")

        self.MANGAID = self.mastar[1].data.field('MANGAID')
        self.NVISITS = self.mastar[1].data.field('NVISITS')
        self.RA = self.mastar[1].data.field('OBJRA')
        self.DEC = self.mastar[1].data.field('OBJDEC')

        self.PARALLAX = self.mastar[1].data.field('PARALLAX')
        self.PARALLAX_ERROR = self.mastar[1].data.field('PARALLAX_ERROR')
        self.PARALLAX_SNR = np.abs(self.PARALLAX / self.PARALLAX_ERROR)
        self.R_EST = self.mastar[1].data.field('R_EST') * u.pc

        self.outputDIR = main_dir + "MaStar_specLum/"
        self.all_filenames = np.genfromtxt(data_dir + "all_MaStar_spec.txt", dtype='str')

        self._makeCompositeSB2(self.file1, self.file2)

    def _interpOntoGrid(self, wavelength, flux): 
        """
        Description:
            A method to put the spectrum flux and variance onto the same
            wavelength grid as the templates (5 km/s equally spaced bins)
        """
        # Interpolate flux and variance onto the wavelength grid
        self.waveGrid = 10**(5*0.43429448190325182 / 299792.458
                             * np.arange(0, 65000) + 3.55)

        interpFlux = np.interp(self.waveGrid, wavelength, flux, right=np.nan, left=np.nan)

        # cut the grids off at 3650 and 10200 like the templates
        startIndex = bisect.bisect_right(self.waveGrid, 3650)
        stopIndex = bisect.bisect_right(self.waveGrid, 10200)

        wavelength = self.waveGrid[startIndex:stopIndex]
        flux = interpFlux[startIndex:stopIndex]

        return wavelength, flux

    def plotCompositeSB2(
            self, saveplot=True,
            plotIndividual=True, filename=None):
        if plotIndividual:
            plt.plot(self.wavelengthComponent1, self.fluxComponent1,
                    label=self.spectype1)
            plt.plot(self.wavelengthComponent2, self.fluxComponent2,
                     label=self.spectype2)
            plt.plot(self.wavelength, self.flux, label='combined')
            plt.legend(loc='best')
            plt.xlabel("Wavelength [$\AA$]")
            plt.ylabel("Luminosity [W]")
            if filename is None:
                filename = self.spectype1 + "+" + self.spectype2 + ".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
            plt.clf()
            plt.close()
        else:          
            plt.plot(self.wavelength, self.flux)
            plt.xlabel("Wavelength [$\AA$]")
            plt.ylabel("Luminosity [W]")
            if filename is None:
                filename = self.spectype1 + "+" + self.spectype2 + ".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
            plt.clf()
            plt.close()

    def _makeCompositeSB2(self, spectype1, spectype2):
        self.flux = np.zeros(61617)
        self.wavelength = np.zeros(61617)
        self.error = np.zeros(61617)

        spec1 = np.loadtxt(data_dir + "MaStar_specLum_cuts/" + self.file1,
                           skiprows=1, delimiter=",")
        self.wavelengthComponent1 = spec1[:, 0]
        self.fluxComponent1 = spec1[:, 1]
        self.errorComponent1 = spec1[:, 2]

        self.specComponent1 = np.stack(
            (self.wavelengthComponent1,
             self.fluxComponent1, self.errorComponent1), axis=1)

        spec2 = np.loadtxt(data_dir + "MaStar_specLum_cuts/" + self.file2,
                           skiprows=1, delimiter=",")
        self.wavelengthComponent2 = spec2[:, 0]
        self.fluxComponent2 = spec2[:, 1]
        self.errorComponent2 = spec2[:, 2]

        self.specComponent2 = np.stack(
            (self.wavelengthComponent2, self.fluxComponent2, self.errorComponent2), axis=1)

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

    def saveCompositeSB2(self, filename=None):
        filename is None:
        filename = self.spectype1+"+"+self.spectype2+".txt"

        header = "wavelength, Lum[W], Lum_error [W]"
        np.savetxt(filename, self.spec, delimiter=",",
                   header=header, fmt="%5.4f, %0.9e, %0.9e")
