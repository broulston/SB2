import warnings
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import bisect

main_dir = "/Users/benjaminroulston/Dropbox/Research/TDSS/Variable_Stars/WORKING_DIRECTORY/SB2_Composites/"
data_dir ="data/"

class SB2(object):
    """SB2 spectrum class

    """
    def __init__(self, spectype1, spectype2, snrCut=100.0):
        self.snrCut = snrCut
        self.spectype1 = spectype1
        self.spectype2 = spectype2

        self.mastar = fits.open(data_dir+'mastarall-gaia-v2_4_3.fits')
        self.mastarSpec = fits.open(main_dir+'mastar-goodspec-v2_4_3-v1_0_2.fits')

        self.MANGAID = self.mastar[1].data.field('MANGAID')
        self.NVISITS = self.mastar[1].data.field('NVISITS')
        self.RA = self.mastar[1].data.field('OBJRA')
        self.DEC = self.mastar[1].data.field('OBJDEC')

        self.PARALLAX = self.mastar[1].data.field('PARALLAX')
        self.PARALLAX_ERROR = self.mastar[1].data.field('PARALLAX_ERROR')
        self.PARALLAX_SNR = np.abs(self.PARALLAX / self.PARALLAX_ERROR)
        self.R_EST = self.mastar[1].data.field('R_EST') * u.pc

        self.outputDIR = main_dir+"MaStar_specLum/"
        self.all_filenames = np.genfromtxt(data_dir+"all_MaStar_spec.txt", dtype='str')

        self.PyHammerResult = np.genfromtxt(data_dir+"PyHammerResults_MaStar.csv", dtype='str', delimiter=",")
        self.PyHammerResultCondensed = self.PyHammerResult[:,0:4:3]

        for ii in range(self.PyHammerResultCondensed[:,0].size):
            self.PyHammerResultCondensed[ii,0] = self.PyHammerResultCondensed[ii,0].replace("/", " ").split()[-1][:-4]

        self.connectIndex = np.empty(self.all_filenames.size)
        for ii in range(self.connectIndex.size):
            try:
                self.connectIndex[ii] = np.where(self.PyHammerResultCondensed[:,0] == self.all_filenames[ii])[0][0]
            except IndexError:
                self.connectIndex[ii] = np.nan

        self.skipped_spec = np.where(~np.isfinite(self.connectIndex))[0]
        self.connectIndex = self.connectIndex.astype(int)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.highSNR = np.where(self.PARALLAX_SNR > self.snrCut)[0] 

        self.highSNR_specType = np.empty(self.highSNR.shape, dtype="U3")
        for ii in range(self.highSNR.size):
            self.highSNR_specType[ii] = self.PyHammerResultCondensed[self.connectIndex[self.highSNR[ii]],1]

        self._makeCompositeSB2(self.spectype1, self.spectype2)

    def _interpOntoGrid(self, wavelength, flux): 
        """
        Description:
            A method to put the spectrum flux and variance onto the same
            wavelength grid as the templates (5 km/s equally spaced bins)
        """
        # Interpolate flux and variance onto the wavelength grid
        self.waveGrid = 10**(5*0.43429448190325182/299792.458 * np.arange(0,65000) + 3.55)

        interpFlux = np.interp(self.waveGrid, wavelength, flux, right=np.nan, left=np.nan)

        #cut the grids off at 3650 and 10200 like the templates
        startIndex = bisect.bisect_right(self.waveGrid, 3650)
        stopIndex = bisect.bisect_right(self.waveGrid, 10200)

        wavelength = self.waveGrid[startIndex:stopIndex]
        flux = interpFlux[startIndex:stopIndex]

        return wavelength, flux

    def _getSpecTypesIndex(self, spectype1, spectype2):
        ErrorMessage1 = None
        ErrorMessage2 = None
        foundspecType_index1 = None
        foundspecType_index2 = None
        highSNR_index1 = np.where(self.highSNR_specType == spectype1)[0]
        if highSNR_index1.size >= 1:
            foundspecType_index1 = self.highSNR[highSNR_index1[0]]
        else:
            ErrorMessage1 = "No found matches to SpecType1"
        highSNR_index2 = np.where(self.highSNR_specType == spectype2)[0]
        if highSNR_index2.size >= 1:
            foundspecType_index2 = self.highSNR[highSNR_index2[0]]
        else:
            ErrorMessage2 = "No found matches to SpecType2"
        return foundspecType_index1, foundspecType_index2, ErrorMessage1, ErrorMessage2

    def plotCompositeSB2(self, saveplot=True, plotIndividual=True, filename=None):  
        if plotIndividual:
            plt.plot(self.wavelengthComponent1, self.fluxComponent1, label=self.spectype1)
            plt.plot(self.wavelengthComponent2, self.fluxComponent2, label=self.spectype2)
            plt.plot(self.wavelength, self.flux, label='combined')
            plt.legend(loc='best')
            plt.xlabel('Wavelength [$\AA$]')
            plt.ylabel('Luminosity [W]')
            if filename == None:
                filename = self.spectype1+"+"+self.spectype2+".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
            plt.clf()
            plt.close()
        else:          
            plt.plot(self.wavelength, self.flux)
            plt.xlabel('Wavelength [$\AA$]')
            plt.ylabel('Luminosity [W]')
            if filename == None:
                filename = self.spectype1+"+"+self.spectype2+".eps"
            if saveplot:
                plt.savefig(filename, dpi=600)
            else:
                plt.show()
            plt.clf()
            plt.close()

    def _makeCompositeSB2(self, spectype1, spectype2):
        self.foundspecType_index1, self.foundspecType_index2, self.ErrorMessage1, self.ErrorMessage2 = self._getSpecTypesIndex(self.spectype1, self.spectype2) 
        if self.foundspecType_index1==None or self.foundspecType_index2==None:
            print(self.ErrorMessage1)
            print(self.ErrorMessage2)
        else:
            self.flux = np.zeros(61617)
            self.wavelength = np.zeros(61617)
            for ii in [self.foundspecType_index1, self.foundspecType_index2]:
                MaStarIndex = ii# np.where(all_filenames == PyHammerFilenames[ii])[0]
                current_filename = self.all_filenames[MaStarIndex]
                current_spec = np.loadtxt(self.outputDIR+current_filename+".txt", skiprows=1, delimiter=",")

                if ii == self.foundspecType_index1:
                    self.wavelengthComponent1 = current_spec[:,0]
                    self.fluxComponent1 = current_spec[:,1]
                    self.errorComponent1 = current_spec[:,2]
                    self.error = self.errorComponent1
                    self.specComponent1 = np.stack((self.wavelengthComponent1, self.fluxComponent1, self.errorComponent1), axis=1)
                else:
                    self.wavelengthComponent2 = current_spec[:,0]
                    self.fluxComponent2 = current_spec[:,1]
                    self.errorComponent2 = current_spec[:,2]
                    self.specComponent2 = np.stack((self.wavelengthComponent2, self.fluxComponent2, self.errorComponent2), axis=1)

                self.wavelength, current_spec_interpFlux = self._interpOntoGrid(current_spec[:,0], current_spec[:,1])
                self.wavelength, current_spec_interpError = self._interpOntoGrid(current_spec[:,0], current_spec[:,2])
                if ii == self.foundspecType_index1:
                    self.error = current_spec_interpError

                self.flux += current_spec_interpFlux
                self.error = np.sqrt(self.error**2 + current_spec_interpError**2)
                self.spec = np.stack((self.wavelength, self.flux, self.error), axis=1)

    def saveCompositeSB2(self, filename=None):
            if filename == None:
                filename = self.spectype1+"+"+self.spectype2+".txt"

            header = "wavelength, Lum[W], Lum_error [W]"
            np.savetxt(filename, self.spec, delimiter=",", header=header, fmt="%5.4f, %0.9e, %0.9e") 










