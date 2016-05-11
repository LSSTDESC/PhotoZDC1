"""Classes and functions that manipulate the SEDs and the filters

  Classes:
  SED:                holds flux vs wavelength data
  Filter:             holds transmission vs wavelength data
  MaskSEDs            mask emission line regions of rest-frame SEDs
  ContinuumExtractor  extract continuum from SEDs
  EmissionLine: holds emission line model (TO BE IMPLEMENTED)
  
  Helper functions:
  createSedDict:        Reads list of SEDs from a file, then write the data into a dictionary (of SED objects)
  createFilterDict:     Reads list of filters from a file, then writes their data into a dictionary 
                        (of Filter objects)
  orderFiltersByLamEff: Order the filters in a dictionary by their effective wavelength, return ordered list
                        of strings of their names
  

"""
import os
import scipy.interpolate as interp
import scipy.integrate as integ
from sklearn.decomposition import PCA as sklPCA
import scipy.signal as signal
import numpy as np
import itertools
import collections

import sedMapper



class SED(object):
    """Spectral Energy Distribution
    
       ISSUE: change implementation so SED has given a wavelength grid that it is returned on. This
       is for treatment of emission lines which would require high resolution vs continuum of SED which 
       only requires low resolution. Then integration of the SED can just be a simple trapz sum integral
       as long as it can deal with un-even sampling.
    """

    def __init__(self, waveLengths, fluxes, igmModel=None, fluxUnits="wl"):
        """Initialise the spectral energy distribution

        @param waveLengths  in Angstroms 
        @param fluxes       in wavelength or frequency units 
        @param igmModel     model to use for IGM (or None)
        @param fluxUnits    units flux is in, either wavelength=wl or frequency=fq

        """        
 
        # keep track of original range and resolution of SED 
        self.lamMin = waveLengths[0]
        self.lamMax = waveLengths[-1]
        self.nLam = len(waveLengths)
        self.waveLengths = waveLengths
        
        if (self.lamMin<0.):
            raise ValueError("ERROR! wavelengths cannot be less than zero!")
            
        if (len(np.where(fluxes<0.)[0])>0):
            raise ValueError("ERROR! cannot have negative fluxes")
            
        # record the flux units 
        self.fluxUnits = fluxUnits
        
        # IGM model
        self.igmModel = igmModel

        # make interp object
        self.sed = interp.InterpolatedUnivariateSpline(waveLengths, fluxes, k=1)
        

    def __str__(self):
        """Return custom string representation of SED"""
        
        str_msg ='\n  SED Object: \n'
        str_msg += '  SED wavelength range = {0.lamMin!s} to {0.lamMax!s} angstroms\n'.format(self)
        str_msg += '  SED wavelength resolution = {0.nLam!s}\n'.format(self)
        str_msg += '  SED flux units = {0.fluxUnits!s}\n'.format(self)
        str_msg += '  IGM model = {0.igmModel!s}\n'.format(self)
        return str_msg
        

    def getFlux(self, lam, z=0.):
        """Return flux at wavelength lambda. If redshift = 0, lambda is the rest-frame
           If z>0 lambda is the observed frame 
           
           @param lam    wavelength in angstroms
           @param z      redshift of spectrum
        """
        if np.any(lam<0.):
            raise ValueError("ERROR! wavelength cannot be less than zero!")
        
        if (z<0.):
            raise ValueError("ERROR! redshift cannot be less than zero")

        lamFrame = lam/(1.+z) # does nothing if z=0 (i.e. this is the rest-frame)
        flux = self.sed(lamFrame)
        
        if np.any(flux<0.): #this is actually not needed
            flux[flux<0.] = 0.
            #print "Warning: negative flux was interpolated, set to zero instead"
            
        # add igm here
        if (self.igmModel != None):
            lamMeters = lam*1e-10
            flux *= self.igmModel.getObserverFrameTransmission(lamMeters, z)
        
        return flux


    def getSedData(self, **kwargs):
        """Return wavelength and flux data for wavelength grid supplied (e.g. ready for plotting)
           
           kwargs can be:
           @param lamMin       minimum wavelength in angstroms
           @param lamMax       maximum wavelength in angstroms
           @param nLam         number of points in wavelength grid
           @param wavelengths  wavelength grid
           @param z            redshift of spectrum
           
        """
        # Check if wavelengths are supplied as grid or min, max, num range
        if 'wavelengths' in kwargs:
            wavelengths = kwargs['wavelengths']
        elif ('lamMin' in kwargs and 'lamMax' in kwargs and 'nLam' in kwargs):
            wavelengths = np.linspace(kwargs['lamMin'], kwargs['lamMax'], kwargs['nLam'])
        else:
            raise ValueError("ERROR! wavelength grid incorrectly specified")
            
        # Check for redshift
        z = 0
        if 'z' in kwargs:
            z = kwargs['z']
            
        fluxes = self.getFlux(wavelengths, z)
        
        if 'wavelengths' in kwargs:
            return fluxes
        else:
            return wavelengths, fluxes
        

    def setEmLine(self, emLineModel):
        """Add emission lines to the spectrum """
 
        # remake interp object so now includes emission lines as described by emLineModel
        raise RuntimeError("ERROR! Emission line model not implemented yet")


    def returnSedRange(self):
        """Return range of rest-frame wavelengths interpolation is valid within"""
        return self.lamMin, self.lamMax


    def returnSedResolution(self):
        """Return number of wavelengths that had defined flux values between lamMin and lamMax"""
        return self.nLam


    def returnFluxUnits(self):
        return self.fluxUnits
        


class Filter(object):
    """ Filter transmission function
    
    
    """
    
    def __init__(self, waveLengths, transmission):
        """Initialise the filter transmission function

        @param waveLengths  
        @param transmission

        """        
 
        # keep track of original range and resolution of transmission function 
        self.lamMin = waveLengths[0]
        self.lamMax = waveLengths[-1]
        self.nLam = len(waveLengths)   
		#this could come handy if we need to defer the construction of the interpolator
        self.wavelengths = waveLengths
        self.transmission = transmission
        #cache these values for later reuse. Note that this will fail if the filter
        #support is disconnected... need to build a safer interface?
        self.integral1 = integ.trapz(transmission*waveLengths,waveLengths)
        self.integral2 = integ.trapz(transmission/waveLengths,waveLengths)
        self.effLambda = np.sqrt(self.integral1/self.integral2)
        
        # check wavelength and transmission domains are valid
        if (self.lamMin<0.):
            raise ValueError("ERROR! wavelengths cannot be less than zero!")
            
        if np.any(transmission<0.):
            raise ValueError("ERROR! cannot have negative transmission")
           
        if np.any(transmission>1.000001): #why not just 1.?
            raise ValueError("ERROR! cannot have transmission > 1") 

        # filter is now represeted as an interpolation object (linear)
        self.filt = interp.InterpolatedUnivariateSpline(waveLengths, transmission, k=1)
        
        
    def __str__(self):
        """Return custom string representation of Filter"""
        
        lamEff = self.getFilterEffectiveWL()
        str_msg ='\n  Filter Object: \n'
        str_msg += '  Filter wavelength range = {0.lamMin!s} to {0.lamMax!s} angstroms\n'.format(self)
        str_msg += '  Filter wavelength resolution = {0.nLam!s}\n'.format(self)
        str_msg += '  Filter effective wavelength = ' + str(int(round(lamEff))) + ' angstroms\n'
        return str_msg
        
        
    def getTrans(self, lam):
        """Return transmission at wavelength lambda."""
        
        if np.any(lam<0.):
            raise ValueError("ERROR! wavelength cannot be less than zero!")

        trans = self.filt(lam)
        
        # protect against interpolation producing negative values
        trans[trans<0] = 0
        return trans
        
    
    def getFilterData(self):
        """Return wavelength and transmission data for filter (e.g. ready for plotting)

        """
        wavelengths = np.linspace(self.lamMin, self.lamMax, self.nLam)
        trans = self.getTrans(wavelengths)
        return wavelengths, trans
        
    
    #def getFilterEffectiveWL_old(self):
    #    """Calculate effective wavelength of the filter"""
    #    # should this be just int F(lam)*lam dlam / int F(lam) dlam ????
    #    
    #    top = integ.quad(self._integrand1, self.lamMin, self.lamMax)[0]
    #    bot = integ.quad(self._integrand2, self.lamMin, self.lamMax)[0]
    #    lamEff = np.sqrt(top/bot)
    #    return lamEff

	#this version of the function runs ~300 times faster on my laptop.
    def getFilterEffectiveWL(self):
        return self.effLambda
        
    def returnFilterRange(self):
        """Return range of filter wavelengths"""
        return self.lamMin, self.lamMax


    def returnFilterResolution(self):
        """Return number of wavelengths that had defined flux values between lamMin and lamMax"""
        return self.nLam
        
        
    def _integrand1(self, lam):
        """Return filter transmission * lambda """
        return self.getTrans(lam)*lam
    
    
    def _integrand2(self, lam):
        """Return filter transmission / lambda """
        return self.getTrans(lam)/lam
        
        
    
class MaskSEDs(object):


    def __init__(self, origSEDs, emission_lines_file):
        """For SEDs in origSEDs mask out wavelength ranges read from emission_lines_file. Wavelength ranges
           correspond to regions where strong nebular emission lines occur in galaxies.
        
           
           @param origSEDs               dictionary containing SED set
           @param emission_lines_file    file containing regions of SED to mask
        """
        
        self.origSEDs = origSEDs
        self.nsed = len(origSEDs)
        self.wlnorm = 5500.  # normalise SEDs to 1 at 5500A
        
        # new dictionary of the masked SEDs
        self.maskedSEDs = {}
        
        # wavelength ranges containing common emission lines in angstroms
        self.el = np.loadtxt(emission_lines_file)
        
        
        
      
    def mask_SEDs(self):
        """Do SED masking """
        
        for ised, (sedname, sed) in enumerate(self.origSEDs.items()):
        #for ised (sedname, sed) in enumerate(self.origSEDs.items()):

            # constant to divide by so SED is normalised to 1 at wlnorm
            flux_norm = sed.getFlux(self.wlnorm)
            
            # get wavelength grid for this SED
            wlgrid = sed.waveLengths
            nwl = len(wlgrid)
            
            # array of masked SED fluxes
            masked_sed = np.zeros((nwl, ))
            #print masked_sed.shape
            #print sedname

            # book-keeping for which emission line masking range we are on
            iLine = 0 
    
            iterator = zip( wlgrid, xrange(nwl) ).__iter__()
            for wl, iwl in iterator:

                # if wl within masking region, perform masking
                if (iLine<len(self.el) and wl>=self.el[iLine,0]):
                
                    ### below works if resolution can be assumed to be constant
                    # get resolution at this wavelength
                    #dwl = wlgrid[iwl+1] - wlgrid[iwl]
                    
                    # number of wavelength grid points within masking range
                    #
                    #nwlr = int( round( (self.el[iLine,1]-self.el[iLine,0])/dwl ) )
            
                    # end wl grid point within masking range
                    #wln = wl + dwl*nwlr
                    
                    
                    # To allow for varying resolution 
                    # get end wl grid point within masking range
                    # get number of wavelength grid points within masking range
                    nwlr = iwl
                    wln = wlgrid[nwlr]
                    wlregion = [wl]
                    while (wln<self.el[iLine,1]):
                        nwlr += 1
                        wln = wlgrid[nwlr]
                        wlregion.append(wln)
                    nwlr -= iwl
            
                    # interpolate across this space
                    y0 = sed.getFlux(wl)/flux_norm
                    yn = sed.getFlux(wln)/flux_norm        
                    #print masked_sed[iwl:iwl+nwlr+1].shape
                    #print self.el[iLine,0],self.el[iLine,1]
                    #print wlregion
                    #print self._mini_interp(wlregion, y0, yn)
                    #print nwlr, iwl, iwl+nwlr+1
                    masked_sed[iwl:iwl+nwlr+1] = self._mini_interp(wlregion, y0, yn)
                    #self._mini_interp(wl, wln, y0, yn, dwl)
            
                    # advance iterator to skip to the end of this masking region
                    iLine += 1
                    self._consume(iterator, nwlr)
            
                else:
                    masked_sed[iwl] = sed.getFlux(wl)/flux_norm
                    
            sed = SED(wlgrid, masked_sed)
            self.maskedSEDs[sedname] = sed

    
            # temporary time saver while testing
            #if (ised>1):
            #    break
                 
      
    def return_masked_SEDs(self):
        """Return newly masked SEDs """        
        return self.maskedSEDs
        
        
    def write_masked_SEDs(self, pathToFile, output_file_ext):
        """ """
        for ised, (sedname, sed) in enumerate(self.maskedSEDs.items()):
    
            fname = pathToFile + sedname + output_file_ext
            
            wlgrid = sed.waveLengths
            
            flux = []
            for wl in wlgrid:
                flux.append(sed.getFlux(wl))
            flux = np.asarray(flux)
    
            X = np.column_stack((wlgrid, flux))
            np.savetxt(fname, X)
            
      
    def set_wavelength_norm(self, wlnorm):
        """Reset wavelength to normalise SEDs to 1 at
        
           @param wlnorm   wavelength in Angstroms
        """
        self.wlnorm = wlnorm
        
        
    
    def _consume(self, iterator, n):
        """eats the next n iterations from iterator"""
        collections.deque(itertools.islice(iterator, n))
    
    
    def _mini_interp(self, xregion, y0, yn):
        """interpolate between x0 and xn, include points at both x0 and xn on return
    # 
    #        x0,y0   start point of interpolation
    #       xn,yn   end point of interpolation
    #       dx      exact spacing between x0 and xn 
    #    """
    
        np = len(xregion)
        x0 = xregion[0]
        xn = xregion[-1]

        # linear interp parameters
        m = (yn-y0)/(xn-x0)
        c = y0 - m*x0
    
        yinterp = []
        for x in xregion:
            yinterp.append(m*x + c)
     
        return yinterp

    
    #def _mini_interp(self, x0, xn, y0, yn, dx):
    #    """interpolate between x0 and xn, include points at both x0 and xn on return
    # 
    #        x0,y0   start point of interpolation
    #       xn,yn   end point of interpolation
    #       dx      exact spacing between x0 and xn 
    #    """
    #    np = (xn-x0)/dx
    #    if (np%2 != 0):
    #        raise ValueError("ERROR! x0 and xn are not evenly spaced")
    #
    #    # linear interp parameters
    #    m = (yn-y0)/(xn-x0)
    #    c = y0 - m*x0
    # 
    #    np = int(np+1)
    #    yinterp = []
    #    for i in range(np):
    #        x = x0 + i*dx
    #        yinterp.append(m*x + c)
    #    
    #    return yinterp
    
class ContinuumExtractor(object):
    """Given SEDs (continuum+emission lines) return continuum part of SEDs only
    
       Prescription:

       1) crude masking by linearly interpolating across emission line regions (MaskSEDs)
       2) PCA a library of synthetic spectra
       3) project crudely masked SEDs onto these PC's
       4) reconstruct the crudely masked SEDs from the PC's
       5) divide out smooth reconstructed SEDs
       6) high pass filter to remove high frequency noise etc
    """
    
    def __init__(self, seds, syntheticSEDs, mask_file, minWavelen=600., maxWavelen=11000., nWavelen=10000,
                       nComps=-1):
        """Initialise PCA 
          
           @param seds          SEDs to find continuum of
           @param sytheticSEDs  Library of synthetic SEDs (no emission lines)
           @param mask_file     file containing wavelength regions to mask
           @param minWavelen    minimum wavelength in Angstroms
           @param maxWavelen    maximum wavelength in Angstroms
           @param nWavelen      number of steps in wavelength grid 
           @param nComps        number of components to keep (if -1 keep all)
        """
        if (nComps<0):
            nComps = len(syntheticSEDs)
        
        self.minWavelen = minWavelen
        self.maxWavelen = maxWavelen
        self.nWavelen = nWavelen
        
        self.sedsContinuum = {}  
        
        ### Find PC's of synthetic spectra
        self._doPCA(nComps, syntheticSEDs)
        
        ### Mask out wavelength regions likely to have emission lines
        self._doMask(seds, mask_file)
        
        
    def getContinua(self, window):
        """Estimate continuum for each SED in library
        
           @param window    width of rolling median window to remove high frequency noise in Angstroms
        """
        
        ### Return list of SED names and SEDs on an array identical to that used to PCA the 
        ### synthetic spectra 
        name_list = self.sedsMasked.keys()
        waveLen, sed_array = sedMapper.get_sed_array(self.sedsMasked, 
                                                        self.minWavelen, self.maxWavelen, self.nWavelen)
                            
        ### Convert window in Angstroms to window in pixels (that is an odd integer)
        window_pix = int(window/(waveLen[1] - waveLen[0]))
        if (window_pix % 2 == 0):
            window_pix += 1
                            
        for i in range(len(self.sedsMasked)):
        
            ### Pick SED out of masked library
            sed = sed_array[i,:]
            
            ### Reconstruct using the eigenspectra of the synthetic library (no emission lines/noise)
            sed_rec = self._reconstructSmoothSpectrum(sed)
            
            ### Return continuum after high-pass filtering to get rid of high frequency noise 
            sed_continuum = self._getContinuum(sed, sed_rec, window_pix)
            
            ### Turn into SED object
            sedContinuum = SED(waveLen, sed_continuum)
            
            ### Add to dictionary
            self.sedsContinuum[name_list[i]] = sedContinuum
            
            
    def returnContinua(self):
        """Return dictionary containing SED continuua"""
        
        if (len(self.sedsContinuum)<1):
            raise ValueError("Error! Continua have not yet been estimated")
        return self.sedsContinuum
        
    
    def _getContinuum(self, sed, sed_rec, window):
        """Return continuum after high-pass filtering to get rid of high frequency noise 
           
           @param sed       original masked SED
           @param sed_rec   SED after reconstruction using eigenspectra of the synthetic library
           @param window    width of rolling median window to remove high frequency noise
        """
        
        ratio = sed/np.squeeze(sed_rec.T)
        smooth = signal.medfilt(ratio, window)
        
        ### Project against smooth being zero
        smooth[ np.where(smooth<1e-12) ] = np.min( smooth[ np.where(smooth>0) ] )
        sed = sed/smooth
        
        ### Protect against negative fluxes
        sed[np.where(sed<0)] = 0.
        
        return sed
    
    
    def _reconstructSmoothSpectrum(self, sed):
        """Reconstruct using the eigenspectra of the synthetic library (no emission lines/noise)
          
           @param sed    SED to reconstruct
        """
        coeffs = np.array(self.specPCA.transform(sed)) 
        sed_rec = np.dot(coeffs, self.specPCA.components_) + self.specPCA.mean_
        norm = np.sum(sed_rec)
        sed_rec /= norm
        return sed_rec
    
        
    def _doMask(self, seds, mask_file):
         """Mask out regions likely to have emission lines
         
            @param seds         SEDs to estimate continuum of
            @param mask_file    list of wavelength regions to mask
         """

         maskSeds = MaskSEDs(seds, mask_file)
         maskSeds.mask_SEDs()
         self.sedsMasked = maskSeds.return_masked_SEDs()
         
        
    def _doPCA(self, nComps, syntheticSEDs):
        """Do PCA on the synthetic spectra
        
           @param nComps           number of components to keep
           @param syntheticSEDs    library of synthetic spectra
        """
    
        waveLen, smooth_spectra = sedMapper.get_sed_array(syntheticSEDs, 
                                                          self.minWavelen, self.maxWavelen, self.nWavelen)
    
        self.specPCA = sklPCA(nComps)
        self.specPCA.fit(smooth_spectra)



class MakeBandpass(object):
    
    def __init__(self, wavelen_min=3000., wavelen_max=12000., wavelen_step=10.):
        """Initialise bandpass as transmission=1 along supplied wavelength grid
        
           @param wavelen_min    minimum wavelength in Angstroms
           @param wavelen_max    maximum wavelength in Angstroms
           @param wavelen_step   step in wavelength in Angstroms
        """
        # initialise wavelength grid
        self.setWavelenLimits(wavelen_min, wavelen_max, wavelen_step)
        self.wavelen = np.arange(self.wavelen_min, self.wavelen_max+self.wavelen_step/2., self.wavelen_step,
                                    dtype='float')
                                    
        # initialise transmission to 1.
        self.trans = np.ones(len(self.wavelen), dtype='float')
    

    def setWavelenLimits(self, wavelen_min, wavelen_max, wavelen_step):
        """
        Set internal records of wavelen limits, _min, _max, _step.
        """
        # If we've been given values for wavelen_min, _max, _step, set them here.
        if wavelen_min is not None:
            self.wavelen_min = wavelen_min
        if wavelen_max is not None:
            self.wavelen_max = wavelen_max
        if wavelen_step is not None:
            self.wavelen_step = wavelen_step
        return
        
        
    def make_band(self, componentList=['detector.dat', 'lens1.dat','lens2.dat', 'lens3.dat',
                                                'm1.dat', 'm2.dat', 'm3.dat','atmos.dat'],
                                                path_to_files = '.'):
                                                
        for component in componentList:
            
            transmission_data = np.loadtxt(path_to_files +"/" + component)
            if (np.max(transmission_data[:,1])>1.):
                print "ARGEHHH max transmission of", component , np.max(transmission_data[:,1])
            trans = interp.InterpolatedUnivariateSpline(transmission_data[:,0], transmission_data[:,1], k=1)
            self.trans *= trans(self.wavelen)
            
            # set to zero all negative transmissions
            self.trans[np.where(self.trans<0)] = 0.
            
            #if (np.max(self.trans)>1.):
            #    print "ARGEHHH"
            
        
    def write_throughput(self, filename, path_to_file):
        """Write throughput to filename """

        trans = np.transpose(np.vstack((self.wavelen, self.trans)))
        np.savetxt(path_to_file + "/" + filename, trans)
        
    
class EmissionLine(object):

    def __init__(self, emissionLinePars):
    
        self.emissionLinePars = emissionLinePars
        print "EmissionLine class not yet implemented"
        pass

        
#### Helper functions ####

def createSedDict(listOfSedsFile, pathToFile="../sed_data/", normMeth='', doPrinting=True):
    """Read file containing list of SEDs to read, then read these SEDs into a dictionary
    
       Dictionary keyword is a string: name of the file that contained the filter (without path or extension)
       Dictionary value is a Filter object: the filter transmission data read from file 
       
    """
    f = open(pathToFile + "/" + listOfSedsFile)
        
    sedDict = {}
    for line in f:
           
        sedData = np.loadtxt(pathToFile + "/" + line.rstrip())
        # sedName = line.rstrip().split('.')[0] -> only works if single '.' in filename delineating extension
        sedName = "_".join(line.rstrip().split('.')[:-1])
        
        if doPrinting:
            print "Adding SED", sedName ,"to dictionary"
            
        # check wavelength grid
        wavelengths = sedData[:,0]
        nlam = len(wavelengths)
        """
        dlam1 = abs(wavelengths[1]-wavelengths[0])
        dlam2 = abs(wavelengths[2]-wavelengths[1])
        
        if (dlam1-dlam2)<1e-6:
            isLinearSpaced = True
            isLogSpaced = False
            print "Wavelengths are linear spaced"
        elif (dlam1/(wavelengths[0]+dlam1/2.) - dlam2/(wavelengths[0]+dlam2/2.))<1e-6:
            isLinearSpaced = False
            isLogSpaced = True
            print "Wavelengths are log spaced"
        else:
            isLinearSpaced = False
            isLogSpaced = False
            print "Wavelengths are log spaced"
        """
            
        # method of normalising the SED
        if normMeth=='':
            # no normalisation done
            flux = sedData[:,1]
        elif normMeth=='sumTo1':
            # normalise so it sums to 1
            flux = sedData[:,1]/np.sum(sedData[:,1])
        elif normMeth=='intTo1':
            # normalise so it (roughly) integrates to 1
            intval = 0
            for i in range(nlam-1):
                dy = abs(sedData[i+1,1] - sedData[i,1])
                dx = wavelengths[i+1] - wavelengths[i]
                intval += dy*dx
            flux = sedData[:,1]/intval
        elif normMeth=='setAt5500':
            # normalise so it is equal to 1 at 5500A
            inorm = np.argmin(abs(wavelengths-5500.))
            flux = sedData[:,1]/sedData[inorm,1]
        elif normMeth=='setAt8000':
            # normalise so it is equal to 1 at 8000A
            inorm = np.argmin(abs(wavelengths-8000.))
            flux = sedData[:,1]/sedData[inorm,1]   
        else:
            raise ValueError("ERROR! unknown normalisation method")
            
        sed = SED(wavelengths, flux)
        sedDict[sedName] = sed
           
    return sedDict
    

def createFilterDict(listOfFiltersFile, pathToFile="../filter_data/"):
    """Read file containing list of filters to read, then read these filters into a dictionary
    
       Dictionary keyword is a string: name of the file that contained the filter (without path or extension)
       Dictionary value is a Filter object: the filter transmission data read from file 
       
    """
    
    f = open(pathToFile + "/" + listOfFiltersFile)
        
    filterDict = {}
    for line in f:
           
        filtData = np.loadtxt(pathToFile + "/" + line.rstrip())
        filtName = line.rstrip().split('.')[0]
        print "Adding filter", filtName ,"to dictionary"
        
        filt = Filter(filtData[:,0], filtData[:,1])
        filterDict[filtName] = filt
           
    return filterDict
    
    
def plotSedFilter(dataDict, ax, isSED=True, lamMin=2500., lamMax=12000., nLam=5000, color='red'):
    """Plots the data in dataDict on matplotlib axes given by ax"""
    
    lamNorm = 5500.
    
    for dataname, data in dataDict.items():
        print dataname,
    
        ### Get array version of data back
        if (isSED):
            wl, data_array = data.getSedData(lamMin=lamMin, lamMax=lamMax, nLam=nLam)
            
            inorm = min(range(len(wl)), key=lambda i: abs(wl[i]-lamNorm))
            data_array/=data_array[inorm]
        else:
            wl, data_array = data.getFilterData()
            print wl[0], wl[-1]
    
        ax.plot(wl, data_array, color=color)

        ax.set_xlabel('wavelength (Angstroms)', fontsize=24)
        if (isSED):
            ax.set_ylabel('flux (f$_\lambda$)', fontsize=24)
            ax.set_yscale('log')
        else:
            ax.set_ylabel('transmission', fontsize=24)
        
        
        
def getNames(filter_or_sed_dict):
    """Return unordered list of all filter names or SED names in the dictionary supplied 
    """
    return filter_or_sed_dict.keys()
    

# Retire this function, use the below 'orderFiltersByLamEff' instead
#def getFilterList(listOfFiltersFile, pathToFile="../filter_data/"):
#    """Read file containing list of filters to read, place the filter names into a list
#       Order filters listed in file is preserved                                                           
#    """
#    return np.genfromtxt(os.path.join(pathToFile,listOfFiltersFile),'str')
    
    
def orderFiltersByLamEff(filterDict):
    """Order the filters in a dictionary by their effective wavelength. Returns list of str's of filter names
       sorted in order of ascending effective wavelenght
       
       @param filterDict    dictionary of filters: keyword=filter name, value = Filter object
    """
    
    # process each entry to get filter name and calculate effective wavelength
    filter_order = []
    filter_effWave = []
    for filtname, filt in filterDict.items():
        
        filter_order.append(filtname)
        filter_effWave.append(filt.getFilterEffectiveWL())
                
    # sort based upon effective wavlength
    filter_order = [name for (lam,name) in sorted(zip(filter_effWave, filter_order))]
    
    return filter_order
    

def writeSEDs(sedDict, outputfile, minWavelen=2999., maxWavelen=12000., nWavelen=10000):
    """
    """

    spectra = []    
    for ised, (sedname, spec) in enumerate(sedDict.items()):
    
        
        print "On SED", ised+1 ,"of", len(sedDict), sedname
    
        # re-grid SEDs onto same wavelengths
        waveLen, fl = spec.getSedData(lamMin=minWavelen, lamMax=maxWavelen, nLam=nWavelen)
        
        if (ised<1):
            spectra.append(waveLen)
        
        # normalise so they sum to 1
        norm = np.sum(fl)
        spectra.append(fl/norm)
        
    spectra = np.asarray(spectra)   
        
    # write to file
    np.savetxt(outputfile, spectra)




