"""Classes and functions that manipulate the SEDs and the filters

  Classes:
  SED:          holds flux vs wavelength data
  Filter:       holds transmission vs wavelength data
  EmissionLine: holds emission line model (TO BE IMPLEMENTED)
  
  Helper functions:
  createSedDict:        Reads list of SEDs from a file, then write the data into a dictionary (of SED objects)
  createFilterDict:     Reads list of filters from a file, then writes their data into a dictionary 
                        (of Filter objects)
  getFilterList:        Reads list of filters from a file and returns the name of each filter as a list
  orderFiltersByLamEff: Order the filters in a dictionary by their effective wavelength, return ordered list
                        of strings of their names
  

"""

import scipy.interpolate as interp
import scipy.integrate as integ
import numpy as np



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
        if (lam<0.):
            raise ValueError("ERROR! wavelength cannot be less than zero!")
        
        if (z<0.):
            raise ValueError("ERROR! redshift cannot be less than zero")

        lamFrame = lam/(1.+z) # does nothing if z=0 (i.e. this is the rest-frame)
        flux = self.sed(lamFrame)
        
        if (flux<0.):
            flux = 0.
            #print "Warning: negative flux was interpolated, set to zero instead"
            
        # add igm here
        if (self.igmModel != None):
            lamMeters = lam*1e-10
            flux *= self.igmModel.getObserverFrameTransmission(lamMeters, z)
        
        return flux


    def getSedData(self, lamMin, lamMax, nLam, z=0.):
        """Return wavelength and flux data for wavelength grid supplied (e.g. ready for plotting)
           
           @param lamMin    minimum wavelength in angstroms
           @param lamMax    maximum wavelength in angstroms
           @param nLam      number of points in wavelength grid
           @param z         redshift of spectrum
        """

        dLam = (lamMax - lamMin)/(nLam - 1.)
        wavelengths = []
        fluxes = []
        for i in xrange(nLam):
        
            lam = lamMin + i*dLam
            flux = self.getFlux(lam, z)
            
            wavelengths.append(lam)
            fluxes.append(flux)
            
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
        
        # check wavelength and transmission domains are valid
        if (self.lamMin<0.):
            raise ValueError("ERROR! wavelengths cannot be less than zero!")
            
        if (len(np.where(transmission<0.)[0])>0):
            raise ValueError("ERROR! cannot have negative transmission")
           
        if (len(np.where(transmission>1.)[0])>0):
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
        
        if (lam<0.):
            raise ValueError("ERROR! wavelength cannot be less than zero!")

        trans = self.filt(lam)
        
        # protect against interpolation producing negative values
        if (trans<0.):
            trans = 0.
            #print "Warning: negative transmission was interpolated, set to zero instead"
        
        return trans
        
    
    def getFilterData(self):
        """Return wavelength and transmisson data for filter (e.g. ready for plotting)

        """

        dLam = (self.lamMax - self.lamMin)/(self.nLam - 1.)
        wavelengths = []
        trans = []
        for i in xrange(self.nLam):
        
            lam = self.lamMin + i*dLam
            t = self.getTrans(lam)
            
            wavelengths.append(lam)
            trans.append(t)
            
        return wavelengths, trans
        
    
    def getFilterEffectiveWL(self):
        """Calculate effective wavelength of the filter"""
        # should this be just int F(lam)*lam dlam / int F(lam) dlam ????
        
        top = integ.quad(self._integrand1, self.lamMin, self.lamMax)[0]
        bot = integ.quad(self._integrand2, self.lamMin, self.lamMax)[0]
        lamEff = np.sqrt(top/bot)
        return lamEff
        
        
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
        

class EmissionLine(object):

    def __init__(self, emissionLinePars):
    
        self.emissionLinePars = emissionLinePars
        print "EmissionLine class not yet implemented"
        pass

        
#### Helper functions ####

def createSedDict(listOfSedsFile, pathToFile):
    """Read file containing list of SEDs to read, then read these SEDs into a dictionary
    
       Dictionary keyword is a string: name of the file that contained the filter (without path or extension)
       Dictionary value is a Filter object: the filter transmission data read from file 
       
    """
    
    f = open(pathToFile + "/" + listOfSedsFile)
        
    sedDict = {}
    for line in f:
           
        sedData = np.loadtxt(pathToFile + "/" + line.rstrip())
        sedName = line.rstrip().split('.')[0]
        print "Adding SED", sedName ,"to dictionary"
        
        sed = SED(sedData[:,0], sedData[:,1])
        sedDict[sedName] = sed
           
    return sedDict
    

def createFilterDict(listOfFiltersFile, pathToFile):
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
    

    
def getFilterList(listOfFiltersFile, pathToFile):
    """Read file containing list of filters to read, place the filter names into a list"""
    
    f = open(pathToFile + "/" + listOfFiltersFile)
        
    filterList = []
    for line in f:

        filtName = line.rstrip().split('.')[0]
        filterList.append(filtName)
           
    return filterList
       
    
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
    
    
