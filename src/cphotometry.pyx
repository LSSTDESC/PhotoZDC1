"""
  
  Classes that perform photometry calculations. All magnitudes are AB system, nothing else is currently
  catered for.
  
  Concept is that each class is instantiated for a given object *type* and *observation system*. I.e.
  object type ->-> SED and observation system ->-> list of filters. This means the same instantiation can
  calculate photometry for varying other properties, e.g. redshift, absolute magnitude, object size. 

  PhotCalcs: holds base photometry calculations, k-correction, color calculation, zero point calculation,
             conversions between flux and magnitude etc

  CalcMag: inherits from PhotCalcs, adds the calculation of magnitudes via m = M + 5*log10(d) + 25 + k(z)

  ObsMag: inherits from CalcMag, adds the photometric error calculations

  Cython does the integration using GSL

"""

from cGslInteg cimport *
from libc.math cimport log10
import scipy.integrate as integ
import math
import time


### Note: had to move integrand functions outside of PhotCalcs class
### It is too hard to supply them to GSL's gsl_function due to the 
### class instance argument that comes with them as a class method

cdef double integrand1(double lam, void * params):
    """Return Sed(lam/(1+z))*Filt(lam)*lam"""
    params_as_object = <object>params
    sed, filt, redshift = params_as_object
    cdef double z = redshift
    return sed.getFlux(lam, z)*filt.getTrans(lam)*lam


cdef double integrand2(double lam, void * params):
    """Return Filt(lam)/lam"""
    params_as_object = <object>params
    filt = params_as_object[0]
    return filt.getTrans(lam)/lam

    

cdef class PhotCalcs: #(object):
    """Base photometric calculations
    
    """
    # AA: made the following definition because it works, might not be best?
    cdef object sed
    cdef object filterDict
    
    # * and ** arguments so that it can accept and ignore extra arguments
    def __cinit__(self, *args): #object sed, object filterDict, ):
        """Initialise photometry calculation
        
           @param sed           SED object (spectral energy distribution)
           @param filterDict    dictionary of filters: keyword=filter filename without path or extension,
                                value=Filter object
        
        """
        self.sed = args[0] #sed
        self.filterDict = args[1] # filterDict
        
    
    def __str__(self):
        """Return custom string representation of PhotCalcs"""
        
        str_msg ='\n  PhotCalcs Object: \n'
        str_msg += "  Contains -\n"
        str_msg += str(self.sed) 
        
        str_msg += "\n  Filters \n"
        for filt in self.filterDict:
            str_msg += "  " + filt
            str_msg += str(self.filterDict[filt])
        return str_msg
    
    
    def kCorrectionXY(self, filtX, filtY, z):
        """Compute k-correction between filter X (obs-frame) and Y (rest-frame) for the SED at redshift z
        
           @param filtX   name (string) of observed-frame filter X
           @param filtY   name (string) of rest-frame filter Y
           @param z       redshift of SED
           
           filtX and filtY must match a keyword in the dictionary filterDict
           
           returns:
        
           kXY = -2.5*log10(
                             int S(lam/(1+z))*X(lam)*lam dlam *
                             int Y(lam)/lam dlam *
                             int S(lam)*Y(lam)*lam dlam *
                             int X(lam)/lam dlam  )
                             
           @warning if the SED has zero flux within the obs-frame filter INFINITY will be returned
        """
        

        if filtX not in self.filterDict:
            emsg = "Filter " + filtX + " is not in the filter dictionary"
            raise LookupError(emsg)
        if filtY not in self.filterDict:
            emsg = "Filter " + filtY + " is not in the filter dictionary"
            raise LookupError(emsg)
            
        # unittest can't see z<1 error thrown by sedFilter.SED.getFlux()
        # probably can fix this? but for now put extra check here
        if (z<0.):
            raise ValueError("ERROR! redshift cannot be less than zero")
    
        # define integral limits by filter extents
        cdef double aX, bX, aY, bY
        aX, bX = self.filterDict[filtX].returnFilterRange()
        aY, bY = self.filterDict[filtY].returnFilterRange()
        #print 'Integrating filter X between', aX ,'and', bX
        #print 'Integrating filter Y between', aY ,'and', bY
        
        
        # these will be the integration results and their errors
        cdef double int1, int2, int3, int4
        cdef double eint1, eint2, eint3, eint4
    
        # set up the memory for the gsl integration workspace
        cdef gsl_integration_workspace * W
        W = gsl_integration_workspace_alloc(1000)
        
        # Integration precision
        cdef double epsabs = 0.
        cdef double epsrel = 1e-2
        
        # GSL function to use for integration
        cdef gsl_function F
        
        # parameters will be stored in pointer p, in a tuple of type object
        cdef object o 
        cdef void* p

        ## Integrand 1: integral of SED over observed-frame filter
        # set parameters, a tuple containing the SED, filter and redshift    
        o = (self.sed, self.filterDict[filtX], z)
        p = <void*>o
        F.function = &integrand1
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aX, bX, epsabs, epsrel, 1000, W, &int1, &eint1)
        
        
        ## Integrand 3: integral of SED over rest-frame filter
        # set parameters, a tuple containing the SED, filter and redshift = 0
        o = (self.sed, self.filterDict[filtY], 0.)
        p = <void*>o
        F.function = &integrand1
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aY, bY, epsabs, epsrel, 1000, W, &int3, &eint3)
        
        if (filtX != filtY):
        
            ## Integrand 2: integral of rest-frame filter
            # set parameters, a tuple containing just the filter
            o = (self.filterDict[filtY], )
            p = <void*>o
            F.function = &integrand2
            F.params = p
        
            # do integration
            gsl_integration_qags(&F, aY, bY, epsabs, epsrel, 1000, W, &int2, &eint2)
        
        
            ## Integrand 4: integral of observed-frame filter
            # set parameters, a tuple containing just the filter
            o = (self.filterDict[filtX], )
            p = <void*>o
            F.function = &integrand2
            F.params = p
        
            # do integration
            gsl_integration_qags(&F, aX, bX, epsabs, epsrel, 1000, W, &int4, &eint4)
        
        else:
            int2 = 1.
            int4 = 1.
        #end_time = time.time()
        #print "time to int1,int2,int3,int4 =", end_time - start_time

        
        # if there is zero flux within the observed-frame filter (e.g. SED completely redshifted out)
        if (int1==0.):
            # zero observed flux so magnitude *should* be infinite
            return float("inf")
            # print 'Zero flux'
            
        # if there is zero flux within rest-frame filter, raise error
        if (int3==0.):
            # not sure of the meaning of this, but it seems possible?
            raise ValueError("ERROR: ?? no flux within rest-frame filter??")
        
        # free the integration workspace
        gsl_integration_workspace_free(W)
        
        return -2.5*log10( 1./(1.+z) * (int1*int2)/(int3*int4) )
    
       
    ## @todo WRITE COMPUTE COLOR
    def computeColor(self, filtX, filtY, z):
        """Compute color (flux in filter X - filter Y) of SED at redshift z, return color in magnitudes
        
           @param filtX   lower wavelength filter
           @param filtY   higher wavelength filter
           @param z       redshift of SED
        """
        
        if filtX not in self.filterDict:
            emsg = "Filter " + filtX + " is not in the filter dictionary"
            raise LookupError(emsg)
        if filtY not in self.filterDict:
            emsg = "Filter " + filtY + " is not in the filter dictionary"
            raise LookupError(emsg)
        if filtX == filtY:
            raise ValueError("ERROR! cannot have color as difference between same filter")
        
        # define integral limits by filter extents
        aX, bX = self.filterDict[filtX].returnFilterRange()
        aY, bY = self.filterDict[filtY].returnFilterRange()
        
        # these will be the integration results and their errors
        cdef double int1, int2, int3, int4
        cdef double eint1, eint2, eint3, eint4
    
        # set up the memory for the gsl integration workspace
        cdef gsl_integration_workspace * W
        W = gsl_integration_workspace_alloc(1000)
        
        # Integration precision
        cdef double epsabs = 0.
        cdef double epsrel = 1e-5
        
        # GSL function to use for integration
        cdef gsl_function F
        
        # parameters will be stored in pointer p, in a tuple of type object
        cdef object o 
        cdef void* p

        ## Integrand 1: integral of SED over filter X
        # set parameters, a tuple containing the SED, filter and redshift    
        o = (self.sed, self.filterDict[filtX], z)
        p = <void*>o
        F.function = &integrand1
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aX, bX, epsabs, epsrel, 1000, W, &int1, &eint1)
        
        
        ## Integrand 2: integral of SED over filter Y
        # set parameters, a tuple containing the SED, filter and redshift = 0
        o = (self.sed, self.filterDict[filtY], z)
        p = <void*>o
        F.function = &integrand1
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aY, bY, epsabs, epsrel, 1000, W, &int2, &eint2)
        

        ## Integrand 3: integral of filter X
        # set parameters, a tuple containing just the filter
        o = (self.filterDict[filtX], )
        p = <void*>o
        F.function = &integrand2
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aX, bX, epsabs, epsrel, 1000, W, &int3, &eint3)
        
        
        ## Integrand 4: integral of filter Y
        # set parameters, a tuple containing just the filter
        o = (self.filterDict[filtY], )
        p = <void*>o
        F.function = &integrand2
        F.params = p
        
        # do integration
        gsl_integration_qags(&F, aY, bY, epsabs, epsrel, 1000, W, &int4, &eint4)
        zp = -2.5*math.log10(int4/int3)
        
        return -2.5*log10(int1/int2) + zp
    
    
    ## I'M UNSURE IF I SHOULD 'UNCALIBRATE' THE FLUX BEFORE CONVERSION??      ##
    ## don't think so as it cancels out when convert back in other direction  ##
    ## (astronomy is stupid and/or hard)                                      ##
    
    def convertMagToFlux(self, mag, filtObs=None):
        """Convert AB magnitude to flux """
        #aX, bX = self.filterDict[filtObs].returnFilterRange()
        Fc = 1. #integ.quad(self._integrand3, aX, bX, args=(filtObs))[0]
        return pow(10, -0.4*mag)*Fc
        
        
    def convertFluxToMag(self, flux, filtObs=None):
        """Convert flux to AB magnitude """
        if (flux<=0):
            raise ValueError("flux cannot be <=0")
        
        #aX, bX = self.filterDict[filtObs].returnFilterRange()
        Fc = 1. #integ.quad(self._integrand3, aX, bX, args=(filtObs))[0]
        return -2.5*math.log10(flux/Fc)
        
        
    def convertFluxAndErrorToMags(self, flux, dFluxOverFlux):
        """Convert flux and fractional flux error to AB magnitude and error """
        if (flux<=0):
            raise ValueError("flux cannot be <=0")
        if (dFluxOverFlux<0):
            raise ValueError("fractional flux error cannot be <=0")
            
        errorMag = dFluxOverFlux/(0.4*math.log(10.))
        mag = -2.5*math.log10(flux)
        return mag, errorMag
        
        
    def convertMagErrorToFracFluxError(self, errorMag):
        """Convert magnitude error into fractional flux error"""
        if (errorMag<0):
            raise ValueError("error on magnitude cannot be <=0")
        dFluxOverFlux = errorMag*(0.4*math.log(10.))
        return dFluxOverFlux
    

    

    
class CalcMag(PhotCalcs):
    """Calculate magnitudes for the given SED at redshift z, with absolute magnitude absMag in all of the
       filters in filterDict
       
       @param sed           SED object (spectral energy distribution)
       @param filterDict    dictionary of filters: keyword=filter filename without path or extension,
                            value=Filter object
       @param cosmoModel    cosmology calculator
    """
    
    #cdef object cosmoModel

    def __init__(self, sed, filterDict, cosmoModel):
    
        #PhotCalcs.__cinit__(self, sed, filterDict)
        #super(PhotCalcs, self).__init__()
        self.cosmoModel = cosmoModel


    def __str__(self):
        """Return custom string representation of CalcMag"""
        
        str_msg ='\n  CalcMag Object: \n'
        str_msg += "  Contains -\n"
        str_msg += PhotCalcs.__str__(self)
        
        str_msg += "\n  Cosmological model: \n"
        str_msg += "  " + str(self.cosmoModel)
        return str_msg
        

    def getMag(self, filtObs, z, absMag):
        """Calculate observed magnitude in filter filtObs for a galaxy at redshift z and with absolute 
           magnitude given by data in absMag
        
        @param filtObs       name of observation filter
        @param z             redshift
        @param absMag        absolute magnitude tuple (filtName, value)
        
        @warning if k-correction is infinite (i.e. by zero flux of SED being within observed-frame filter)
        the magnitude returned will also be infinite
        """
    
        filtRF = absMag[0]
        aM = absMag[1]
        
        # k-correction from rest-frame bad filtRF to observed-frame band
        # this method takes ~all the time of simulating a magnitude (~0.2s/gal)
        kc = self.kCorrectionXY(filtObs, filtRF, z)
        
        
        # luminosity distance
        self.cosmoModel.setEmissionRedShift(z)
        dL = self.cosmoModel.LuminosityDistanceMpc()
        
        # distance modulus; dL is in units of Mpc
        # 1e-5 chosen as limit because 5*log10(1e-5)+25 = 0
        mu = 0. 
        if (dL>1e-5):
            mu = 5.*math.log10(dL) + 25.
        
        # final magnitude
        mag = aM + mu + kc
        
        return mag
        

class ObsMag(CalcMag):
    """Generate observed magnitudes for the given SED at redshift z, with absolute magnitude absMag in all of 
       the filters in filterDict. Uses error model specified in errorModel
       
       @param sed           SED object (spectral energy distribution)
       @param filterDict    dictionary of filters: keyword=filter filename without path or extension,
                                value=Filter object
       @param absMag        absolute magnitude tuple (filtName, value)
       @param z             redshift
       @param cosmoModel    cosmology calculator
       @param errorModel    photometric error model object (includes obs parameters set)
    """

    def __init__(self, sed, filterDict, cosmoModel, errorModel):
    
        CalcMag.__init__(self, sed, filterDict, cosmoModel)
        self.errorModel = errorModel
        
        
    def __str__(self):
        """Return custom string representation of ObsMag"""
        
        str_msg ='\n  ObsMag Object: \n'
        str_msg += "  Contains -\n"
        str_msg += CalcMag.__str__(self)
        
        str_msg += "\n  Error model: \n"
        str_msg += "  " + str(self.errorModel)
        return str_msg
        
    def simulateObservation(self, filtObs, z, absMag, size=-1.):
        """Simulate observed magnitude in filter filtObs for a galaxy at redshift z and with absolute 
           magnitude given by data in absMag
           
        @param filtObs       name of observation filter
        @param absMag        absolute magnitude tuple (filtName, value)
        @param z             redshift
        @param size          angular extent of galaxy in arc-seconds
        """
        
     
        # true magnitude
        mag = self.getMag(filtObs, z, absMag)
        
        
        # if magnitude is infinity, return infinity for observed magnitude and error
        if (mag==float('inf')):
            return float('inf'), float('inf'), mag
        
        # get observed magnitude and error
        obsmag, emag = self.errorModel.getObs(mag, filtObs, self.sed, size)
        
        return obsmag, emag, mag
        
        
        
