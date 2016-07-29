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


"""

import scipy.integrate as integ
import scipy.interpolate as interp
import math
import time
import numpy as np

#import sedFilter

#FAST_INTEG = True

class PhotCalcs(object):
    """Base photometric calculations
    
    """

    def __init__(self, sed, filterDict, FAST_INTEG=True, INTEG_PREC=10000):
        """Initialise photometry calculation
        
           Not that if FAST_INTEG=True then INTEG_PREC>=10000 for LSST precision
        
           @param sed           SED object (spectral energy distribution)
           @param filterDict    dictionary of filters: keyword=filter filename without path or extension,
                                value=Filter object
        
        """
        self.sed = sed
        self.filterDict = filterDict
        self.cache_kcorr = {}
        self.FAST_INTEG = FAST_INTEG
        self.INTEG_PREC = INTEG_PREC # needed if FAST_INTEG
        
        
        #filter_list = sedFilter.orderFiltersByLamEff(self.filterDict)
        #for i in range(len(filter_list)-1):
        #    int3 =self.filterDict[filter_list[i]].integral2
        #    int4 =self.filterDict[filter_list[i+1]].integral2
        #    zp = -2.5*math.log10(int4/int3)
        #    print "zp between filter", filter_list[i] ,"and", filter_list[i+1] ,"=", zp
        #    
        # zp between filter SDSS_u and SDSS_g = -2.01123847391
        # zp between filter SDSS_g and SDSS_r = 0.0155296994297
        # zp between filter SDSS_r and SDSS_i = 0.358714361439
        # zp between filter SDSS_i and SDSS_z = 1.62849985131
        
    
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
            
    
        # define integral limits by filter extents
        aX, bX = self.filterDict[filtX].returnFilterRange()
        aY, bY = self.filterDict[filtY].returnFilterRange()
        #print 'Integrating filter X between', aX ,'and', bX
        #print 'Integrating filter Y between', aY ,'and', bY
        lam = self.filterDict[filtX].wavelengths
        res = self.sed.getFlux(lam, z)*self.filterDict[filtX].getTrans(lam)*lam
        integrand = interp.InterpolatedUnivariateSpline(lam, res, k=1)
        if self.FAST_INTEG:
            x = np.linspace(aX, bX, self.INTEG_PREC)
            y = integrand(x)
            int1 = integ.trapz(y,x)
        else:
            int1 = integ.quad(integrand, aX, bX)[0]

        if filtY in self.cache_kcorr.keys():
            int3 = self.cache_kcorr[filtY]
        else:
            lam = self.filterDict[filtY].wavelengths
            res = self.sed.getFlux(lam, 0)*self.filterDict[filtY].getTrans(lam)*lam
            integrand = interp.InterpolatedUnivariateSpline(lam, res, k=1)
            if self.FAST_INTEG:
                x = np.linspace(aY, bY, self.INTEG_PREC)
                y = integrand(x)
                int3 = integ.trapz(y,x)
            else:
                int3 = integ.quad(integrand, aY, bY)[0]
            self.cache_kcorr[filtY] = int3
        #start_time = time.time()
        # integral of SED over observed-frame filter
        # int S(lam/(1+z))*X(lam)*lam dlam
        #int1 = integ.quad(self._integrand1, aX, bX, args=(z, filtX))[0]
        # integral of SED over rest-frame filter
        # int S(lam)*Y(lam)*lam dlam
        #int3 = integ.quad(self._integrand1, aY, bY, args=(0., filtY))[0]
        
        
        # don't perform these integrals if filter is the same
        if (filtX != filtY):
            # integral of rest-frame filter
            # int Y(lam)/lam dlam
            int2 = self.filterDict[filtY].integral2
            
            # integral of observed-frame filter
            # int X(lam)/lam dlam
            int4 =self.filterDict[filtX].integral2
            
            # check integrals of filters: these should *definitely* never be zero
            if (int2==0. or int4==0.):
                raise ValueError("ERROR: one or more filters integrate to zero")
                
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
        
        return -2.5*math.log10( 1./(1.+z) * (int1*int2)/(int3*int4) )
    
       
    def computeColor(self, filtX, filtY, z=0):
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
        
        
        # int S(lam_obs)*X(lam)*lam dlam
        lam = self.filterDict[filtX].wavelengths
        res = self.sed.getFlux(lam, z)*self.filterDict[filtX].getTrans(lam)*lam
        integrand = interp.InterpolatedUnivariateSpline(lam, res, k=1)
        if self.FAST_INTEG:
            x = np.linspace(aX, bX, self.INTEG_PREC)
            y = integrand(x)
            int1 = integ.trapz(y,x)
        else:
            int1 = integ.quad(integrand, aX, bX)[0]
        
        
        
        # int S(lam_obs)*Y(lam)*lam dlam
        lam = self.filterDict[filtY].wavelengths
        res = self.sed.getFlux(lam, z)*self.filterDict[filtY].getTrans(lam)*lam
        integrand = interp.InterpolatedUnivariateSpline(lam, res, k=1)
        if self.FAST_INTEG:
            x = np.linspace(aY, bY, self.INTEG_PREC)
            y = integrand(x)
            int2 = integ.trapz(y,x)
        else:
            int2 = integ.quad(integrand, aY, bY)[0]
        
        
        # Not sure about this zero-point term? But colors are totally off without it
        # .integral2 = int filter(lam)/lam dlam
        int3 = self.filterDict[filtX].integral2
        int4 = self.filterDict[filtY].integral2
        
        zp = -2.5*math.log10(int4/int3)
        
        
        # zero observed flux in either filter so color should be infinite
        if (int1==0. or int2==0.):
            return float("inf")
            
        return -2.5*math.log10(int1/int2) + zp
    
    
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
    
    
    def _integrand1(self, lam, z, filtname):
        """Return Sed(lam/(1+z))*Filt(lam)*lam"""
        #if (z>2.1):
        #    print lam, lam/(1.+z), self.sed.getFlux(lam, z)
        return self.sed.getFlux(lam, z)*self.filterDict[filtname].getTrans(lam)*lam
    
    
    def _integrand2(self, lam, filtname):
        """Return Filt(lam)/lam"""
        return self.filterDict[filtname].getTrans(lam)/lam
        
        
    #def _integrand3(self, lam, filtname):
    #    """Return fAB(lam)*Filt(lam)*lam: NOT CURRENTLY USED"""
    #    fAB = 3.631e-20 # AB standard source in ergs/s/cm^2/Hz units
    #    return fAB*self.filterDict[filtname].getTrans(lam)*lam
    
    
    
class CalcMag(PhotCalcs):
    """Calculate magnitudes for the given SED at redshift z, with absolute magnitude absMag in all of the
       filters in filterDict
       
       @param sed           SED object (spectral energy distribution)
       @param filterDict    dictionary of filters: keyword=filter filename without path or extension,
                            value=Filter object
       @param cosmoModel    cosmology calculator
    """

    def __init__(self, sed, filterDict, cosmoModel):
    
        PhotCalcs.__init__(self, sed, filterDict)
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
        obsmag, emag = self.errorModel.getObs(mag, filtObs, size)
        
        return obsmag, emag, mag
        
        
        
