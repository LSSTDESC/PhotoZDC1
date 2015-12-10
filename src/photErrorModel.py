"""

  Classes that contain the photometric error model, e.g. can have one for point sources,
  ones for different extended source models, etc


  To use them this will be the process:

  import photErrorModel

  pErr = getattr(photErrorModel, "name_of_model")
  pErr(pos_args).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the error model to be used
  
  The photometric error model classes must all follow the same template (an abstract base class is provided
  to ensure this is followed). They must all have the same methods defined, that take the same format of
  arguments.

"""

import abc
import random
import math
from photometry import PhotCalcs



class BaseErrorModel(object):
    """Abstract base class to set the interface for creating an error model class """
    __metaclass__ = abc.ABCMeta
    
    
    def __init__(self, pars={}):
        """All derived classes should inherit this initalisation method, and not definite it themselves. It 
           checks pars is in the right format and then sets using derived class's setModelPars method """
        if (not isinstance(pars, dict)):
            raise TypeError("Error parameters passed to photometric error model not in dictionary structure")  
        self.setModelPars(pars)
    
    
    @abc.abstractmethod
    def __repr__(self):
        """Set what is printed at the REPL """
        pass
        
    @abc.abstractmethod
    def __str__(self):
        """Set what is printed when calling print with an instance """
        pass

    @abc.abstractmethod
    def setModelPars(self, pars):
        """Set the error model parameters, pars must be a dictionary """
        pass

    @abc.abstractmethod
    def getObs(self, mag, filtObs, size, sed):
        """Return the observed magnitude and error """
        pass
        
        
        
class FractionalErrorModel(BaseErrorModel, PhotCalcs):
    """Simplest error model, just adds a fractional error to the flux """
   
        
    def __repr__(self):
        return 'FractionalErrorModel with f = ' + str(self.fracError)
        
        
    def __str__(self):
        return 'FractionalErrorModel with f = ' + str(self.fracError)
        
        
    def setModelPars(self, pars):
        if "fracError" not in pars:
            self.fracError = 0.1
        else:
            self.fracError = pars["fracError"]
        
        if "minFlux" not in pars:
            self.minFlux = 2.5e-40 # equivalent to mag=99
        else: 
            self.minFlux = pars["minFlux"]
            
        
    def getObs(self, mag, filtObs=None, size=None, sed=None):
        # this should probs be a **args thing?
        
        if (mag==float('inf')):
            raise ValueError("ERROR! magnitude is infinity, cannot calculate errors")
        
        # convert mag to flux
        flux = self.convertMagToFlux(mag)
        
        # add Gaussian flux error of size fracError*flux
        errorFlux = self.fracError*flux
        fluxObs = random.gauss(mu=flux, sigma=errorFlux)
        
        # protect against negative fluxes
        if (fluxObs<self.minFlux):
            fluxObs = self.minFlux
        
        # convert flux and error back to mags
        dFluxOverFlux = errorFlux/flux
        obsmag, errorMag = self.convertFluxAndErrorToMags(fluxObs, dFluxOverFlux)
        
        return obsmag, errorMag
        
        
class LSSTErrorModel(BaseErrorModel, PhotCalcs):
    """LSST error model, median seeing conditions """
   
        
    def __repr__(self):
        printMsg = '\n LSSTErrorModel parameters:\n'
        printMsg +=' Exposure time = ' + str(self.tvis) + ' s\n'
        printMsg +=' Number of years of observations = ' + str(self.nYrObs) + '\n'
        printMsg +=' Number of visits per year per band: ' + str(self.nVisYr) + '\n'
        printMsg +=' Systematic error = ' + str(self.sigmaSys) + ' mag\n'
        printMsg +=' Airmass = ' + str(self.airMass) + '\n'
        printMsg +=' Sky brightness per band: ' + str(self.msky) + ' (mag)\n'
        printMsg +=' Seeing per band: ' + str(self.msky) + ' (arcsec)\n'
        printMsg +=' gamma per band: ' + str(self.gamma) + '\n'
        printMsg +=' Cm per band: ' + str(self.Cm) + '\n'
        printMsg +=' Extinction coeff. per band: ' + str(self.km) + '\n'
        printMsg +=' Extended source model: add ' + str(self.extendedSource) + ' mag to 5-sigma depth for'
        printMsg +=' point sources\n'
        printMsg +=' Minimum allowed observed flux = ' + str(self.minFlux) + '\n'
        return printMsg
        
        
    def __str__(self):
        printMsg = '\n LSSTErrorModel parameters:\n'
        printMsg +=' Exposure time = ' + str(self.tvis) + ' s\n'
        printMsg +=' Number of years of observations = ' + str(self.nYrObs) + '\n'
        printMsg +=' Number of visits per year per band: ' + str(self.nVisYr) + '\n'
        printMsg +=' Systematic error = ' + str(self.sigmaSys) + ' mag\n'
        printMsg +=' Airmass = ' + str(self.airMass) + '\n'
        printMsg +=' Sky brightness per band: ' + str(self.msky) + ' (mag)\n'
        printMsg +=' Seeing per band: ' + str(self.msky) + ' (arcsec)\n'
        printMsg +=' gamma per band: ' + str(self.gamma) + '\n'
        printMsg +=' Cm per band: ' + str(self.Cm) + '\n'
        printMsg +=' Extinction coeff. per band: ' + str(self.km) + '\n'
        printMsg +=' Extended source model: add ' + str(self.extendedSource) + ' mag to 5-sigma depth for'
        printMsg +=' point sources\n'
        printMsg +=' Minimum allowed observed flux = ' + str(self.minFlux) + '\n'
        return printMsg
        
        
    def setModelPars(self, pars):
        """Set all parameters by default to median observation values (table 3.2 in Science Book) """
        
        # exposure time
        if "tvis" not in pars:
            self.tvis = 30.  
        else:
            self.tvis = pars["tvis"]
        
        # expected irreducible error
        if "sigmaSys" not in pars:
            self.sigmaSys = 0.0025  
        else:
            self.sigmaSys = pars["sigmaSys"]
            
        # number of years of observations
        if "nYrObs" not in pars:
            self.nYrObs = 1.
        else:
            self.nYrObs = pars["nYrObs"]
            
        # number of visits per year
        # change to 1 in each filter to get SINGLE VISIT errors (along with having self.nYrObs=1)
        if "nVisYr" not in pars:
            self.nVisYr = {'LSST_u':6,'LSST_g':8,'LSST_r':18,'LSST_i':18,'LSST_z':16,'LSST_y':16}
        else:
            tmp = pars["nVisYr"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "nVisYr parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.nVisYr = tmp
        
        # band dependent parameter
        if "gamma" not in pars:
            self.gamma = {'LSST_u':0.037,'LSST_g':0.038,'LSST_r':0.039,'LSST_i':0.039,'LSST_z':0.040,
                          'LSST_y':0.040}
        else:
            tmp = pars["gamma"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "gamma parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.gamma = tmp
            
        # band dependent parameter
        if "Cm" not in pars:
            self.Cm = {'LSST_u':23.60,'LSST_g':24.57,'LSST_r':24.57,'LSST_i':24.47,'LSST_z':24.19,
                       'LSST_y':23.74}
        else:
            tmp = pars["Cm"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "Cm parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.Cm = tmp
            
        # sky brightness
        if "msky" not in pars:
            self.msky = {'LSST_u':21.8,'LSST_g':22.0,'LSST_r':21.3,'LSST_i':20.0,'LSST_z':19.1,
                         'LSST_y':17.5}
        else:
            tmp = pars["msky"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "msky parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.msky = tmp
            
        # seeing
        if "theta" not in pars:
            self.theta = {'LSST_u':0.77,'LSST_g':0.73,'LSST_r':0.70,'LSST_i':0.67,'LSST_z':0.65,
                          'LSST_y':0.63}
        else:
            tmp = pars["theta"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "theta parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.theta = tmp
        
        # extinction coefficient
        if "km" not in pars:
            self.km = {'LSST_u':0.48,'LSST_g':0.21,'LSST_r':0.10,'LSST_i':0.07,'LSST_z':0.06,
                       'LSST_y':0.06}
        else:
            tmp = pars["km"]
            if (not isinstance(tmp, dict) or len(tmp)!=6):
                emsg = "km parameter has wrong type (type=" + str(type(tmp)) + ") or length (length="
                emsg +=str(len(tmp)) + ")"
                raise ValueError(emsg)
            self.km = tmp
            
        # air mass
        if "airMass" not in pars:
            self.airMass = 1.2 
        else:
            self.airMass = pars["airMass"]
            
        # extended source model: simple constant added to m5 (makes m5 fainter)
        if "extendedSource" not in pars:
            self.extendedSource = 0.
        else:
            self.extendedSource = pars["extendedSource"]
        
        # set  minimum allowed flux    
        if "minFlux" not in pars:
            self.minFlux = 2.5e-40 # equivalent to mag=99
        else: 
            self.minFlux = pars["minFlux"]
        
        
    def getObs(self, mag, filtObs, size=None, sed=None):
        # this should probs be a **args thing?
        
        if (mag==float('inf')):
            raise ValueError("ERROR! magnitude is infinity, cannot calculate errors")
        
        # flux and error on flux
        flux, errorFlux = self.getFluxAndError(mag, filtObs)
        
        # add Gaussian flux error
        fluxObs = random.gauss(mu=flux, sigma=errorFlux)
        
        # protect against negative fluxes
        # this sets everything with a very tiny or negative flux to have a mag~99
        if (fluxObs<self.minFlux):
            fluxObs = self.minFlux
        
        # convert flux and error back to mags
        dFluxOverFlux = errorFlux/flux
        obsmag, errorMag = self.convertFluxAndErrorToMags(fluxObs, dFluxOverFlux)
        
        return obsmag, errorMag
        
        
    def getFluxAndError(self, mag, filtObs):
        """Return the flux, and photometric error in flux for filter filtObs"""
        
        # from LSST equations (no randomness)
        errorMag = self.getMagError(mag, filtObs)
        
        # dFluxOverFlux = errorMag*(0.4*math.log(10.))
        # because 0.4*math.log(10.) = 0.921, errorMag ~= dFluxOverFlux
        dFluxOverFlux = self.convertMagErrorToFracFluxError(errorMag) 
        
        # flux = pow(10, -0.4*mag)
        flux = self.convertMagToFlux(mag) 
        
        errorFlux = flux*dFluxOverFlux
        #print "eMag = ", errorMag , "fluxTrue =", flux , "dFlux/Flux =", dFluxOverFlux ,"eFlux =", errorFlux,
        
        return flux, errorFlux
        
        
    def getMagError(self, mag, filtObs):
        """Return the photometric error in magnitudes for filter filtObs"""
        
        # 5-sigma depth with correction for extended sources (default no correction)
        # extended sources have **shallower/brighter** 5-sigma depth
        m5 = self._m5PointSources(self.Cm[filtObs], self.msky[filtObs], self.theta[filtObs], self.tvis, 
                                 self.km[filtObs], self.airMass) - self.extendedSource
        x = self._x(mag, m5)
        nStackedObs = self.nVisYr[filtObs]*self.nYrObs
        sigmaRand = self._sigmaRand(self.gamma[filtObs], x, nStackedObs)
        sigmaTotal = self._errorTotal(self.sigmaSys, sigmaRand)
        return sigmaTotal
        
        
    def _errorTotal(self, sigmaSys, sigmaRand):
        """Total photmetric error random + systematic, eqn 3.1 in Science Book"""
        return math.sqrt(sigmaSys*sigmaSys + sigmaRand*sigmaRand)
        
        
    def _sigmaRand(self, gamma, x, nStackedObs):
        """Random photometric error as a function of magnitude per visit, eqn 3.2 in Science Book
           
           gamma depends on sky brightness, readout noise, and other factors
           x = 10^(0.4*(m-m5))
           m5 is 5-sigma depth for point sources in a given band
           nStackedObs is the number of stacked observations
        """
        sigmaRandSqSingleExp = (0.04 - gamma)*x + gamma*x*x
        sigmaRand = math.sqrt(sigmaRandSqSingleExp/nStackedObs) # actually = N/S

        # following is the full expression not included in the Science Book:
        sigmaRand = 2.5*math.log10(1. + sigmaRand)  # ~= 1.0857*N/S (stops being good approx after N/S>2

        # The Science Book eqn is only valid for high signal to noise objects which is based upon the
        # approximation that sigmaRand ~= N/S
        
        return sigmaRand
        
    def _m5PointSources(self, Cm, msky, theta, tvis, km, X):
        """5-sigma depth for point sources (all pars are band dependent except for tvis), eqn 3.3 in SciBook
        
        Cm depends on the overall throughput of the instrument (and is determined using the ETC)
        msky is the sky brightness (mag/arcsec^2)
        theta is the seeing (FWHM in arcsec)
        tvis is the exposure time in seconds (2*15s)
        km is the atmospheric extinction coefficient
        X is airmass
        """
        return Cm + 0.5*(msky-21.) + 2.5*math.log10(0.7/theta) + 1.25*math.log10(tvis/30.) - km*(X-1.)
        
        
    def _x(self, m, m5):
        """Definition of variable x"""
        return pow(10., 0.4*(m-m5))
        
        
        
