"""
  Classes that contain the IGM model, e.g. can have one for Madau, Empirical, None


  To use them this will be the process:

  import igmModel

  igm = getattr(igmModel, "name_of_model")
  igm(pos_args).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the igm model to be used
  
  The IGM classes must all follow the same template (will need
  to make an abstract base class to ensure this is followed). They 
  must all take the same number of arguments and have the same methods defined.

"""

import abc
import numpy as np
from genAstro import AtomicCalcs


class BaseIGMModel(object):
    """Abstract base class to set the interface for creating an IGM model class """
    __metaclass__ = abc.ABCMeta
    
    
    #def __init__(self, pars={}):
    #    """All derived classes should inherit this initalisation method, and not definite it themselves. It 
    #       checks pars is in the right format and then sets using derived class's setModelPars method """
    #    if (not isinstance(pars, dict)):
    #        raise TypeError("Error parameters passed to photometric error model not in dictionary structure")  
    #    self.setModelPars(pars)
    
    def setUp(self, pars={}):
        """All derived classes should inherit this setUp method, and not definite it themselves, and include
           it in the __init__ definition. It checks pars is in the right format and then sets using derived
           class's setModelPars method """
        if (not isinstance(pars, dict)):
            raise TypeError("Error parameters passed to photometric error model not in dictionary structure")  
        self.setModelPars(pars)
    
    
    @abc.abstractmethod
    def __str__(self):
        """Set what is printed when calling print with an instance """
        pass
        
    @abc.abstractmethod
    def setModelPars(self, pars):
        """Set the IGM model parameters, pars must be a dictionary """
        pass
        
    @abc.abstractmethod
    def getObserverFrameTransmission(self, lamObs, zSource):
        """ """
        pass
        
        
    @abc.abstractmethod
    def getRestFrameTransmission(self, lamEm, zSource):
        """ """
        pass
        
        
class MadauIGM(BaseIGMModel, AtomicCalcs):
    """Madau absorption model for the IGM. See Madau 1995"""
    
    def __init__(self, pars={}):
        self.setUp(pars)
        AtomicCalcs.__init__(self)
    
    def __str__(self):
        str_msg = "  MadauIGM object:\n"
        str_msg += "  Computing Lyman series up to " + str(self.lineNumStop) + "\n"
        if (self.doLyCont):
            str_msg += "  Including Lyman continuum contribution\n"
        else: 
            str_msg += "  Not including Lyman continuum contribution\n"
        return str_msg
    
    
    def setModelPars(self, pars): 
        
        # coefficients for Lyman-alpha, beta, gamma, delta
        self.Acoeffs = [0.0036, 1.7e-3, 1.2e-3, 9.3e-4]
        
        # line numbers of each Lyman transmission
        self.lineNums = [2,3,4,5]
        
        # maximum line number taken account of
        self.lineNumMax = self.lineNums[-1]
        
        # by default compute effect of all above lines
        if "lineNumStop" not in pars:
            self.lineNumStop = self.lineNumMax
        else:
            self.lineNumStop = pars["lineNumStop"]
            if (self.lineNumStop<2 or self.lineNumStop>self.lineNumMax):
                raise ValueError("ERROR! maximum Lyman line to compute is too high or low")
            
        # by default include Lyman continuum
        if "doLyCont" not in pars:
            self.doLyCont = True
        else:
            self.doLyCont = pars["doLyCont"]
    
    
    def getObserverFrameTransmission(self, lamObs, zSource):
        """Calculate the tranmission in the observer's frame at wavelength lamObs for an object at redshift
           zSource
           
           @param lamObs   observed wavelength in meters
           @param zSource  redshift of object
        """
        tau = self._getObserverFrameOpticalDepth(lamObs, zSource)
        return np.exp(-tau)
        
        
    def getRestFrameTransmission(self, lamEm, zSource):
        """Calculate the tranmission in rest frame at wavelength lamEm for an object at redshift zSource 
        
           @param lamEm    rest-frame wavelength in meters
           @param zSource  redshift of object
        """
        
        lamObs = lamEm*(1. + zSource)
        tau = self._getObserverFrameOpticalDepth(lamObs, zSource)
        return np.exp(-tau)
        
        
    def _getObserverFrameOpticalDepth(self, lamObs, zSource):
        """Calculate the observer frame optical depth
        
           @param lamObs   observed wavelength in meters
           @param zSource  redshift of object
        """
        
        if (zSource<0.):
            raise ValueError("ERROR! redshift is negative")
        if (lamObs<0.):
            raise ValueError("ERROR! wavelength is negative")
        
        lamEm = lamObs/(1. + zSource)
        lamLyAlpha = self.getWavelengthLymanSeries(2)
        lamLyLim = self.WAVE_LYMANLIM_METERS
        lamLyLineMax = self.getWavelengthLymanSeries(self.lineNumStop)
        
        tau = 0.
        
        # Regime 1: lamEm > lamLyAlpha -> no IGM absorption
        if (lamEm>lamLyAlpha):
            #print "in Regime 1"
            return tau
        else:
        
            # pre-calc all lyman line contributions to tau (eqn 15)
            tau_lines = []
            for line in range(2, self.lineNumStop+1):
                tau_lines.append(self._tauEffLine(lamObs, line))
        
            # Regime 2: lamLyLim <= lamEm <= lamLyAlpha -> Lyman line absorption
            if (lamEm>=lamLyLim and lamEm<=lamLyAlpha):
                #print "in Regime 2"
                # eqn 15
                line = 2
                while (lamEm<=self.getWavelengthLymanSeries(line) and line<=self.lineNumStop):
                    tau += tau_lines[line-2]
                    line+=1
                    
            # Regime 3: lamEm < lamLyLim -> Lyman continuum absorption AND Lyman line absorption
            elif (self.doLyCont):  
                #print "in Regime 3"
                # eqn 15: add contribution from ALL Lyman lines
                tau += sum(tau_lines) 
                
                # eqn 16: add contribution from Lyman continuum
                zmin = lamObs/lamLyLim - 1. # min redshift an IGM absorber could have
                if (zmin<0.):
                    zmin = 0.
                tau += self._tauEffContinuum(zSource, zmin)
                
        return tau
                
                
    def _tauEffLine(self, lamObs, lineNum):
        """Equation 15 in Madau 1995: contribution of Lyman line number lineNum to optical depth
        
           @param lamObs   observed wavelength in meters
           @param lineNum  Lyman line number
        """
        
        lamLine = self.getWavelengthLymanSeries(lineNum)
        tau = self.Acoeffs[lineNum-2]*pow(lamObs/lamLine, 3.46)
        return tau
        
        
    def _tauEffContinuum(self, zSource, zAbsorber):
        """Approximate integration of equation 16 in Madau 1995: contribution of Lyman continuum absorption
        
           @param zSource    redshift of object
           @param zAbsorber  redshift of nearest Lyman cloud
        """
        
        tau = 0
        xe = 1. + zSource;
        xc = 1. + zAbsorber;
        term1 = 0.25 * (xc*xc*xc) * (pow(xe, 0.46) - pow(xc, 0.46))
        term2 = 9.4 * pow(xc, 1.5) * (pow(xe, 0.18) - pow(xc, 0.18))
        term3 = -0.7 * (xc*xc*xc) * (pow(xc, -1.32) - pow(xe, -1.32))
        term4 = -0.023 * (pow(xe, 1.68) - pow(xc, 1.68))
        tau = term1 + term2 + term3 + term4;
        return tau
        
