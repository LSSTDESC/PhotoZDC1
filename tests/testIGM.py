import unittest
import os
import numpy as np
import genAstro
import igmModel


class TestAtomicCalcs(unittest.TestCase):
    
    
    def test_lyman_series_wavelengths(self):
        """Check the expected Lyman line wavelengths are returned"""
        
        # wavelengths of lyman lines in meters
        lymanWavelengths = 1e-9*np.array([121.6,102.6,97.3,95.0,93.8,93.1,92.6,92.3,92.1,91.9])
        
        atomicCalcs = genAstro.AtomicCalcs()
        
        for line in range(2, len(lymanWavelengths)+2):
            
            ll = atomicCalcs.getWavelengthLymanSeries(line)
            #print lymanWavelengths[line-2],ll
            self.assertAlmostEqual(lymanWavelengths[line-2],ll)
            
            
        
class TestMadauIGM(unittest.TestCase):

    def test_error_when_linemax_invalid(self):
        """Check errors are thrown in line max<2 or >5"""
        pars={}
        
        pars["lineNumStop"] = 1
        self.assertRaises(ValueError, lambda: igmModel.MadauIGM(pars))
        
        pars["lineNumStop"] = 6
        self.assertRaises(ValueError, lambda: igmModel.MadauIGM(pars))
        
        
    def test_error_thrown_negative_z(self):
        """Check error is thrown when negative redshift supplied"""
        
        m = igmModel.MadauIGM()
        lam = 3.5e-7
        z = -1.
        self.assertRaises(ValueError, lambda: m.getObserverFrameTransmission(lam, z))
        
        
    def test_error_thrown_negative_lam(self):
        """Check error is thrown when negative wavelength supplied"""
        
        m = igmModel.MadauIGM()
        lam = -3.5e-7
        z = 1.
        self.assertRaises(ValueError, lambda: m.getObserverFrameTransmission(lam, z))
        
    
    def test_tranmission_always_0_to_1(self):
        """Check transmission always between zero and 1 for decent range of variables"""
        
        for line in range(2,6):
            for lam in range(1,1000000,1000):
                for z in range(0,10,1):
                    
                    pars={}
                    pars["lineNumStop"] = line
                    m = igmModel.MadauIGM(pars)
                    
                    lam *= 1e-10 # convert from A to m
                    t = m.getObserverFrameTransmission(lam, z)
                    self.assertTrue(t>=0 and t<=1)
                    
        
if __name__ == '__main__':
    unittest.main()
