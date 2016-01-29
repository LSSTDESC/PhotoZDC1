""" Run tests for sedFilter.Filter class
"""
import unittest
import os
import numpy as np
import sedFilter


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), '../filter_data/LSST_i.res')

class TestFilter(unittest.TestCase):


    def setUp(self):
        """Load in the test filter"""
        self.testdata = np.loadtxt(TESTDATA_FILENAME)
       
       
    def tearDown(self):
        """Delete test spectrum"""
        del self.testdata
    
    
    def test_min_wavelength(self):
        """Check the expected minimum wavelength of the filter is returned"""
        
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        lmin, lmax = s.returnFilterRange()
        self.assertEqual(self.testdata[0,0], lmin)
        
        
    def test_max_wavelength(self):
        """Check the expected maximum wavelength of the filter is returned"""
        
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        lmin, lmax = s.returnFilterRange()
        self.assertEqual(self.testdata[-1,0], lmax)
        
        
    def test_resolution_wavelength(self):
        """Check the expected wavelength resolution of the filter is returned"""
        
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        nl = s.returnFilterResolution()
        self.assertEqual(len(self.testdata[:,0]), nl)
        
        
    def test_no_negative_wl_on_read(self):
        """Check that error is thrown if a negative wavelength is read from the data"""
        
        tmpdata = np.copy(self.testdata)
        tmpdata[0,0] = -1.
        
        self.assertRaises(ValueError, lambda: sedFilter.Filter(tmpdata[:,0], tmpdata[:,1]))
        
        
    def test_no_negative_tr_on_read(self):
        """Check that error is thrown if a negative transmission is read from the data"""
        
        tmpdata = np.copy(self.testdata)
        tmpdata[0,1] = -1.
        
        self.assertRaises(ValueError, lambda: sedFilter.Filter(tmpdata[:,0], tmpdata[:,1]))
        
    
    def test_no_tr_above_1_on_read(self):
        """Check that error is thrown if a transmission greater than 1 is read from the data"""
        
        tmpdata = np.copy(self.testdata)
        tmpdata[0,1] = 1.4
        
        self.assertRaises(ValueError, lambda: sedFilter.Filter(tmpdata[:,0], tmpdata[:,1]))
        

    def test_neg_wl_response(self):
        """Check that error is thrown if a negative wavelength is given to getTrans method"""
    
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        self.assertRaises(ValueError, lambda: s.getTrans(-1.))
        
        
    def test_no_neg_trans_returned(self):
        """Over huge wavelength domain check no negative transmissions are returned"""
        
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        lamMax = 10000000000000
        nlam = 10000
        dlam = lamMax/(nlam-1)
        
        for n in xrange(nlam):
            self.assertGreaterEqual(s.getTrans(n*dlam), 0.)
            
            
    def test_no_tr_above_1_returned(self):
        """Over huge wavelength domain check no transmissions greater than 1 are returned"""
        
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        
        lamMax = 10000000000000
        nlam = 10000
        dlam = lamMax/(nlam-1)
        
        for n in xrange(nlam):
            self.assertLessEqual(s.getTrans(n*dlam), 1.)
            
            
    def test_eff_lam_reasonable(self):
        """Check effective wavelength lies within wavelength range of filter"""
                
        s = sedFilter.Filter(self.testdata[:,0], self.testdata[:,1])
        lamEff = s.getFilterEffectiveWL()
        lmin, lmax = s.returnFilterRange()
        
        self.assertTrue(lamEff>lmin and lamEff<lmax)
        

if __name__ == '__main__':
    unittest.main()

