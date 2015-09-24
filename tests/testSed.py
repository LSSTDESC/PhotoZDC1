import unittest
import os
import numpy as np
import sedFilter


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), '../sed_data/El_B2004a.sed')

class TestSed(unittest.TestCase):
    # Can we set this up so the test loops over a bunch of different SEDs and performs these same
    # tests on all of them?


    def setUp(self):
        """Load in the test spectrum"""
        self.testdata = np.loadtxt(TESTDATA_FILENAME)
       
       
    def tearDown(self):
        """Delete test spectrum"""
        del self.testdata
    
    
    def test_min_wavelength(self):
        """Check the expected minimum wavelength of the SED is returned"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        lmin,lmax = s.returnSedRange()
        self.assertEqual(self.testdata[0,0], lmin)
        
        
    def test_max_wavelength(self):
        """Check the expected maximum wavelength of the SED is returned"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        lmin,lmax = s.returnSedRange()
        self.assertEqual(self.testdata[-1,0], lmax)
        
        
    def test_resolution_wavelength(self):
        """Check the expected wavelength resolution of the SED is returned"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        nl = s.returnSedResolution()
        self.assertEqual(len(self.testdata[:,0]), nl)
        
        
    def test_no_negative_wl_on_read(self):
        """Check that error is thrown if a negative wavelength is read from the data"""
        
        tmpdata = np.copy(self.testdata)
        tmpdata[0,0] = -1.
        
        self.assertRaises(ValueError, lambda: sedFilter.SED(tmpdata[:,0], tmpdata[:,1]))
        
        
    def test_no_negative_fl_on_read(self):
        """Check that error is thrown if a negative flux is read from the data"""
        
        tmpdata = np.copy(self.testdata)
        tmpdata[0,1] = -1.
        
        self.assertRaises(ValueError, lambda: sedFilter.SED(tmpdata[:,0], tmpdata[:,1]))


    def test_flux_units(self):
        """Check flux units are 'wl' (no other option implemented yet)"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        self.assertEqual(s.returnFluxUnits(),"wl")
        

    def test_neg_wl_response(self):
        """Check that error is thrown if a negative wavelength is given to getFlux method"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        wl = -1.
        self.assertRaises(ValueError, lambda: s.getFlux(wl))
        
        
    def test_neg_z_response(self):
        """Check that error is thrown if a negative redshift is given to getFlux method"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        z = -1.
        self.assertRaises(ValueError, lambda: s.getFlux(self.testdata[5,0], z))
        
        
    def test_no_neg_flux_returned(self):
        """Over huge wavelength, redshift domain check no negative fluxes are returned"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        
        lamMax = 10000000000000
        nlam = 10000
        dlam = lamMax/(nlam-1)
        
        dz = 0.5
        ndz = 10
        for i in xrange(ndz):
            z = i*dz
            
            for n in xrange(nlam):
                self.assertGreaterEqual(s.getFlux(n*dlam, z), 0.)
                
                
    def test_interpolation(self):
        """Check interpolation by feeding original wavelength grid back to interpolation function"""
        
        wl = self.testdata[:,0]
        fl = self.testdata[:,1]
        s = sedFilter.SED(wl, fl)
        
        for i in xrange(len(wl)):
            self.assertAlmostEqual(fl[i], s.getFlux(wl[i]))
            
        
    def test_interpolation_again(self):
        """Check interpolation by checking flux values of wls between any two wl in the grid """
        
        wl = self.testdata[:,0]
        fl = self.testdata[:,1]
        s = sedFilter.SED(wl, fl)
        
        for i in xrange(len(wl)-1):
            
            # mid point wl
            wlt = 0.5*(wl[i] + wl[i+1])
            
            # fluxes around mid point
            flBounds = [fl[i], fl[i+1]]
            flMax = max(flBounds)
            flMin = min(flBounds)

            self.assertTrue( s.getFlux(wlt)<=flMax and s.getFlux(wlt)>=flMin )

                
                
    def test_set_em_line(self):
        """Check setEmLine throws error as not implemented yet"""
        
        s = sedFilter.SED(self.testdata[:,0], self.testdata[:,1])
        s.setEmLine(4.)


if __name__ == '__main__':
    unittest.main()
