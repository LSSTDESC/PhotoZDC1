import unittest
import os
import numpy as np
import cphotometry as cphot
import sedFilter
import cosmo
import photErrorModel as pem

# Single SED for testing
fname = os.path.join(os.path.dirname(__file__), '../sed_data/El_B2004a.sed')
seddata = np.loadtxt(fname)
TEST_SED = sedFilter.SED(seddata[:,0], seddata[:,1])

# List of filters for testing
listOfFilters = 'LSST.filters'
pathToFilters = '../filter_data'
TEST_FILTERDICT = sedFilter.createFilterDict(listOfFilters, pathToFilters)
FILTERLIST = sedFilter.orderFiltersByLamEff(TEST_FILTERDICT)

# List of SEDs for testing
SEDLIST = []
sedfiles = ['../sed_data/El_B2004a.sed','../sed_data/Sbc_B2004a.sed','../sed_data/Scd_B2004a.sed',
           '../sed_data/Im_B2004a.sed']
for name in sedfiles:
    seddata = np.loadtxt(name)
    SEDLIST.append(sedFilter.SED(seddata[:,0], seddata[:,1]))

# Cosmological model
COSMOMODEL = cosmo.cosmologyCalculator()


class TestPhotCalcs(unittest.TestCase):


    def setUp(self):
        """Load in the test SED and filters"""
        self.testsed = TEST_SED
        self.testfilters = TEST_FILTERDICT
        self.filterlist = FILTERLIST
       
       
    def tearDown(self):
        """Delete test SED and filters"""
        del self.testsed
        del self.testfilters
        
        
    def test_throw_no_filter_error(self):
        """Test throws error when filter name not in filter dictionary is supplied """
    
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        
        z = 1.
        self.assertRaises(LookupError, lambda: p.kCorrectionXY("junk", self.filterlist[0], z) )
        self.assertRaises(LookupError, lambda: p.kCorrectionXY(self.filterlist[0], "junk", z) )
        
        
    def test_throw_neg_z_error(self):
        """Test throws error with negative redshift """
        
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        z = -1.
        self.assertRaises(ValueError, lambda: p.kCorrectionXY(self.filterlist[0], self.filterlist[1], z) )
        
        
    def test_kcorrection_zero(self):
        """Test k-correction = 0 when observed-frame filter = rest-frame filter and z=0 """
        
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        self.assertEqual(p.kCorrectionXY(self.filterlist[0], self.filterlist[0], 0), 0.)
        
        
    def test_convert_mag_to_flux_positive(self):
        """Test converting any reasonable magnitude to a flux returns a positive flux"""
        
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        
        minMag = -100.
        maxMag = 100.
        nMag = 100
        dmag = (maxMag - minMag)/(nMag - 1.)
        
        for i in xrange(nMag):
            mag = minMag + i*dmag
            self.assertGreater(p.convertMagToFlux(mag, self.filterlist[0]), 0.)
            
            
    def test_convert_flux_to_mag(self):
        """Test cannot pass zero or negative flux to any flux to mag converter"""
            
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        flux = 0.
        self.assertRaises( ValueError, lambda: p.convertFluxToMag(flux, self.filterlist[0]) )
        dFluxOverFlux = 0.1
        self.assertRaises( ValueError, lambda: p.convertFluxAndErrorToMags(flux, dFluxOverFlux) )
        flux = 1.
        dFluxOverFlux = -0.1
        self.assertRaises( ValueError, lambda: p.convertFluxAndErrorToMags(flux, dFluxOverFlux) )
        
        
    def test_convert_mag_error_to_flux(self):
        """Test cannot pass negative magnitude error to flux converter"""
        
        p = cphot.PhotCalcs(self.testsed, self.testfilters)
        errorMag = -1.
        self.assertRaises( ValueError, lambda: p.convertMagErrorToFracFluxError(errorMag) )
    
        

class TestCalcMag(unittest.TestCase):


    def setUp(self):
        """Load in the test SED and filters"""
        self.testsed = TEST_SED
        self.testfilters = TEST_FILTERDICT
        self.filterlist = FILTERLIST
        self.cosmomodel = COSMOMODEL
        self.sedlist = SEDLIST
       
       
    def tearDown(self):
        """Delete test SED and filters"""
        del self.testsed
        del self.testfilters
        
   
    def test_throw_no_filter_error(self):
        """Test throws error when filter name not in filter dictionary is supplied """
        
        p = cphot.CalcMag(self.testsed, self.testfilters, self.cosmomodel)
        
        z = 1.
        absMag = (self.filterlist[0], -19)
        self.assertRaises(LookupError, lambda: p.getMag("junk", z,  absMag) )
        
        absMag = ("junk", -19)
        self.assertRaises(LookupError, lambda: p.getMag(self.filterlist[0], z,  absMag) )
        
    def test_throw_neg_z_error(self):
        """Test throws error with negative redshift """
        
        p = cphot.CalcMag(self.testsed, self.testfilters, self.cosmomodel)
        z = -1.
        absMag = (self.filterlist[0], -19)
        self.assertRaises(ValueError, lambda: p.getMag(self.filterlist[0], z, absMag) )
        
    
    def test_reasonable_magnitude_range(self):
        """Test that for a bunch of different SEDs with small/large redshifts and abs. mags a reasonable
           magnitude is calculated"""
    
        minAmag = -50.
        maxAmag = 50.
        zMin = 0.
        zMax = 10.
        
        minMag = 1e10
        maxMag = -1e10
        for sed in self.sedlist:
            
            p = cphot.CalcMag(sed, self.testfilters, self.cosmomodel)
            
            for filt in self.filterlist:
            
                absMag = (filt, minAmag)
                m1 = p.getMag(filt, zMin, absMag)
            
                absMag = (filt, maxAmag)
                m2 = p.getMag(filt, zMax, absMag)
            
                if (m1>m2):
                    raise ValueError("this shouldn't happen")
                
                if (m1<minMag):
                    minMag = m1
                if (m2<minMag):
                    minMag = m2
                if (m1>maxMag):
                    maxMag = m1
                if (m2>maxMag):
                    maxMag = m2
                
        self.assertTrue( maxMag<=float("inf") and minMag>=2.*minAmag )
        
        
        
if __name__ == '__main__':

    unittest.main()
