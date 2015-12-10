"""Demo of how to use the cython PhotCalcs class


"""

import cphotometry as cphot # import cphotometry instead of photometry
import photometry as phot   # this would be the regular photometry class import
import numpy as np
import sedFilter
import time
import cosmo

# filter set of choice is the LSST filters
listOfFilters = 'LSST.filters'
pathToFilters = '../filter_data'


# create a dictionary: keyword is filter name, value is a Filter object
filterDict = sedFilter.createFilterDict(listOfFilters, pathToFilters)


# load in elliptical SED (1st col: wavelength in A, 2nd col: fluxes in wl units)
seddata = np.loadtxt('../sed_data/El_B2004a.sed')
testsed = sedFilter.SED(seddata[:,0], seddata[:,1])


# need cosmological model (with default parameter values)
cosmomodel = cosmo.cosmologyCalculator() 


# Notice that the API for both the python and cython modules is the same: only difference is the import
calcMagCython = cphot.CalcMag(testsed, filterDict, cosmomodel)
calcMagPython =  phot.CalcMag(testsed, filterDict, cosmomodel)


# calculate magntidude in the LSST g band for a galaxy with M_r=-20 at different redshifts
nz = 100
for i in range(nz):
    z = i*0.01

    start_time = time.time()
    cython_result = calcMagCython.getMag("LSST_g", z, ("LSST_r", -20.))
    end_time = time.time()
    cython_time = end_time - start_time

    start_time = time.time()
    python_result = calcMagPython.getMag("LSST_g", z, ("LSST_r", -20.))
    end_time = time.time()
    python_time = end_time - start_time
    
    msg = "cython/GSL = {0:.4f}, python/scipy.quad = {1:.4f}: cython is {2:.1f} times faster, "
    msg +="accurate to within {3:.4f} mags"
    
    print msg.format(cython_result, python_result, python_time/cython_time, abs(cython_result-python_result))

    
