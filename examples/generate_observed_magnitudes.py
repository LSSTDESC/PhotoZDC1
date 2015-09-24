"""Example of how to generate observed magnitudes given a certain filter set and SED

   Performs a simplifed mini-simulation of elliptical galaxies
   Uses both photometric error model that is a simple fraction of the flux and the LSST error model


"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ
import random
import time

import sedFilter
import photometry as phot
import photErrorModel as ephot
import cosmo




#### Bunch of helper functions to simulate a redshift distribution
def nzfunc(z):
    """Simple n(z) function: n(z)~z^2*exp(-z/z0)"""
    z0 = 0.4
    return z*z*np.exp(-z/z0)
    
def normalize_nz():
    """Finds normalisation constant of n(z)~z^2*exp(-z/z0) by integrating it """
    a = 0.
    b = 10.
    norm = integ.quad(nzfunc, a, b)[0]
    return norm
    
def n_of_z(z):
    """Simple n(z) function: n(z)~z^2*exp(-z/z0), with correct normalisation"""
    return nzfunc(z)/normalize_nz()
    
def get_max_of_nofz():
    """Find maximum value of n(z) via brute force"""
    zval = 0.
    zinc = 0.01
    nzvalold = n_of_z(zval)
    
    zval += zinc
    nzval = n_of_z(zval)
    
    while (nzval>nzvalold):
    
        zval += zinc
        nzvalold = nzval
        nzval = n_of_z(zval)

    return nzval
    
def draw_z(maxval, z):
    """Draw a redshift according to n(z) probability distribution"""
    r = random.uniform(0.,1.)
    nz_at_z = n_of_z(z)
    
    
    if (r < (nz_at_z/maxval)):
       return True
    else:
        return False    


#### Real code starts here


print "\n\n\n ***Warning*** --- this program takes a few minutes to run --- \n\n"

# filter set of choice is the LSST filters
listOfFilters = 'LSST.filters'
pathToFilters = '../filter_data'


# create a dictionary: keyword is filter name, value is a Filter object
filterDict = sedFilter.createFilterDict(listOfFilters, pathToFilters)


# return the filter names
filterList = sedFilter.getFilterList(listOfFilters, pathToFilters)
print 'Filter list = ', filterList
nFilter = len(filterList)


# load in elliptical SED (1st col: wavelength in A, 2nd col: fluxes in wl units)
seddata = np.loadtxt('../sed_data/El_B2004a.sed')


# turn into SED object
sed = sedFilter.SED(seddata[:,0], seddata[:,1])


# instantiate photometry calculations

# need cosmological model
h=1.
omegamat=0.3
omegaDE=0.7
wX=-1.
wXa=0.
cosmoModel = cosmo.cosmologyCalculator(h, omegamat, omegaDE, wX, wXa)


# need error model

# fractional flux error = 10%
pars = {}
pars["fracError"] = 0.1
errorModel = ephot.FractionalErrorModel(pars)

# LSST errors, default parameters (=median observing conditions and 1yr of observations)
LsstErrorModel = ephot.LSSTErrorModel()

# Now can instantiate ObsMag class
p = phot.ObsMag(sed, filterDict, cosmoModel, errorModel)
pLSST = phot.ObsMag(sed, filterDict, cosmoModel, LsstErrorModel)


# absolute magnitude to define luminosity of galaxy
# tuple of filter absolute magnitude is defined in and the value of the magnitude
absMag = (filterList[3], -19.)


# number of galaxies to simulate, max redshift to simulate up to
nz = 500
zmax = 2.
maxnz = get_max_of_nofz() # value of n(z) at its maximum

mag = np.zeros([nz, nFilter*3+1])
magLSST = np.zeros([nz, nFilter*3+1])

numGals = 0
while (numGals<nz):

   z = random.uniform(0.,zmax)
   doKeep = draw_z(maxnz, z)
   if (not doKeep):
       continue
   
   start_time = time.time()
   print 'simulating galaxy at', z
   mag[numGals,0] = z
   magLSST[numGals,0] = z
   
   for j in xrange(nFilter):
   
       # fractional flux error
       # observed magnitude, magnitude error, true magnitude
       m, em, mt = p.simulateObservation(filterList[j], z, absMag)
       
       ii = 1 + j*3
       mag[numGals,ii] = m
       mag[numGals,ii+1] = em
       mag[numGals,ii+2] = mt
       
       # LSST 1 yr error
       # observed magnitude, magnitude error, true magnitude
       m, em, mt = pLSST.simulateObservation(filterList[j], z, absMag)
       
       magLSST[numGals,ii] = m
       magLSST[numGals,ii+1] = em
       magLSST[numGals,ii+2] = mt
   
   end_time = time.time()
   print "Time to simulate galaxy in LSST filters (twice) = ", end_time-start_time, 's'
   numGals+=1
   
# if want to save magnitude data to a file
np.savetxt("obsmag_example.txt", mag)
np.savetxt("obsmag_example_LSST.txt", magLSST)


fobs = len(magLSST[magLSST[:,10]>40,10])/float(len(mag[:,10]))
print 'Fraction of galaxies not observed by LSST =', fobs

fig = plt.figure(figsize=(20,10))

# Histogram of redshift
ax = fig.add_subplot(221)
ax.hist(mag[:,0],20)
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$n(z)$', fontsize=24)

# Histogram of i band magnitudes
ax = fig.add_subplot(222)
ax.hist(mag[:,10], 20, normed=True, histtype='stepfilled', facecolor = "none", label='10% fractional errors')
ax.hist(magLSST[:,10], 20, normed=True, histtype='stepfilled', facecolor = "none", edgecolor='red', 
        label='LSST errors')
ax.set_xlabel('$m_i$', fontsize=24)
ax.set_ylabel('$n(m)$', fontsize=24)
ax.set_title('Fraction of galaxies not observed by LSST = ' + str(fobs))
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper right',prop={'size':12})

# i vs z
ax = fig.add_subplot(223)
ax.plot(mag[:,0], mag[:,10], color="black", marker='.', linestyle='none', label='10% fractional errors')
ax.plot(mag[:,0], magLSST[:,10], color="red", marker='.', linestyle='none', label='LSST errors')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$m_i$', fontsize=24)
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper right',prop={'size':12})
plt.ylim([15., 35.])


# g-r vs r-i
ax = fig.add_subplot(224)
ax.plot(mag[:,4]-mag[:,7], mag[:,7]-mag[:,10], color="black", marker='.', linestyle='none', 
        label='10% fractional errors')
ax.plot(magLSST[:,4]-magLSST[:,7], magLSST[:,7]-magLSST[:,10], color="red", marker='.', linestyle='none',
        label='LSST errors')
ax.set_xlabel('$g-r$', fontsize=24)
ax.set_ylabel('$r-i$', fontsize=24)
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper right',prop={'size':12})
plt.xlim([0, 3.])
plt.ylim([-0.5, 3.])


plt.show(block=True)
