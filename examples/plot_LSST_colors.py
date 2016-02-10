"""Example script that plots colors of a galaxy
   
"""

import os
import numpy as np
import sedFilter
import photometry as phot
import matplotlib.pyplot as plt


# filter set of choice
listOfFilters = 'LSST.filters'
pathToFilters = '../filter_data'


# creates a dictionary: keyword is filter name, value is a Filter object
filterDict = sedFilter.createFilterDict(listOfFilters, pathToFilters)


# return the filter names
filterList = sedFilter.orderFiltersByLamEff(filterDict)
print 'Filter list = ', filterList
nFilter = len(filterList)


# redshift grid to calculate colors at
zmin = 0.
zmax = 2.6
nz = 100
dz = (zmax-zmin)/(nz-1.)


# load in CWW templates
sedList = ['../sed_data/El_B2004a.sed','../sed_data/Sbc_B2004a.sed','../sed_data/Scd_B2004a.sed',
           '../sed_data/Im_B2004a.sed']


# calculate and store g-r, r-i colorw with z
zdata = np.zeros([nz,])
g_minus_r = np.zeros([nz,])
r_minus_i = np.zeros([nz,])


# load in elliptical CWW template
sedname = '../sed_data/El_B2004a.sed'
seddata = np.loadtxt(sedname)
sed = sedFilter.SED(seddata[:,0], seddata[:,1])


# instantiate photometry calculations
p = phot.PhotCalcs(sed, filterDict)


# loop over redshifts
for i in xrange(nz):

    z = zmin + i*dz
    zdata[i] = z
   
    g_minus_r[i] = p.computeColor("LSST_g", "LSST_r", z)
    r_minus_i[i] = p.computeColor("LSST_r", "LSST_i", z)
    


fig = plt.figure(figsize=(20,10))
fig.suptitle('Elliptical galaxy', fontsize=24)

# plot of g-r vs z
ax = fig.add_subplot(131)
ax.plot(zdata, g_minus_r, color='black')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$g-r$', fontsize=24)

# plot of r-i vs z
ax = fig.add_subplot(132)
ax.plot(zdata, r_minus_i, color='black')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$r-i$', fontsize=24)

# plot of g-r vs r-i
ax = fig.add_subplot(133)
sc = ax.scatter(g_minus_r, r_minus_i, 20., zdata, marker='o')
ax.set_xlabel('$g-r$', fontsize=24)
ax.set_ylabel('$r-i$', fontsize=24)
plt.colorbar(sc)

plt.show(block=True)

