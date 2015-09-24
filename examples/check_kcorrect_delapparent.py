"""Example script that uses the k-correction calculation in photometry.PhotCalcs

   Reproduces figure 9 of de Lapparent et al 2004:
   The ESO-Sculptor Survey, Evolution at z~0.1 to 0.5 
   
"""

import os
import numpy as np
import sedFilter
import photometry as phot
import matplotlib.pyplot as plt


# filter set of choice is the Johnson filters: these may differ from exact ones used in de Lapparent
listOfFilters = 'Johnson.filters'
pathToFilters = '../filter_data'


# creates a dictionary: keyword is filter name, value is a Filter object
filterDict = sedFilter.createFilterDict(listOfFilters, pathToFilters)


# return the filter names
filterList = sedFilter.getFilterList(listOfFilters, pathToFilters)
print 'Filter list = ', filterList
nFilter = len(filterList)


# redshift grid to calculate k-corrections for
zmin = 0.
zmax = 2.6
nz = 10
dz = (zmax-zmin)/(nz-1.)


# load in CWW templates
sedList = ['../sed_data/El_B2004a.sed','../sed_data/Sbc_B2004a.sed','../sed_data/Scd_B2004a.sed',
           '../sed_data/Im_B2004a.sed']


# loop over SEDs, redshifts, filters calculating k-correction for each
zdata = np.zeros([nz,])
kcorrdata = np.zeros([len(sedList), nz, nFilter])
ised = 0
for sedname in sedList:

    # read SED data from file
    name =  sedname.split('/')[-1].split('_')[0]
    seddata = np.loadtxt(sedname)


    # turn into SED object
    sed = sedFilter.SED(seddata[:,0], seddata[:,1])


    # instantiate photometry calculations
    p = phot.PhotCalcs(sed, filterDict)


    # loop over redshifts
    for i in xrange(nz):

        z = zmin + i*dz
        zdata[i] = z
   
        # loop over BVR filters
        for j in xrange(nFilter):

            # k-correction is from filterList[j] in rest-frame of SED to filterList[j] in observed-frame of SED
            # here k-correction is from filter in observed frame to same filter in rest frame
            k = p.kCorrectionXY(filterList[j], filterList[j], z)
            kcorrdata[ised, i, j] = k
   
   
    # if want to save k-correction data to a file
    # fname = 'kcorr_sed' + name + ".txt"
    # np.savetxt(fname, kcorrdata)
    
    ised +=1


# make equivalent plot to figure 9 in de Lapparent et al

fig = plt.figure(figsize=(20,10))

# plot of R band
ax = fig.add_subplot(311)
# El
ax.plot(zdata, kcorrdata[0,:, 2], linestyle='dashed', color='black',label='E')
ax.plot(zdata, kcorrdata[1,:, 2], linestyle='dashed', color='black',label='Sbc')
ax.plot(zdata, kcorrdata[2,:, 2], linestyle='dotted', color='black',label='Scd')
ax.plot(zdata, kcorrdata[3,:, 2], linestyle='dashed', color='black',label='Mag Irr')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$k(R)$', fontsize=24)
plt.xlim([0.,2.6])
plt.ylim([-1., 5.6])
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper left',prop={'size':12})


# plot of V band
ax = fig.add_subplot(312)
# El
ax.plot(zdata, kcorrdata[0,:, 1], linestyle='dashed', color='black',label='E')
ax.plot(zdata, kcorrdata[1,:, 1], linestyle='dashed', color='black',label='Sbc')
ax.plot(zdata, kcorrdata[2,:, 1], linestyle='dotted', color='black',label='Scd')
ax.plot(zdata, kcorrdata[3,:, 1], linestyle='dashed', color='black',label='Mag Irr')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$k(V)$', fontsize=24)
plt.xlim([0.,2.6])
plt.ylim([-1., 5.6])
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper left',prop={'size':12})


# plot of B band
ax = fig.add_subplot(313)
# El
ax.plot(zdata, kcorrdata[0,:, 0], linestyle='dashed', color='black',label='E')
ax.plot(zdata, kcorrdata[1,:, 0], linestyle='dashed', color='black',label='Sbc')
ax.plot(zdata, kcorrdata[2,:, 0], linestyle='dotted', color='black',label='Scd')
ax.plot(zdata, kcorrdata[3,:, 0], linestyle='dashed', color='black',label='Mag Irr')
ax.set_xlabel('redshift', fontsize=24)
ax.set_ylabel('$k(B)$', fontsize=24)
plt.xlim([0.,2.6])
plt.ylim([-1., 5.6])
handles, labels = ax.get_legend_handles_labels()
ax.legend(labels)
ax.legend(loc='upper left',prop={'size':12})


plt.show(block=True)

