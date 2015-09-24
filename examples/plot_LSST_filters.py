"""Very simple demo program that reads in and then plots the LSST filters


"""

import sedFilter
import matplotlib.pyplot as plt

# Read in filters into a dictionary (key=filter name, value=Filter object)
listOfFiltersFile = "LSST.filters"
pathToFile = "../filter_data/"
filterDict = sedFilter.createFilterDict(listOfFiltersFile, pathToFile)

# start figure
fig = plt.figure(figsize=(20,20))

ifilt = 1
for filtname, filt in filterDict.iteritems():

    wavelengths, transmissions = filt.getFilterData()
    
    # plot
    ax = fig.add_subplot(2,3,ifilt)
    ax.plot(wavelengths, transmissions, linestyle='solid', color='black', label=filtname)
    ax.set_xlabel('wavelength (angstroms)', fontsize=24)
    ax.set_ylabel('transmission', fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(labels)
    ax.legend(loc='lower right',prop={'size':12})
    
    ifilt+=1

plt.show(block=True)
