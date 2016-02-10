"""Compare E2V and ITL filters


"""

import sedFilter
import matplotlib.pyplot as plt


### Read in ITL filters
listOfFiltersFile = "LSST_ITL.filters"
pathToFile = "../filter_data/"
itl = sedFilter.createFilterDict(listOfFiltersFile, pathToFile)
itl_list = sedFilter.orderFiltersByLamEff(itl)

### Read in E2V filters
listOfFiltersFile = "LSST_E2V.filters"
pathToFile = "../filter_data/"
e2v = sedFilter.createFilterDict(listOfFiltersFile, pathToFile)
e2v_list = sedFilter.orderFiltersByLamEff(e2v)


### start figure
fig = plt.figure(figsize=(20,20))

ifilt = 1
for itlname,e2vname in zip(itl_list,e2v_list):

    wavelengths, transmissions = itl[itlname].getFilterData()
    wavelengths_e2v, transmissions_e2v = e2v[e2vname].getFilterData()
    
    # plot
    ax = fig.add_subplot(2,3,ifilt)
    ax.plot(wavelengths, transmissions, linestyle='solid', color='black', label=e2vname)
    ax.plot(wavelengths_e2v, transmissions_e2v, linestyle='dashed', color='red')
    ax.set_xlabel('wavelength (angstroms)', fontsize=24)
    ax.set_ylabel('transmission', fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(labels)
    ax.legend(loc='lower right',prop={'size':12})
    
    ifilt+=1

plt.show(block=True)
