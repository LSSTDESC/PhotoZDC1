"""Very simple demo program that reads in and then plots the CWWK templates


"""
import numpy as np
import sedFilter
import igmModel
import matplotlib.pyplot as plt


# location and file names of CWWK SEDs
sed_location = '../sed_data/'
sed_names = ['El_B2004a', 'Sbc_B2004a', 'Scd_B2004a', 'Im_B2004a', 'SB2_B2004a', 'SB3_B2004a']
sed_file_ext = '.sed'

# wavelength grid to plot
lamMin = 1000.
lamMax = 10000.
nLam = 1000

# start figure
fig = plt.figure(figsize=(20,20))

ised = 1
for name in sed_names:
    
    # create SED object
    fname = sed_location + name + sed_file_ext
    f = np.loadtxt(fname)
    s = sedFilter.SED(f[:,0], f[:,1])
    
    # return flux data for wavelength grid
    wavelengths, fluxes = s.getSedData(lamMin=lamMin, lamMax=lamMax, nLam=nLam)
    
    # plot
    ax = fig.add_subplot(2,3,ised)
    ax.plot(wavelengths, fluxes, linestyle='solid', color='black', label=name)
    ax.set_xlabel('wavelength (angstroms)', fontsize=24)
    ax.set_ylabel('flux', fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(labels)
    ax.legend(loc='lower right',prop={'size':12})
    
    ised+=1

plt.show(block=True)
