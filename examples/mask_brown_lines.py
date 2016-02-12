"""Produce a set of Brown spectra with major emission lines masked out

   Code takes a few minutes to run
"""


import sedFilter
import numpy as np
import itertools
import collections
import math
import matplotlib.pyplot as plt


# File containing wavelength regions to mask out
emission_lines_file = '../eml_data/emission_lines.dat'


# Read in Brown SEDs
listOfSedsFile = "brown.seds"
pathToFile = "/mnt/drive2/repos/PhotoZDC1/sed_data/"
brownSEDs = sedFilter.createSedDict(listOfSedsFile, pathToFile)


# Class that does the masking (linear interpolation across wavelength region)
msed = sedFilter.MaskSEDs(brownSEDs, emission_lines_file)
msed.mask_SEDs()


# return SEDs with lines masked
masked_seds = msed.return_masked_SEDs()


# Plot first 10 spectra with and without masking
nMax = 10

# wavelength grid
wlmin = 3000
wlmax = 12000
nlam = 10000

# normalisation wavelength
wnorm = 5500.

for i, (sedname, sed) in enumerate(brownSEDs.items()):

    sed_masked = masked_seds[sedname]
    
    wlgrid, fl_masked = sed_masked.getSedData(lamMin=wlmin, lamMax=wlmax, nLam=nlam)
    wlgrid, fl_orig = sed.getSedData(lamMin=wlmin, lamMax=wlmax, nLam=nlam)
    
    inorm = (np.abs(wlgrid-wnorm)).argmin()

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.plot(wlgrid, fl_masked/fl_masked[inorm], color='r', label="masked")
    ax.plot(wlgrid, fl_orig/fl_orig[inorm], color='black', linestyle='dotted', label="original")
    ax.set_xlabel('Wavelength (Angstroms)', fontsize=24)
    ax.set_ylabel('Flux', fontsize=24)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(labels, loc='lower right')

    #ax.set_xlim(100, 7000)
    
    if (i>nMax-2):
        break

plt.show(block=True)

