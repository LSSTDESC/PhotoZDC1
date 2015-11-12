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
# wltrans = transition wavelength from fine (dwfine) to coarse (dwcoarse) wavelength resolution (in Ang.) 
msed = sedFilter.MaskSEDs(brownSEDs, emission_lines_file, wltrans=10000, dwfine=0.1, dwcoarse=100.)

# return SEDs without lines masked
unmasked_seds = msed.return_unmasked_SEDs()

# return SEDs with lines masked
masked_seds = msed.return_masked_SEDs()

# return wavelength grid common to both sets of SEDs abvoe
wlgrid = msed.return_wl_grid()

# Plot first 10 spectra with and without masking
nMax = 10
for i in range(nMax):

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.plot(wlgrid, masked_seds[i, :], color='r', label="masked")
    ax.plot(wlgrid, unmasked_seds[i, :], color='black', linestyle='dotted', label="original")
    ax.set_xlabel('Wavelength (Angstroms)', fontsize=24)
    ax.set_ylabel('Flux', fontsize=24)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(labels, loc='lower right')

    ax.set_xlim(100, 7000)
    
plt.show(block=True)

