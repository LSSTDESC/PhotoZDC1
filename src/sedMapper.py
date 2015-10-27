"""

  Classes that return SED "parameters" given input galaxy properties (e.g. z, phot_sim, phys_pars)

  To use them this will be the process:

  import sedMapper

  sedMap = getattr(sedMapper, "name_of_model")
  sedMap(pos_args).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the SED mapper model to be used
  
  The  model classes must all follow the same template (will need
  to make an abstract base class to ensure this is followed). They 
  must all take the same number of arguments and have the same methods defined.
"""

from sklearn.cluster import KMeans
import numpy as np
import pandas as pd

import photometry as phot


def get_sed_colors(sedDict, filterDict):
    """Calculate the colors for all the SEDs in sedDict given the filters in filterDict
    
    """
    
    ncolors = len(filterDict) - 1
    nseds = len(sedDict)
    
    # process to get filter order and color names
    filter_order = []
    filter_effWave = []
    for filtname, filt in filterDict.items():
        filter_order.append(filtname)
        filter_effWave.append(filt.getFilterEffectiveWL())
        
    # sort based upon effective wavelength
    filter_order = [name for (lam,name) in sorted(zip(filter_effWave, filter_order))]
    
    # get names of colors
    color_names = []
    for i in range(ncolors):
        color_names.append(str(filter_order[i]) + "-" + str(filter_order[i+1]) )
    
    # calculate SED colors
    sed_colors = np.zeros((nseds, ncolors))
    i=0
    for sedname, sed in sedDict.items():
    
        print "Calculating colors for SED:", sedname
        p = phot.PhotCalcs(sed, filterDict)
  
        for j in range(ncolors):
        
            sed_colors[i,j] = p.computeColor(filter_order[j], filter_order[j+1], 0.)

        i+=1
    
    
    
    # convert to dataframe and return
    return pd.DataFrame(sed_colors, columns=color_names)


def check_color_match(galaxy_colors, sed_colors, name_mapping, nstd=3.):
    """Check typical distance between galaxy colors and colors of the SED set
    
       @param galaxy_colors   pd dataframe of galaxy colors (row=galaxy, column=color)
       @param sed_colors      pd dataframe of SED colors (row=SED, column=color: same order as galaxy_colors)
       @param name_mapping    dictionary mapping filter names to catalog column names
       @param nstd            number of standard deviation to check an SED is within range of
    """
    
    ncolors = galaxy_colors.shape[1]
    nseds = sed_colors.shape[0]
    print "Number of colors =", ncolors
    print "Number of SEDs =", nseds
    
    
    mean_colors = galaxy_colors.mean(axis=0)
    std_colors = galaxy_colors.std(axis=0)
    
    poor_match = [ [] for c in range(ncolors) ]
    for i in range(ncolors):
            
            
            # define acceptable range as being within some n stds of the mean color of all galaxies
            max_color = mean_colors[i] + nstd*std_colors[i]
            min_color = mean_colors[i] - nstd*std_colors[i]
            
            # check all SEDs are within this range and append to this color's list if not
            for j in range(nseds):
                
                if (sed_colors.iloc[j][i]>max_color or sed_colors.iloc[j][i]<min_color):
                    poor_match[i].append(j)
                
    
    # define "bad color" as one where there are no SEDs with colors within n-std of its mean color
    ibad_color = []
    i = 0
    for pm in poor_match:
        if (len(pm)==nseds):
            ibad_color.append(i) 
        i+=1
    if (len(ibad_color)>1):
        print ibad_color,"colors are a poor match to the SEDs"
    else:
        print "All colors are OK matches to SEDs"
        
    return poor_match
    
    
def perform_color_match(galaxy_colors, sed_colors, poor_match, tol=1):
    """Match each galaxy to an SED, ignoring all colors that are a bad match
    
       @param galaxy_colors   pd dataframe of galaxy colors (row=galaxy, column=color)
       @param sed_colors      pd dataframe of SED colors (row=SED, column=color: same order as galaxy_colors)
       @param poor_match      list of SEDs for each color that *don't* match that color well
       @param tol             minimum number of SEDs that a color must have an OK match with
                              tol=0: (least stringent) don't care if SED colors match galaxy colors at all
                              tol=1: at least one SED must match each galaxy color OK
                              tol=2: (more stringent) at least two SEDs must match each galaxy color OK
    """
    
    ncolors = galaxy_colors.shape[1]
    nseds = sed_colors.shape[0]
    
    # int array of which colors to select
    good_colors = []
    i = 0
    for pm in poor_match:
        if ( (nseds-len(pm))>=tol):
            good_colors.append(i)
        i+=1
    
    # std of each color
    std_colors = galaxy_colors[good_colors].std(axis=0)
    
    # One iteration of K-means will simply find the nearest SED to each galaxy color
    # normalise by std of each galaxy color to weight the distances appropriately
    color_cluster = KMeans(n_clusters=nseds, init=(sed_colors[good_colors]/std_colors), n_init=1, max_iter=1)
    color_cluster.fit((galaxy_colors[good_colors]/std_colors))
    
    # if you take labels out after one iteration then KMeans just finds all points closest to each SED color
    sed_label = color_cluster.labels_
    
    return sed_label
    
