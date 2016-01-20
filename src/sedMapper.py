"""

  Classes that return SED "parameters", or simply an SED, given input galaxy properties 
  (e.g. z, phot_sim, phys_pars)

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
import os.path
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA as sklPCA
from sklearn.gaussian_process import GaussianProcess
import numpy as np
import pandas as pd
import time

import photometry as phot
import sedFilter as sedFilter


class PcaGaussianProc(object):
    """Decompose the spectral library into eigenvalues and eigenvectors. Calculate colors for each SED in 
       given library. Train relation between each distinct combination of eigenvectors (that makes up an SED) 
       and its colors via Gaussian Process (basically a fancy regression method). Supplying a new (unseen) set
       of colors, trained method can then return the predicted set of eigenvector combinations and so a new 
       SED is reconstructed and returned
    """

    def __init__(self, sedDict, filterDict, color_file, ncomp, minWavelen=2999., maxWavelen=12000.,
                       nWavelen=10000, corr_type='cubic', theta0=0.2):
        """Initialise PCA+GP calculation 
        
           @param sedDict       dictionary of SEDs
           @param filterDict    dictionary of filters
           @param color_file    file to save SED colors to (if doesn't exist) or read colors from (if exists)
           @param ncomp         number of principal components to use
           @param minWavelen    minimum wavelength of wavelength grid
           @param maxWavelen    maximum wavelength of wavelength grid
           @param nWavelen      number of points in wavelength grid
           @param corr_type     type of covariance function (GP parameter)
           @param theta0        parameters associated with covariance function (GP parameter)
        """ 
        # first turn SED set into a set of spectra on a common wavelength grid
        # and calculate their colors within supplied filter set
        waveLen, spectra, colors = get_sed_array(sedDict, filterDict, color_file, minWavelen, maxWavelen,
                                                 nWavelen)
        
        self._waveLen = waveLen
        self._spectra = spectra
        self._colors = colors
         
        # perform PCA and train GP
        self._corr_type = corr_type
        self._theta0 = theta0
        self._doPCA(ncomp, self._spectra)
        self._trainGP(self._colors)
        
        
    def generateSpectrum(self, colors):
        """Generate a new SED object given colors supplied
        
        """
        # Generate new set of eigenvalues for supplied colors
        eigenvals_generated = self._gp.predict(colors)

        # Reconstruct SED and normalise so it sums to 1
        spec_rec = np.dot(eigenvals_generated, self.eigenspectra) + self.meanSpec
        norm = np.sum(spec_rec)
        spec_rec /= norm
        
        # Protect against zero fluxes
        spec_rec[np.where(spec_rec<0)] = 0.
        
        # Recreate SED object
        sed_rec = sedFilter.SED(self._waveLen, spec_rec)
        
        return sed_rec
        
        
    def reTrainGP(self, ncomp, index_to_remove=float('inf')):
        """For debugging/testing/checking
        
           Re-does the PCA+GP but after removing one of the spectra, the one indexed by 'index_to_remove'
           If index_to_remove=float('inf') it trains on ALL spectra, i.e. it returns the GP to 
           initially trained state
           
        """
        # indices of every SED
        all_indices = range(len(self._spectra))

        # remove index with value `index`
        colors_test = self._colors[np.delete(all_indices, index_to_remove),:]
        spectra_test = self._spectra[np.delete(all_indices, index_to_remove),:]
        
        # redo PCA
        self._doPCA(ncomp, spectra_test)
        
        # redo GP
        self._trainGP(colors_test)
        
    
    def _doPCA(self, ncomp, spectra):
        """
        """
    
        specPCA = sklPCA(ncomp)
        specPCA.fit(spectra)
        self.meanSpec = specPCA.mean_
        self.eigenspectra = specPCA.components_
        self.eigenvalue_coeffs = np.array(specPCA.transform(spectra))
        
        print "Mean spectrum shape:", self.meanSpec.shape
        print "Eigenspectra shape:", self.eigenspectra.shape
        print "Eigenvalues shape:", self.eigenvalue_coeffs.shape
        
        
    def _trainGP(self, colors):
        """
        """
        self._gp = GaussianProcess(corr = self._corr_type, theta0=self._theta0)

        # BK note:
        # Make sure we only include unique color values
        # Used method for taking unique rows in array found here:
        # http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        data = colors
        find_unique = np.ascontiguousarray(data).view(np.dtype((np.void, data.dtype.itemsize*data.shape[1])))
        unique_idx = np.unique(find_unique, return_index=True)[1]
        print "Number of unique colors", len(unique_idx), "total number of colors =", len(colors)

        # Train and predict eigenvalues for this color set
        self._gp.fit(colors[unique_idx], self.eigenvalue_coeffs[unique_idx, :])
        
        
    
        
        
    
        
        
##### Helper functions


def get_sed_array(sedDict, filterDict, color_file, minWavelen=2999., maxWavelen=12000., nWavelen=10000):
    """Return array of SEDs on same wavelength grid (along with colors as defined by filterDict and the 
       wavelength grid)
    
       @param sedDict       dictionary of SEDs
       @param filterDict    dictionary of filters
       @param color_file    file to save SED colors to or read colors from (if exists)
       @param minWavelen    minimum wavelength of wavelength grid
       @param maxWavelen    maximum wavelength of wavelength grid
       @param nWavelen      number of points in wavelength grid
       
    """
    
    # sort based upon effective wavelength
    filter_order = sedFilter.orderFiltersByLamEff(filterDict)
    
    
    # check if file exists and need to calculate colors
    isFileExist = os.path.isfile(color_file)
    if (isFileExist):
        print "\nColors already computed"
    else:
        print "\nComputing colors ... "


    # loop over each SED
    nSED = len(sedDict)
    spectra = []
    colors = []
    for ised, (sedname, spec) in enumerate(sedDict.items()):
    
        
        print "On SED", ised+1 ,"of", nSED, sedname
    
        # re-grid SEDs onto same wavelengths
        waveLen, fl = spec.getSedData(minWavelen, maxWavelen, nWavelen)
        
        
        # normalise so they sum to 1
        norm = np.sum(fl)
        spectra.append(fl/norm)
        
        
        # calculate or read colors
        cs = []
        if (isFileExist):
        
            # reading colors
            colors_in_file = np.loadtxt(color_file)
            cs = colors_in_file[ised,:]
    
        else:
        
            # calculating colors
            spec = sedFilter.SED(waveLen, fl/norm)
            pcalcs = phot.PhotCalcs(spec, filterDict)
    
            # in each filter
            for i in range(len(filterDict)-1):
                color = pcalcs.computeColor(filter_order[i], filter_order[i+1])
                if (color == float('inf')):
                    color = 99.
                cs.append(color)
        
        # store colors for this SED
        colors.append(cs)
        
    
    # conver to np arrays for ease
    spectra = np.array(spectra)
    colors = np.array(colors)
    
    
    # if had to calculate, save colors to file to re-use
    if (not isFileExist):
        print "Saving colors to file for future use"
        np.savetxt(color_file, colors)

    return waveLen, spectra, colors


def get_sed_colors(sedDict, filterDict):
    """Calculate the colors for all the SEDs in sedDict given the filters in filterDict, return as pandas
       data frame
    
       @param sedDict       dictionary of SEDs
       @param filterDict    dictionary of filters
    """
    
    ncolors = len(filterDict) - 1
    nseds = len(sedDict)

    # sort based upon effective wavelength
    filter_order = sedFilter.orderFiltersByLamEff(filterDict)
    
    # get names of colors
    color_names = []
    for i in range(ncolors):
        color_names.append(str(filter_order[i]) + "-" + str(filter_order[i+1]) )
    
    # calculate SED colors
    sed_colors = np.zeros((nseds, ncolors))
    sed_names = []
    i=0
    tot_time = 0.
    for sedname, sed in sedDict.items():
    
        print "Calculating colors for SED:", sedname
        sed_names.append(sedname)
        p = phot.PhotCalcs(sed, filterDict)
  
        start_time = time.time()
        for j in range(ncolors):
        
            sed_colors[i,j] = p.computeColor(filter_order[j], filter_order[j+1], 0.)
        end_time = time.time()
        print "Took", end_time - start_time, "to compute", ncolors, "colors"
        tot_time += (end_time - start_time)
        i+=1
    print "Total time to compute colors for SEDs =", tot_time
    
    # convert to dataframe and return
    return pd.DataFrame(sed_colors, columns=color_names, index=sed_names)



        

def check_color_match(galaxy_colors, sed_colors, nstd=3.):
    """Check typical distance between galaxy colors and colors of the SED set
    
       @param galaxy_colors   pd dataframe of galaxy colors (row=galaxy, column=color)
       @param sed_colors      pd dataframe of SED colors (row=SED, column=color: same order as galaxy_colors)
       @param nstd            number of standard deviation to check an SED is within range of
    """
    
    ncolors = galaxy_colors.shape[1]
    nseds = sed_colors.shape[0]
    print "Number of colors =", ncolors
    print "Number of SEDs =", nseds
    
    # mean and std of each galaxy color in simulation
    mean_colors = galaxy_colors.mean(axis=0)
    std_colors = galaxy_colors.std(axis=0)
    
    # poor_match is a list of SEDs that don't match each color
    poor_match = [ [] for c in range(ncolors) ]
    for i in range(ncolors):
            
            
            # define acceptable range as being within some n stds of the mean color of all galaxies
            max_color = mean_colors[i] + nstd*std_colors[i]
            min_color = mean_colors[i] - nstd*std_colors[i]
            #print "Range SED must lie within:", min_color ,"to", max_color
            
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
    print sed_colors.columns
    
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
    
    return sed_colors.index.values[sed_label]
    
    
