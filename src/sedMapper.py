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
                       nWavelen=10000, nfit=-1, corr_type='cubic', theta0=0.2):
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
        waveLen, spectra, colors = get_sed_array(sedDict, minWavelen, maxWavelen, nWavelen, 
                                                 filterDict, color_file)
        
        # in case not all eigenvalues are used in GP mapping
        #self.templateFitter = BestFitTemplateSpectrum(sedDict)
        
        self._waveLen = waveLen
        self._spectra = spectra
        self._colors = colors
         
        # parameters of GP covariance function
        self._corr_type = corr_type
        self._theta0 = theta0
        
        # perform PCA and train GP
        self._doPCA(ncomp, self._spectra)
        self._trainGP(self._colors, nfit)
        
        self._isRetrained = False
        
        
    def generateSpectrum(self, colors):
        """Generate a new SED object given colors supplied
        
        """
        # Generate new set of eigenvalues for supplied colors
        eigenvals_generated = self._gp.predict(colors)
        
        # check if not all eigenvalues were used in GP
        nGP = eigenvals_generated.shape[1]
        if (nGP < self.eigenspectra.shape[0] and not self._isRetrained):
        
            if (self._isRetrained):
                raise ValueError("ERROR! cannot have ")
            print nGP, self.eigenspectra.shape[0]
        
            # find SED with closest colors
            ifit = self._fitColors(colors)
            
            #eigenvals_generated = eigenvals_generated + self.eigenvalue_coeffs[ifit, nGP+1,:]
        
            add_on = np.reshape(self.eigenvalue_coeffs[ifit, nGP:], (1,len(self.eigenvalue_coeffs[ifit, nGP:])))
            eigenvals_generated = np.concatenate((eigenvals_generated, add_on), axis=1)
            #raise ValueError("Error! not yet implemented")

        # Reconstruct SED and normalise so it sums to 1
        spec_rec = np.dot(eigenvals_generated, self.eigenspectra) + self.meanSpec
        norm = np.sum(spec_rec)
        spec_rec /= norm
        
        # Protect against zero fluxes
        spec_rec[np.where(spec_rec<0)] = 0.
        
        # Recreate SED object
        sed_rec = sedFilter.SED(self._waveLen, spec_rec)
        
        return sed_rec
        
        
    def reTrainGP(self, ncomp, index_to_remove=float('inf'), nfit=-1):
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
        self._trainGP(colors_test, nfit)
        
        self._isRetrained = True
        
    
    def _doPCA(self, ncomp, spectra):
        """PCA the SED set in array
        
           @param ncomp    number of principle components to keep
           @param spectra  array of SEDs
        """
    
        specPCA = sklPCA(ncomp)
        specPCA.fit(spectra)
        self.meanSpec = specPCA.mean_
        self.eigenspectra = specPCA.components_
        self.eigenvalue_coeffs = np.array(specPCA.transform(spectra))
        
        print "Mean spectrum shape:", self.meanSpec.shape
        print "Eigenspectra shape:", self.eigenspectra.shape
        print "Eigenvalues shape:", self.eigenvalue_coeffs.shape
        
        
    def _trainGP(self, colors, nfit=-1):
        """Train the mapping between eigenvalues and color via Gaussian Process
        
           @param colors    array of galaxy colors
           @param nfit      number of eigenvalues to use in the mapping (if =-1 use all)
        """
        self._gp = GaussianProcess(corr = self._corr_type, theta0=self._theta0)
        
        if (nfit<0):
            nfit = self.eigenvalue_coeffs.shape[1]
            print "Using all", nfit ,"eigenvalues in GP"

        # BK note:
        # Make sure we only include unique color values
        # Used method for taking unique rows in array found here:
        # http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        data = colors
        find_unique = np.ascontiguousarray(data).view(np.dtype((np.void, data.dtype.itemsize*data.shape[1])))
        unique_idx = np.unique(find_unique, return_index=True)[1]
        print "Number of unique colors in SED set", len(unique_idx), "total number of SEDs =", len(colors)

        # Train and predict eigenvalues for this color set
        self._gp.fit(colors[unique_idx], self.eigenvalue_coeffs[unique_idx, :nfit])
        
        
    def _fitColors(self, colors):
        """Find best fit SED colors to galaxy colors
        """
        c = np.reshape(colors, (1,len(colors)))
        delta_color = self._colors - c
        delta_color_sq = delta_color*delta_color
        #mean_dcsq = np.mean(delta_color_sq, axis=1)
        #std_dcsq = np.std(delta_color_sq, axis=1)
        
        rms_colors = np.sqrt(np.sum(delta_color_sq, 1))
        
        """
        import matplotlib.pyplot as plt
        
        
        ## can make this vectorised?
        std_lim = 50.
        n_outlier = 0
        rms_colors = []
        for i in range(len(delta_color)):

            isort = np.argsort(abs(delta_color_sq[i,:]-mean_dcsq[i]))
            
            id_remove = isort[-n_outlier:][np.where(delta_color_sq[i, isort[-n_outlier:]] > std_lim*std_dcsq[i])[0]]
            ir = np.delete(range(len(isort)), id_remove)
            print len(ir)
            rms_colors.append(np.sqrt(np.sum(delta_color_sq[i, ir])/float(len(ir))))
            
            
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            ax.hist(abs(delta_color_sq[i,:]-mean_dcsq[i]))
            print std_lim*std_dcsq[i], delta_color_sq[i, isort[-10:]]
            #print len(ir) 
            #print np.max(delta_color_sq[i, ir])
            
            print "keeping", len(delta_color_sq[i,ir])
            """
        #print rms_colors
        ifit = np.argmin(rms_colors)
        #print "ifit =", ifit
        return ifit
        
        
        
        
class BestFitTemplateSpectrum(object):
    """Given a trial spectrum, find the best-fit spectrum within some template set 
       @warning: this may be garbage"""
    
    def __init__(self, sedDict):
        """
        @param sedDict    dictionary containing trial template set
        """
        self.sedDict = sedDict
        
    
    def fitSpectrum(self, waveLen, fLambda, inorm):
        """Find best-fit template spectrum to supplied spectrum fLambda(waveLen)
        
            @param waveLen    wavelength grid of supplied spectrum (in Angstroms)
            @param fLambda    fluxes (wavelength units) at each wavelength of wavelength grid
            @param inorm      index of wavelength to normalise at
        """
        # normalise so spectrum sums to 1
        fLambda /= fLambda[inorm] #np.sum(fLambda)
        
        min_rms = 1e19
        sed_fit = ""
        for i, (sedname, sed) in enumerate(self.sedDict.items()):
        
            # retrieve SED flux and normalise so sums to 1
            fl = sed.getSedData(wavelengths=waveLen)
            fl /= fl[inorm] #np.sum(fl)
            
            rms = np.sqrt(np.sum((fl - fLambda)*(fl - fLambda)))
            
            if (rms<min_rms):
                min_rms = rms
                sed_fit = sedname
                
        if (sed_fit==""):
            raise ValueError("Error!: no spectrum was fitted")
            
        #print "Spectrum", sed_fit ,"was the best fit"
        return sed_fit
        
     
        
class NearestNeighbor(object):
    """Find nth nearest neighbors in color space to galaxy given an SED library, 
       weight each appropriately and make composite spectrum.
       
       Could be extended to work in e.g. PCA eigenvalue space
    
    """

    def __init__(self, sedDict, filterDict, ipivot=2, nnn=7, nWavelen=10000, wlnorm=8000):
    #sed_lib_data, minWavelen, maxWavelen, nWavelen, nnn=7, wlnorm=8000):
        """
            @param sedDict     Dictionary of SEDs to use  
            @param filterDict  Filter set to define colors
            @param ipivot      Index of filter to reference ALL colors to (if -1 just does usual)
            @param nnn         Number of nearest neighbors to use 
            @param nWavelen    Resolution of composite spectrum
            @param wlnorm      Wavelength to normalise to 1 at
        """
        
        # Calculate colors of all SEDs in library
        self._sedDict = sedDict
        self._sed_lib_data  = get_sed_colors(self._sedDict, filterDict, ipivot, doPrinting=False)
        
        self._minWavelen, self._maxWavelen = sedFilter.getWaveLenRangeOfLib(self._sedDict)
        self._nWavelen = nWavelen
        self._nnn = nnn
        self._wlnorm = wlnorm
        
    
    def get_composite_spectrum(self, data_vector):
        """Main method: given a set of galaxy colors* return the composite spectrum 
           derived from nth nearest neighbors in the SED library
           
           @param data_vector    array of galaxy colors (*could be extended to e.g. PCA eigenvalues)
        
        """
        # first find the nearest neigbors
        neighbors = self.find_nearest_neighbors(data_vector)
        
        # then turn into composite spectrum
        waveLen, sed_comp, sed_all = self.make_composite(neighbors)
        
        return waveLen, sed_comp, sed_all
        
        
    def find_nearest_neighbors(self, data_vector):
        """Return the indices and distances of the SEDs in the SED library of the nnn nearest neighbors 
        
        """
        # check data vector has enough features
        if (len(data_vector) !=  self._sed_lib_data.shape[1]):
            raise ValueError("Error! Data vector has wrong number of features")
        
        # make data vector the same size as sed_lib
        data_vector_matrix = np.zeros(self._sed_lib_data.shape)
        for i in range(self._sed_lib_data.shape[0]):
            data_vector_matrix[i,:] = data_vector
        
        
        # subtract from each other, square, then sum across columns to get distance squared to each SED
        ddatasq = (self._sed_lib_data - data_vector_matrix)*(self._sed_lib_data - data_vector_matrix)
        distance_sq = np.sum(ddatasq, axis=1)

        # sort by ascending distance        
        isort = np.argsort(distance_sq)
        dsort = np.sort(distance_sq)
        
        ### The weight kernel is currently 1/d^2 but this can be changed
        
        # weights will be 1/distance^2, make sure they sum to 1
        suminvsq = 0.
        eps = 1e-10
        for i in range(self._nnn):
            suminvsq += 1./(distance_sq[isort[i]]+eps)
        # put result into list of tuples: first is index of neighbor SED, second is its weight
        
        neighbors = []
        for i in range(self._nnn):
            weight = 1./(distance_sq[isort[i]]+eps)
            #print weight, suminvsq
            neighbors.append( (isort[i], weight/suminvsq, distance_sq[isort[i]]) )
        
        return neighbors
        
        
    def make_composite(self, neighbors):
        """Using the list of SEDs sorted by ascending distance along with their weights, 
           create composite spectrum
        
           @param neighbors    list of: (index of SED in library, weight of SED, distance sq away)
                               sorted by ascending distance away
        """
    
        sednames = self._sedDict.keys()
        dWavelen = (self._maxWavelen-self._minWavelen)/(self._nWavelen-1.)
        
        ### get the num_nearest_neighbors SEDs
        sed_comp = np.zeros((self._nWavelen,))
        sed_all = np.zeros((self._nWavelen, self._nnn))
        for i, neighbor in enumerate(neighbors):
    
            sedid = neighbor[0]
            frac_sed = neighbor[1]
    
            # get SED
            spec = self._sedDict[sednames[sedid]]
            waveLen, fl = spec.getSedData(lamMin=self._minWavelen, lamMax=self._maxWavelen, 
                                          nLam=self._nWavelen)
        
            # normalise so flux = 1 at wlnorm
            inorm = np.argmin(abs(waveLen-self._wlnorm))
            flnorm = fl[inorm]
            fl /= flnorm
        
            # add to sed_comp
            sed_comp += fl*frac_sed
        
            # all the seds
            sed_all[:,i] = fl/np.sum(fl*dWavelen)
    
    
        # normalise so integrates to 1
        norm = np.sum(sed_comp*dWavelen)
        sed_comp/=norm
            
        return waveLen, sed_comp, sed_all
    
    
    
##### Helper functions


def get_sed_array(sedDict, minWavelen=2999., maxWavelen=12000., nWavelen=10000, 
                  filterDict=None, color_file=None):
    """Return array of SEDs on same wavelength grid (optionally along with colors as defined by filterDict)
    
       If computing colors, first orders filters by effective wavelength, then a color is:
       color_{i, i+1} = mag_filter_i - mag_filter_i+1
    
       @param sedDict       dictionary of SEDs
       @param minWavelen    minimum wavelength of wavelength grid
       @param maxWavelen    maximum wavelength of wavelength grid
       @param nWavelen      number of points in wavelength grid
       @param filterDict    dictionary of filters
       @param color_file    file to save SED colors to or read colors from (if exists)
       
    """
    
    doColors = True
    if (filterDict==None or color_file==None):
        doColors = False
    
    
    isFileExist = False
    if doColors:
    
        # sort based upon effective wavelength
        filter_order = sedFilter.orderFiltersByLamEff(filterDict)
    
    
        # check if file exists and need to calculate colors
        isFileExist = os.path.isfile(color_file)
        if (isFileExist):
            print "\nColors already computed,",
        else:
            print "\nComputing colors,",
    print "placing SEDs in array ..." 

    # loop over each SED
    nSED = len(sedDict)
    spectra = []
    colors = []
    for ised, (sedname, spec) in enumerate(sedDict.items()):
    
        
        print "On SED", ised+1 ,"of", nSED, sedname
    
        # re-grid SEDs onto same wavelengths
        waveLen, fl = spec.getSedData(lamMin=minWavelen, lamMax=maxWavelen, nLam=nWavelen)
        
        
        # normalise so they sum to 1
        norm = np.sum(fl)
        spectra.append(fl/norm)
        
        
        if doColors: 
            # calculate or read colors
            cs = []
            if (isFileExist):
        
                # reading colors
                colors_in_file = np.loadtxt(color_file)
                cs = colors_in_file[ised,:]
    
            else:
        
                # calculating colors
                spec = sedFilter.SED(waveLen, fl)#/norm)
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
    if (not isFileExist and doColors):
        print "Saving colors to file for future use"
        np.savetxt(color_file, colors)

    if doColors:
        return waveLen, spectra, colors
    else:
        return waveLen, spectra
        


def get_sed_colors(sedDict, filterDict, ipivot=-1, doPrinting=True):
    """Calculate the colors for all the SEDs in sedDict given the filters in filterDict, return as pandas
       data frame
    
       @param sedDict       dictionary of SEDs
       @param filterDict    dictionary of filters
       @param ipivot        index of filter to reference ALL colors to (if -1 just does usual)
       @param doPrinting 
    """
    
    nfilters = len(filterDict)
    ncolors = nfilters - 1
    nseds = len(sedDict)
    
    if ipivot>ncolors:
        raise ValueError("Error! pivot filter outside range")

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
    
        if doPrinting:
            print "Calculating colors for SED:", sedname
        sed_names.append(sedname)
        p = phot.PhotCalcs(sed, filterDict)
  
        start_time = time.time()
        
        if (ipivot>=0):
            # all colors in reference to a pivot filter, e.g. u-r, g-r, r-i, r-z, r-y
            #ii = 0
            colors = []
            for j in range(nfilters):
                if (j<ipivot):
                    colors.append(p.computeColor(filter_order[j], filter_order[ipivot], 0.))
                    #sed_colors[i,ii] = p.computeColor(filter_order[j], filter_order[ipivot], 0.)
                    if (i<1):
                        #print ii, 
                        print "Doing", filter_order[j] ,"-", filter_order[ipivot], 
                        print p.computeColor(filter_order[j], filter_order[ipivot], 0.)
                    #ii=+1
                elif (j>ipivot):
                    #sed_colors[i,ii] = p.computeColor(filter_order[ipivot], filter_order[j], 0.)
                    colors.append(p.computeColor(filter_order[ipivot], filter_order[j], 0.))
                    if (i<1):
                        #print ii, 
                        print "Doing", filter_order[ipivot] ,"-", filter_order[j],
                        print p.computeColor(filter_order[ipivot], filter_order[j], 0.)
                    #ii=+1
                # note that nothing is done when filter index j = ipivot    
            if (i<1):
                print colors, len(colors)
            for j in range(ncolors):
                sed_colors[i,j] = colors[j]  
        else:
            # traditional color definition: e.g. u-g, g-r, r-i etc
            for j in range(ncolors):
                sed_colors[i,j] = p.computeColor(filter_order[j], filter_order[j+1], 0.)
            
        end_time = time.time()
        
        if doPrinting:
            print "Took", end_time - start_time, "to compute", ncolors, "colors"
        
        tot_time += (end_time - start_time)
        
        i+=1
    if doPrinting:
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
    
    
