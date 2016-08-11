"""Modified from Zongge's Liu's code based upon the 
   Beck et al method: train model relating continuum 
   PC's to log[equivalent width] of emission lines
   
   Zongge applied model trained on SDSS data to match 
   emission lines to Brown et al Atlas spectra


"""


from pylab import *
import numpy as np 
import scipy.spatial as ssp
from scipy.spatial.distance import pdist
from scipy import interp
import filereading as reader
import scipy.ndimage.filters as flt
import time


### Various readers for reading in the SDSS data ###
### This is the SDSS file made by Beck et al:  
### http://www.vo.elte.hu/papers/2015/emissionlines/data/galaxy_sample_description.txt

def read_logEW(gal_sample_file):
    """Read logEW from SDSS file"""

    # column indices of the log(EW) measured from the spectra
    inds = 19 + 4*np.arange(10)
    
    # log(EW) measured from the spectra
    ems = np.loadtxt(gal_sample_file)[:,inds]
    
    return ems
    

def read_galNames(gal_sample_file):
    """Read SDSS specObjID from file"""
    names = np.loadtxt(gal_sample_file)[:,0]
    return names
    

### Functions to process the spectra ###

def norm_rules(xhat, yhat):
    """An ancillary function that normalizes the spectra according to Beck's method
       @param xhat    wavelength
       @parma yhat    flux
    """

    # list of wl's within which to find median flux
    medlist = [[4250,4300], [4600,4800], [5400,5500], [5600,5800]]
     
     # store median flux here
    medrec = np.zeros(len(medlist))
    for jj in np.arange(len(medlist)):
          
          m = medlist[jj]
          mmin = m[0]
          mmax = m[1]
          
          # indices of xhat wl's that are between mmin and mmax
          sinds = np.logical_and((xhat>mmin) , (xhat<mmax))
          
          # return median flux in this range
          medrec[jj] = np.median(yhat[sinds])
          
    # find mean of median fluxes
    tarmean = np.mean(medrec)
     
    # new flux is 1/mean(median_fluxes)*original_flux
    yhat = 1.0/tarmean * yhat
    return yhat



def inp_spec(wl_grid_file, gal_sample_file, output_dir):
    """Process the original spectra (in same dir as gal_sample_file), 
       separate the continua and degrade their resolution. 
       Save them into folder given by output_dir
       
       @param wl_grid_file       file containing wavelength grid to be intepolated
                                 E.g. given the SDSS spectra, we may want to interpolate
                                 along the wavelength of Brown SEDs
       @param gal_sample_file    13k properties measured from SDSS spectra
       @param output_dir         directory to save processed spectra to

    """
    
    # get path to galaxy spectra data
    path_to_data = '/'.join(gal_sample_file.split('/')[:-1]) + '/'
    
    # load the galaxy ID names
    names = read_galNames(gal_sample_file)

    # total number of SDSS spectra
    nn = len(names) # 13788
     
    # load in wavelengths to interpolate on to?
    xhat = np.loadtxt(wl_grid_file)     
    nx = len(xhat) # 1099
    lx = np.max(xhat) - np.min(xhat) # 3030.44
     
    # why? something about convolving spectrum?
    sigma = 8.0/2.5*nx/lx # ~1.16
    # the mask size is 20A
     
    # ?? areas to mask ??
    filer = [[3716,3738],[4330,4350],[4851,4871],[4949,4969],
             [4997,5017],[6540,6600],[6706,6726],[6721,6741]]
     
     
    # peakws undefined! snmatrix, nsmatrix unused
    #snmatrix = np.zeros((nn, len(peakws)))
    #nsmatrix = np.zeros((nn, len(peakws)))

    # loop over each SDSS spectrum
    for ii in 1+np.arange(nn):
     
        if(ii%1000==0):
            print "finishing ",100.0*float(ii)/nn,"%"

        # create name of spectrum file to read in
        fstr = 'spec_' + '%020d'%(names[ii-1]) + '.txt'
          
        # create name of spectrum file to save to
        fsavstr = 'spec_' + '%06d'%(ii) + '.txt'
          
        # full path of file to read in
        sampledir = path_to_data + fstr
        
        # open file (speed up by factor 2 using genfromtxt)
        cpspec = np.genfromtxt(sampledir)
        #with open(sampledir) as f:
        #      
        #    # read in each non-comment line
        #    lines = [line for line in f if not line.startswith('#')]
        #    FH = np.loadtxt(lines)
        #    cpspec = FH
          
        # grab flux values of spectrum in 2nd column and Gaussian filter
        # http://www.vo.elte.hu/papers/2015/emissionlines/data/spectra_description.txt
        cpspec[:,1] = flt.gaussian_filter1d(cpspec[:,1], sigma=sigma)
        
        
        # ?? loop over areas to mask ??
        for f in filer:
            fmin = f[0]
            fmax = f[1]
               
            # return indices of everything outside of mask area
            sinds = np.logical_or((cpspec[:,0]<fmin) , (cpspec[:,0]>fmax))
            cpspec = cpspec[sinds,:]     
               
            # this variable is undefined! also not used?
            #finp = finp[sinds]     

          
        # run a median filter over gaussian convolved and masked spectrum
        cpspec[:,1] = flt.median_filter(cpspec[:,1], size=10)

        yout = cpspec[:,1] # flux column (now median filtered, gaussian-convolved and masked)
        xin = cpspec[:,0]  # wavelength column from file
          
        # interpolate/extrpolate onto xhat from xin 
        # interp function from scipy, it extrapolates left and right values as constant beyond the range
        yhat = interp(xhat, xin, yout)
          
        # normalize the spectrum according to Beck's method
        yhat = norm_rules(xhat, yhat)
          
        # save the spectrum in single file
        np.savetxt(output_dir + fsavstr,
                   np.hstack((np.transpose(np.array([xhat])),np.transpose(np.array([yhat])))))

    return nn


def sum_spec(wl_grid_file, spectra_dir, nn, outfile):
    """Read the processed/normalized continua in the spectra_dir 
       and save them as one matrix n*m matrix, n is the number of 
       spectra, m is the number of wl grid points 

       @param wl_grid_file       file containing wavelength grid to be intepolated
                                 E.g. given the SDSS spectra, we may want to interpolate
                                 along the wavelength of Brown SEDs
       @param spectra_dir        directory containing normalised continua
       @param nn                 number of normalised continua
       @param outfile            file continua are saved to
    """
     
    # wavelengths spectra were interpolated to
    xhat = np.loadtxt(wl_grid_file) 
    
    # number of wavelengths
    nx = len(xhat)
     
    # spec_all array is (number of spectra, number of wavelengths)
    spec_all = np.zeros((nn,nx))
     
    # loop over spectra
    for ii in 1+np.arange(nn):
        
        if(ii%1000==0):
            print "finishing ",100.0*float(ii)/nn,"%"

        # make file name
        fstr = 'spec_' + '%06d'%(ii) + '.txt'
          
        # full file path
        sampledir = spectra_dir + fstr
          
        # open file (speed up by factor 2 using genfromtxt)
        FH = np.genfromtxt(sampledir)
        #with open(sampledir) as f:
        #      
        #    # read in non-comment lines
        #    lines = [line for line in f if not line.startswith('#')]
        #    FH = np.loadtxt(lines)
              
        # store in spec_all array
        spec_all[ii-1,:] = FH[:,1] # 2nd column is the one with the flux (1st is wl)
          
    # save array
    np.savetxt(outfile, spec_all)
     


def compute_eig_vec(continua_file, outroot, ndims=50):
    """Compute the eigenvector for the continua matrix read from continua_file

       @param continua_file     file containing n continua (rows) at m wavelengths (cols)
       @param outroot           file root to save mean spectrum, eigenvectors and PC's to
       @param ndims             number of PC's to keep
    """

    # read in the continua
    # n_spectra x n_wlgrid
    X = np.loadtxt(continua_file)
    print 'loading finished'
     
    # note this will be a shallow copy, also variable not used
    Xorg = X
     
    # mean spectrum
    xmean = np.mean(X, axis=0)
    np.savetxt(outroot + '_specMean.txt', xmean)
     
    # centerize the spectra
    X  = X - xmean
     
    print 'starting SVD'
    # U has shape (n_spectra, n_spectra)
    # V has shape (n_wlgrid, n_wlgrid)
    # s is the singular values sorted in descending order (n_wlgrid)
    # U,s,V: X = UsV^T
    U, s, V = np.linalg.svd(X, full_matrices=True)
    print 'finishing'
     
    # singular values
    # np.savetxt(outroot + '_singVals.txt', s)
     
    # v is the principle direction.
    np.savetxt(outroot + '_eigenV.txt', V)
     
    # projected centered spectra onto eigenvectors
    proj = np.dot(X, np.transpose(V))
    
    # take only first 50 PC's
    proj = proj[:,0:ndims]
    np.savetxt(outroot + "_" + str(ndims) + "PCs.txt", proj)

    return proj



### Following 3 functions are used as tools local linear regression. 
### We can also use scipy/sklearn to avoid writing our own.


def k_guassian(x, x0, c, a=1.0):
    """Define guassian kernel: exp{(x-x0)^2/2c^2}"""
    dif= x - x0
    return a * np.exp(dif * dif.T/(-2.0 * c**2))



def get_weights(inputs, datapoint, c=1.0):
    """Calculate the weight matrix for between certain datapoint 
       and the inputs data (or, the training data X)

       @param inputs       n*m matrix, n is the total number of data, m is the dimension of the data
       @param datapoint    query test datapoint
       @param c            parameter used in Gaussian kernel (not currently implemented)
       
       The weights matrix is returned
    """

    # interpret input as matrix
    x = np.mat(inputs)
    n_rows = x.shape[0]
     
    # only valid for Gaussian kernel weight
    # create a small value to avoid underflow
    # small_underflow = np.nextafter(0,1)
     
    # set up weight matrix as identity matrix
    weights = np.mat(np.eye(n_rows))
     
    for i in xrange(n_rows):
        # weights[i, i] = k_guassian(datapoint, x[i], c) + small_underflow
         
        # weight = 1/d
        weights[i,i] = 1/np.sqrt((datapoint-x[i])*(datapoint-x[i]).T)
        
    return weights



def lwr_fun(training_inputs, training_outputs, datapoint, c=1.0):
    """Local weighted regression function

       @param training_inputs      training data X
       @param training_outputs     training data Y
       @param datapoint            query test datapoint X
       @param c                    parameter used in Gaussian kernel (not currently implemented)
    """
    # c parameter not currently in use (relevant for Gaussian kernel weight)

    # diagonal matrix of weights (n_inputs, n_inputs)
    # weight is 1/distance between training input and datapoint
    weights = get_weights(training_inputs, datapoint, c=c)
    
    # n_inputs
    nn = weights.shape[0]
    
    # sum up weights and find max weight
    sumw = 0
    maxw = 0
    for i in xrange(nn):
         
         if weights[i,i]>maxw:
              maxw = weights[i,i]
              
         sumw = sumw + weights[i,i]
    
    # divide weights by their sum (so they sum to 1)
    for i in xrange(nn):
         weights[i,i] = weights[i,i]/sumw
         
    # treat as matrices
    x = np.mat(training_inputs)
    y = np.mat(training_outputs)
    
    # weighted covariance matrix
    xt = x.T*(weights*x)
    
    # if all weight in ~1 input
    if (np.max(weights)>0.999):
    
        # don't need this variable
         ind_1D = np.argmax(weights) 
         
         # return 2D index of max element in weights matrix
         ind_2D = np.unravel_index(np.argmax(weights),(nn,nn))
         
         # then just use that training input as the prediction
         yhat = y[ind_2D[0]]
    else:
        # (Moore-Penrose) pseudo-inverse of training inputs weighted covariance matrix
         xtI = np.linalg.pinv(xt)
         
         # linear regression
         betas = xtI * (x.T * (weights * y))   
         yhat = datapoint * betas

    # yhat: predicted test datapoint Y
    # sumw: not important, just sum of the weights
    return yhat, sumw



### Class for assigning emission lines ###
    
class EmLinePaste(object):

    def __init__(self, continua_pc_file, eigenV_file, mean_spec_file, 
                       gal_sample_file, wl_grid_file, npc=5):
        """Initialise emission line paster
        
           @param continua_pc_file   file containing PC's for each training spectrum
           @param eigenV_file        file containing eigenvectors from training spectra
           @param mean_spec_file     file containing mean spectrum of training spectra
           @param gal_sample_file    file containing measured data from training spectra
           @param wl_grid_file       file containing wavelength grid to place query spectrum onto
           @param npc                number of PC's to keep
        """
    
        self.continua_pc = np.loadtxt(continua_pc_file)[:,:npc]
        
        # array of ones size (n_spec, 6)
        tempPCA = np.ones((self.continua_pc.shape[0],self.continua_pc.shape[1]+1))
        # fill array up to last column with continumm PCs
        tempPCA[:,0:self.continua_pc.shape[1]] = self.continua_pc
        # basically just added column of 1's onto end of self.continua_pc array
        self.continua_pc = tempPCA
        
        
        self.eigenV = np.loadtxt(eigenV_file)
        self.mean_spec = np.loadtxt(mean_spec_file)
        self.ems = read_logEW(gal_sample_file)
        self.wl = np.loadtxt(wl_grid_file)
        self.npc = npc
        
        # kd-tree for quick nearest-neighbor continuum PC lookup
        self.mytree = ssp.KDTree(self.continua_pc)
        
        
    def add_emission_lines(self, wl, flux, kNN=30):
        """Add emission lines predicted from continua PC's to spectrum
        
           @param wl      wavelengths 
           @param flux    fluxes
           @param kNN     k-nearest neighbors to use in model
           

        """
        # re-grid onto same wavelength grid as training spectra
        flux = self._re_grid_wl(wl, flux)
        
        # get continuum PC's
        contPCs = self._get_query_spec_pcs(flux)
        
        # predict logEW's
        logEWs = self._pred_logEW(contPCs, kNN)
        
        # add logEW lines to spectrum method
        
        
    def _re_grid_wl(self, wl, flux):
        """Re-grid onto same wavelength grid as training spectra
           
           @param wl      wavelengths 
           @param flux    fluxes
        """
        
        sed_func = interp(self.wl, wl, flux)
        return sed_func
        
        
    def _get_query_spec_pcs(self, xsed):
        """Get the PC's of the query spectrum after projecting onto eigenvectors
           and keeping only npc PC's
        
           @param xsed    query spectrum
        """

        xsed  = xsed - self.mean_spec
        proj = np.dot(xsed, np.transpose(self.eigenV))[:self.npc]
        return proj
        
        
    def _pred_logEW(self, contPCs, kNN):
        """Predict logEW given continuum PC's of query spectrum
        
           @param contPCs    the 5 continuum PC's for this query spectrum
           @param kNN        number of k nearest neighbors to use in model
        """

        # add one onto end of contPCs
        tempPCA = ones((len(contPCs)+1, ))
        tempPCA[:len(contPCs)] = contPCs
        contPCs = tempPCA
        
        # query kd-tree for nearest neigbors
        # k=number of nearest neighbors to return 
        # just store indices of nearest neighbors
        # +1 because by construction first neighbor will be the query point
        _, tempind = self.mytree.query(contPCs, kNN+1)
          
        # continuum PC's of 30 nearest neighbors: TRAINING INPUTS
        nb_contPC = self.continua_pc[tempind[1:kNN+1],:]
          
        # logEW's of 30 nearest neighbors: TRAINING OUTPUTS
        nb_ems = self.ems[tempind[1:kNN+1],:]
          
        # perform local linear regression
        estem, estwt = lwr_fun(nb_contPC, nb_ems, contPCs)
        
        return estem, estwt
        
        
    def _add_line(self):
        """ I think to convert from log EW to a line you have to choose a width
            and then figure out the matching flux since log EW is the equivalent 
            width of the spectrum that would need to be completely set to zero 
            (from the continuum, I think?) in order to have the equivalent to the 
            line so, you work out the flux, then figure out what height the 
            equivalent line of width sigma = X AA would need to be to equal 
            that same amount of flux it's a goofy system it was designed for 
            absorption lines, right? so I always found it confusing when you use 
            it for emission lines though def double check me, I'm saying this 
            from memory I remember finding it hard to locate a good reference 
            for this in any textbook/article I think you're expected to pick 
            it up on the streets from people you know who process spectra
        """
        raise NotImplementedError
            
    

"""
main function for regression.
galaxy_sample.txt: a file provided by Beck which gives the logEW for SDSS galaxies
est_ems.txt: estimation for the emission lines
real_ems.txt: real value for the emission lines
test_ems.txt: Beck's estimation for the emission lines, use as a sanity check 
estwt.txt: not important, saved weight
"""

def s():

    # column indices to read in
    cols = np.arange(59)
     
    # read in properties measured from SDSS spectra
    # http://www.vo.elte.hu/papers/2015/emissionlines/data/galaxy_sample_description.txt
    datadir = 'galaxy_sample.txt'
    subsample = reader.read_file(datadir, cols)
     
    # compute the (first 50) PC's of the SDSS spectra
    contPCA = compute_eig_vec() # size (n_spec, 50)
    # take the first 5?
    contPCA = contPCA[:,0:5] # size (n_spec, 5)

    # array of ones size (n_spec, 6)
    tempPCA = np.ones((contPCA.shape[0],contPCA.shape[1]+1))
    # fill array up to last column with continumm PCs
    tempPCA[:,0:contPCA.shape[1]] = contPCA
    # basically just added column of 1's onto end of contPCA array
    contPCA = tempPCA
    
    # column indices of the log(EW) measured from the spectra
    inds = 19 + 4*np.arange(10)
    # log(EW) measured from the spectra
    ems = subsample[:,inds]
    # log(EW) estimated from the spectra by Beck et al
    testem = subsample[:,inds+1]
     
    # number of spectra?
    nn = len(contPCA)
     
    # kd-tree for quick nearest-neighbor continuum PC lookup
    mytree = ssp.KDTree(contPCA)
     
    i = 0  # unused variable?
    k = 30 # number of nearest neighbors to return
     
    estem = np.zeros((nn,10))
    estwt = np.zeros(nn)

    # loop over each spectrum
    for i in np.arange(nn):
     
        # the 5 continuum PC's for this spectrum
        tempPCA = contPCA[i,:]
          
        # query kd-tree for nearest neigbors
        # k=number of nearest neighbors to return 
        # just store indices of nearest neighbors
        # +1 because by construction first neighbor will be the query point
        _,tempind = mytree.query(tempPCA, k+1)
          
        # continuum PC's of 30 nearest neighbors: TRAINING INPUTS
        nb_contPCA = contPCA[tempind[1:k+1],:]
          
        # logEW's of 30 nearest neighbors: TRAINING OUTPUTS
        nb_ems = ems[tempind[1:k+1],:]
          
        # perform local linear regression
        estem[i],estwt[i] = lwr_fun(nb_contPCA, nb_ems, tempPCA, c=0.1)
          
        #estem[i] = rhks_predict(nb_contPCA, nb_ems, tempPCA, l=0.1)
          
        if(i%1000==0):
            print "finishing ",100.0*float(i)/nn,"%"
     
    # save estimated logEW's
    np.savetxt('est_ems.txt',estem)
     
    # save measured logEW's
    np.savetxt('real_ems.txt',ems)
     
    # save Beck et al logEW estimations
    np.savetxt('test_ems.txt',testem)
    print np.mean(estwt)










