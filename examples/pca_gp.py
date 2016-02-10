"""Demo of how to generate a new SED from a continuous distribution given any set of colors via the
   PCA+Gaussian Process method

"""
import sys
import getopt
from sklearn.decomposition import PCA as sklPCA
import matplotlib.pyplot as plt
import numpy as np
import sedFilter
import photometry as phot
import sedMapper



def usage(val):
  print ''
  print ' Usage instructions: '
  print ' pca_gp.py -s save_stem -p performance_lim -c color_file -f filter_list -g covfunc,pars'
  print ''
  print ' save_stem          output will be saved to files beginning `save_stem`                             '
  print ' performance_lim    classed as good if num=`performance_lim` colors are better than LSST sys err    '
  print ' color_file         colors will be read from or saved to this file                                  '
  print ' filter_list        list of filters to use                                                          '
  print ' cov_func,pars      type of covariance function to use in GP and its parameters                     '
  print ''
  sys.exit(2)



def main(argv):
    
    save_stem = 'new_lsst' # files will be saved to filenames beginning `save_stem`
    perf_lim = 3           # performance limit: min number of colors that should reach LSST sys err
    color_file = "../tmp/brown_colors_lsst.txt"  # File to contain colors or to read colors from 
    listOfFilters = 'LSST.filters'               # Filter set to use                           
    corr_type = 'cubic'    # type of covariance function to use in GP
    theta0 = 0.2           # parameters for GP covariance function
                 
    try:
        opts, args = getopt.getopt(argv,"hs:p:c:f:g:")
    except getopt.GetoptError as err: # if include option that's not there
        usage(2)
      
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt in ("-s"):
            save_stem = arg
        elif opt in ("-p"):
            perf_lim = int(arg)
        elif opt in ("-c"):
            color_file = arg
        elif opt in ("-f"):
            listOfFilters = arg
        elif opt in ("-g"):
            corr_type = arg.split(',')[0]
            theta0 = float(arg.split(',')[1])
            
    print '\n Command line arguments:'
    print ' Saving to files ... ', save_stem
    print ' Reading/saving colors from/to file', color_file
    print ' Using', listOfFilters ,'filter set'
    print ' At least', perf_lim ,'colors must meet LSST sys err to be `good`'
    print ' Covariance function will be', corr_type ,'with parameter', theta0
    print ''


    ### Read SEDs into a dictionary
    listOfSeds = 'brown_masked.seds'                             
    pathToSEDs = '/mnt/drive2/repos/PhotoZDC1/sed_data/'
    sedDict = sedFilter.createSedDict(listOfSeds, pathToSEDs)
    nSED = len(sedDict)
    print "Number of SEDs =", nSED


    ### Filter set to calculate colors
    pathToFilters = '/mnt/drive2/repos/PhotoZDC1/filter_data/'
    filterDict = sedFilter.createFilterDict(listOfFilters, pathToFilters)
    filterList = sedFilter.orderFiltersByLamEff(filterDict)
    nFilters = len(filterList)
    print "Number of filters =", nFilters


    ### Wavelength grid to do PCA on
    minWavelen = 2999.
    maxWavelen = 12000.
    nWavelen = 10000

                 
    ### Do PCA and train GP
    ncomp = nSED
    
    pcaGP = sedMapper.PcaGaussianProc(sedDict, filterDict, color_file, ncomp, 
                                      minWavelen, maxWavelen, nWavelen,
                                      corr_type, theta0)
    colors = pcaGP._colors
    spectra = pcaGP._spectra
    waveLen = pcaGP._waveLen
    meanSpectrum = pcaGP.meanSpec
    projected_all = pcaGP.eigenvalue_coeffs
    print "... done\n"


    ### Leave out each SED in turn
    delta_mag = np.zeros((nSED,nFilters))
    perf = []
    for i, (sedname, spec) in enumerate(sedDict.items()):
    
        print "\nOn SED", i+1 ,"of", nSED
    

        ### Retrain GP with SED removed
        nc = nSED-1
        pcaGP.reTrainGP(nc, i)
    
    
        ### Reconstruct SED
        sed_rec = pcaGP.generateSpectrum(colors[i,:])
        
    
        ### Calculate colors of reconstructed SED
        pcalcs = phot.PhotCalcs(sed_rec, filterDict)
        cnt = 0
        isBad = False

        for j in range(nFilters-1):
            cs = pcalcs.computeColor(filterList[j], filterList[j+1])
        
            delta_mag[i,j] = cs-colors[i,j]
            if (j<6):
                print "(", cs, colors[i,j], delta_mag[i,j],")"
            if (abs(delta_mag[i,j])<0.005):
                cnt+=1
            if (abs(delta_mag[i,j])>0.05):
                isBad = True
        print ""


        ### Get array version of SED back
        wl, spec_rec = sed_rec.getSedData(minWavelen, maxWavelen, nWavelen)

    
        ### Plot
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        ax.plot(waveLen, spectra[i,:], color='blue', label='true')
        ax.plot(wl, spec_rec, color='red', linestyle='dashed', label='estimated')
        ax.plot(waveLen, meanSpectrum, color='black', linestyle='dotted', label='mean')
        ax.set_xlabel('wavelength (angstroms)', fontsize=24)
        ax.set_ylabel('flux', fontsize=24)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(loc='lower right', prop={'size':12})
        ax.set_title(sedname, fontsize=24) 
        
        annotate =  "Mean $\Delta$ color = {0:.5f} \n".format(np.mean(delta_mag[i,:]))
        annotate += "Stdn $\Delta$ color = {0:.5f} ".format(np.std(delta_mag[i,:]))
        y1, y2 = ax.get_ylim()
        ax.text(9000, 0.9*y2, annotate, fontsize=12)
        plt.savefig(save_stem + '_' + 'bad_' + sedname + '.png')
        #plt.show(block=True)
        
        
        ### Performance check
        print cnt,"colors within LSST systematic error"
        perf.append(cnt)
    

    perf = np.asarray(perf)

    ### Save results
    np.savetxt(save_stem + '_deltamag.txt', delta_mag)

     
    ### Plot eigenvalue 1 vs eigenvalue 2
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)    
    ax.plot(projected_all[:, 0], projected_all[:, 1], linestyle='none', marker='o', color='blue', label='good')
    ax.plot(projected_all[np.where(perf<perf_lim), 0], projected_all[np.where(perf<perf_lim), 1],
            linestyle='none', marker='o', color='red', label='bad')
    ax.set_xlabel('eigenvalue 1', fontsize=24)
    ax.set_ylabel('eigenvalue 2', fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:4], labels[:4], loc='lower right', prop={'size':12})


    ### Histogram of number of colors per SED better than LSST systematic error
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.hist(perf,20, normed=False, histtype='stepfilled')
    ax.set_xlabel('number of colors better than sys error', fontsize=24)
    plt.savefig(save_stem + '_' + 'perf.png')
    plt.show(block=True)
    

    ### Histogram of delta-mags
    for j in range(nFilters-1):
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
    
        dmag_finite = delta_mag[np.where(abs(delta_mag[:,j])<50),j].T
    
        ax.hist(dmag_finite, 20, normed=False, histtype='stepfilled')
        ax.set_xlabel('$\Delta$color$_{' + str(j) + "}$", fontsize=24)
        plt.savefig(save_stem + '_color' + str(j) + '.png')

 
    plt.show(block=True)


if __name__ == "__main__":
    main(sys.argv[1:])
