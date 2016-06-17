import pandas as pd
import numpy as np
import astropy.io.fits as fits
import h5py
import matplotlib.pyplot as plt
import plot_color



### Hardcode LSST requirements
"""
From the LSST Science Book, section 3.8 (outlier rate edited)

For this high S/N gold standard subset (i<25.3) over the redshift interval, 0< z <3, the photometric redshift requirements are:

(1) The root-mean-square scatter in photometric redshifts, sigma_z/(1+z) must be smaller than 0.05, 
    with a goal of 0.02.
    [root-mean-square/sigma_z now defined as IQR/1.34896, i.e. the approx of std for a normal distribution]
(2) The fraction of 3-sigma outliers at all redshifts must be below 10%.
    [an outlier has: whichever is smaller: abs(ez)>3-sigma or abs(ez)=0.06]
(3) The bias in e_z = (z_p-z)/(1+z) must be below 0.003 (or 0.01 for combined analyses of weak lensing and 
    baryon acoustic oscillations); the uncertainty in sigma_z/(1+z) must also be known to similar accuracy.
    [bias is defined as median]
"""
lsst_sigma_max = 0.05
lsst_sigma_goal = 0.02
lsst_bias_max = 0.003
lsst_outlier_max = 0.1


### Reader for text file
def read_txt(filename, sep):
    """Read photo-z point estimate results from text file. 
       File can be \t or ' ' separated
       Place data into pandas dataframe indexed by objid
    
    """
    if sep=='\t':
        delim_whitespace = False
    else:
        delim_whitespace = True
    
    # read text file
    df = pd.read_table(filename, delim_whitespace=delim_whitespace)

    if df.columns[0][0]=='#':
        raise ValueError("Error! text file header must NOT start with #")
    
    # set the objid column as the dataframe index    
    if "objid" not in df.columns:
        raise ValueError("Error! objid column not found in text file")
    
    df = df.set_index("objid")
    
    #print "First entries of data read in"
    #print df.head(), "\n"

    
    return df
    

### Reader for FITS binary file
def read_fits(filename):
    
    hdulist = fits.open(filename)

    cols = hdulist[1].columns.names
    data = hdulist[1].data
    
    if "objid" not in cols:
        raise ValueError("Error! objid column not found in FITS file")

    # because of the following error:
    # ValueError: Big-endian buffer not supported on little-endian compiler
    # have to unpack the FITS file and manually stack it back into an array
    # before adding to pandas dataframe
    x = data[cols[0]]
    for i in range(1, len(cols)):
        x = np.vstack((x, data[cols[i]]))

    data_array = x.T
    data_array.shape
    
    # create dataframe 
    df = pd.DataFrame(data_array, columns=cols)
    
    # set the objid column as the dataframe index    
    df = df.set_index("objid")
    
    return df
    
    
### Reader for HDF5
def read_hdf5(filename):
    raise NotImplementedError
    

### Generic reader function
def read_photoz_point_estimates(filename, sep='\t'):
    """Read photo-z point estimate results from file. 
       Place data into pandas dataframe indexed by objid.

       @param filename    name (and relative path) of file containing data
       @param delimiter   delimiter of file (relevant if text or csv)
       
       filename can be a:
       - text file (extension either .txt or .dat)
       - csv file (extension .csv)
       - FITS Binary Table (extension .fits or .fit)
       - HDF5 file (extension HDF5)
       
    """
    # parse file extension
    file_ext = filename.split('.')[-1]
    
    # choose correct reader
    if file_ext in ['dat', 'txt']:
        df = read_txt(filename, sep)
    elif file_ext == 'csv':
        df = read_txt(filename, sep=',')
    elif file_ext in ['fit','fits']:
        df = read_fits(filename)
    elif file_ext in 'hdf5':
        df = read_hdf5(filename)
    else:
        print filename
        raise ValueError("Error! unknown file type")
        
    return df
    
    
### Dataframe manipulation
def join_dataframes(photoz_df, truth_df, join_type="left"):
    """Join the photo-z estimates and truth dataframes together
    
       @param photoz_df    dataframe containing photo-z point estimates
       @param truth_df     dataframe containing true redshift data
       @param join_type    like SQL, left=left join etc
    
    """
    return photoz_df.join(truth_df, how=join_type, lsuffix="_pz", rsuffix="_tr")
    

def add_ez_column(df, true_z, photo_z="z_photo", col_name="ez"):
    """Add calculated column of ez=(zp-zt)/(1+zt) to dataframe containing photo-z estimates
    
       @param df          dataframe containing photo-z estimates
       @param true_z      either: Series containing true z or column name (str) of true z in df
       @param photo_z     column name (str) of photo-z estimates in df
       @param col_name    column name (str) to name new ez column
    
       
    """
    
    if isinstance(true_z, basestring):
        ez = (df[photo_z] - df[true_z])/(1. + df[true_z])
    else:
        ez = (df[photo_z] - true_z)/(1. + true_z)
        
    df[col_name] = ez
    
    
### Basic stats
def global_stats(photo_z, true_z, selection=None):
    """Return the standard deviation, median, IQR and IQR-sigma of ez=(zp-zs)/(1+zs)
    
       @param photo_z    Series containing point photo-z estimates
       @param true_z     Series containing true (or spec) redshifts 
       @param selection  Series containing bool values
       
       The selection is performed on both Series objects. All must contain same indices
       
       Returned:
       std[ez], median[ez], iqr[ez], iqr[ez]/1.349
       
       and if selection is not None also number of objects meeting selection criteria

    """
    if (selection is not None):
        if ( len(photo_z)!=len(true_z) or len(selection)!=len(true_z) ):
            raise ValueError("Error: photo-z/true-z/selection are unalignable")


    if (selection is None): 
        ez = (photo_z - true_z)/(1. + true_z)
    else:
        ez = (photo_z[selection] - true_z[selection])/(1. + true_z[selection])
    
    global_sigma = np.std(ez)
    global_bias = np.median(ez)
    
    # interquartile range
    x75, x25 = np.percentile(ez,[75.,25.])
    global_iqr = x75 - x25
    
    # relation between inter-quartile range and Normal sigma
    global_iqrsig = global_iqr/1.349 
    

    print "global sigma: {0:.3f} global bias: {1:.3f}".format(global_sigma, global_bias)
    print "global IQR: {0:.3f} global IQR-sigma: {1:.3f}".format(global_iqr, global_iqrsig)
    nselect = len(ez)
    if (selection is not None):
        print "number selected =", nselect, "\n"
    else:
        print "total number =", nselect, "\n"
    
    if (selection is not None):
        return global_sigma, global_bias, global_iqr, global_iqrsig, nselect
    else:
        return global_sigma, global_bias, global_iqr, global_iqrsig
    

### LSST requirement statistic definitions
def bias(ez):
    """Bias is the median e_z = (z_p-z)/(1+z)
    """
    return np.median(ez)
    
    
def sigma(ez):
    """Sigma is IQR/1.34896, which is equivalent to the 
       standard deviation for a normal distribution
    """
    x75, x25 = np.percentile(ez,[75.,25.])
    return (x75 - x25)/1.34896


def outlier(ez):
    """Outliers are galaxies with either abs(ez)>3*sigma or abs(ez)>0.06
       whichever of 3*sigma or 0.06 is smaller. The outlier rate is the 
       fraction of galaxies that are outliers
    """
    if (3.*sigma(ez)<0.06):
        outlier_limit = 3.*sigma(ez)
    else:
        outlier_limit = 0.06
    
    return len(ez[abs(ez)>outlier_limit])/float(len(ez))

def catastrophic_outlier(ez):
    """Catastrophic outliers are galaxies with abs(ez)>0.15. The 
       catastrophic outlier rate is the fraction of galaxies that 
       are catastrophic outliers
    """
    return len(ez[abs(ez)>0.15])/float(len(ez))
    

#### Output results
def output_stats_vs_redshift(outfile, df, photo_z, true_z, selection=None, zbins=[0,2,0.2]):
    """Output statistics as a function of redshift to a file
     
       @param outfile      file to output statistics to
       @param df           dataframe containing photo-z data 
       @param photo_z      column name (str) containing photo-z estimates
       @param true_z       column name (str) containing true redshifts
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
    """

    # unpack redshift binning
    zmin = zbins[0]
    zmax = zbins[1]
    dbin = zbins[2]
    zvals = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # add column that indicates redshift bin "category" for each object
    df["zbin"] = pd.cut(df[true_z], np.arange(zmin, zmax+dbin, dbin))
    
    # groupby redshift bin and calculate stats
    if (selection is not None):
        biases = df.loc[selection,:].groupby("zbin")["ez"].apply(bias)
        sigmas = df.loc[selection,:].groupby("zbin")["ez"].apply(sigma)
        outliers = df.loc[selection,:].groupby("zbin")["ez"].apply(outlier)
        cat_outliers = df.loc[selection,:].groupby("zbin")["ez"].apply(catastrophic_outlier)
    else:
        biases = df.groupby("zbin")["ez"].apply(bias)
        sigmas = df.groupby("zbin")["ez"].apply(sigma)
        outliers = df.groupby("zbin")["ez"].apply(outlier)
        cat_outliers = df.groupby("zbin")["ez"].apply(catastrophic_outlier)

    f = open(outfile, 'w')
    for z, b, s, o, c in zip(zvals, biases, sigmas, outliers, cat_outliers):
        f.write(str(z) + "  " + str(b) + "  " + str(s) + "  " + str(o) + "  " + str(c) + "\n")
    f.close()
    
    
#### Plotting

def histogram_ez(ez, ax, selection=None, plot_stats=False, **kwargs):
    """Histogram of ez=(zp-z)/(1+z)
       
       @param ez           Series of ez values
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param plot_stats   plot median, IQR, standard deviation, outliers
       @param **kwargs     matplotlib.pyplot.hist keyword arguments
       
       ez and selection must contain same indices
    """
    if (selection is not None):
        ax.hist(ez[selection], **kwargs)
    else:
        ax.hist(ez, **kwargs)
        
    if (plot_stats):
        ylims = ax.get_ylim()
        
        if (selection is not None):
            sig = np.std(ez[selection])
            med = np.median(ez[selection])
            x75, x25 = np.percentile(ez[selection],[75.,25.])
        else:
            sig = np.std(ez)
            med = np.median(ez)
            x75, x25 = np.percentile(ez,[75.,25.])
            
        eta = 0.15
        
        ax.plot([sig, sig], ylims, color='black')
        ax.plot([-sig, -sig], ylims, color='black')
        ax.plot([med, med], ylims, color='black')
        ax.plot([x25, x25], ylims, color='black', linestyle='dashed')
        ax.plot([x75, x75], ylims, color='black', linestyle='dashed')
        ax.plot([eta, eta], ylims, color='black', linestyle='dotted')
        ax.plot([-eta, -eta], ylims, color='black', linestyle='dotted')
        
        
def photoz_vs_z(photo_z, true_z, ax, selection=None, plot_sigs=False, guide_color='red', **kwargs):
    """Scatter plot of photo-z vs true-z with y=x line
       
       @param photo_z      Series containing point photo-z estimates
       @param true_z       Series containing true (or spec) redshifts 
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param plot_sigs    Plot lines of standard deviation
       @param guide_color  color of the guide lines (y=x and sigma if being plotted)
       @param **kwargs     matplotlib.pyplot.plot keyword arguments
       
       The selection is performed on both Series objects. All must contain same indices
    
    """
    # plot zp vs z, performing selection on the data if required
    # also calculate std in (zp-z)/(1+z) in case needed
    if (selection is not None):
        ax.plot(true_z[selection], photo_z[selection], **kwargs)
        sig = np.std((photo_z[selection]-true_z[selection])/(1.+true_z[selection]))
    else:
        ax.plot(true_z, photo_z, **kwargs)
        sig = np.std((photo_z-true_z)/(1.+true_z))
        
    # add y=x line and make sure axes start close to smallest photo-z
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    
    minv = np.min([xlims[0], ylims[0]])
    maxv = np.min([xlims[1], ylims[1]])
    ax.plot([minv, maxv], [minv, maxv], color=guide_color)
    
    # make sure x-axis starts close to where y-axis starts
    minv = np.min(photo_z)
    
    ax.set_xlim([minv, xlims[1]])
    ax.set_ylim([minv, ylims[1]])
    
    # add lines showing the std in (zp-z)/(1+z)
    if plot_sigs:
        xlims = ax.get_xlim()
        zvals = np.arange(xlims[0], xlims[1], 0.1)
        zplus = zvals + sig*(1+zvals)
        zminus = zvals - sig*(1+zvals)
        ax.plot(zvals, zplus, color=guide_color, linestyle='dashed')
        ax.plot(zvals, zminus, color=guide_color, linestyle='dashed')


def one_stat_vs_redshift(df, photo_z, true_z,  ax, stat='SIGMA', selection=None, zbins=[0,2,0.2], **kwargs):
    """Plot sigma, bias, OR outlier rate vs redshift
     
       @param df           dataframe containing photo-z data 
       @param photo_z      column name (str) containing photo-z estimates
       @param true_z       column name (str) containing true redshifts
       @param ax           axes handle to plot on
       @param stat         statistic to plot, one of: SIGMA, BIAS, OUTLIER, CATASTROPHIC
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
       @param **kwargs     matplotlib.pyplot.plot keyword arguments
    """
    if stat not in ['SIGMA', 'BIAS', 'OUTLIER', 'CATASTROPHIC']:
        raise ValueError("Error! unknown stat type: " + stat)
    
    
    # unpack redshift binning
    zmin = zbins[0]
    zmax = zbins[1]
    dbin = zbins[2]
    zvals = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # add column that indicates redshift bin "category" for each object
    df["zbin"] = pd.cut(df[true_z], np.arange(zmin, zmax+dbin, dbin))
    
    # groupby redshift bin and calculate stats
    if (selection is not None):
        if stat=='BIAS':
            statvals = df.loc[selection,:].groupby("zbin")["ez"].apply(bias)
        elif stat =='SIGMA':
            statvals = df.loc[selection,:].groupby("zbin")["ez"].apply(sigma)
        elif stat=='OUTLIER':
            statvals = df.loc[selection,:].groupby("zbin")["ez"].apply(outlier)
        elif stat=='CATASTROPHIC':
            statvals = df.loc[selection,:].groupby("zbin")["ez"].apply(catastrophic_outlier)
    else:
        if stat=='BIAS':
            statvals = df.groupby("zbin")["ez"].apply(bias)
        elif stat =='SIGMA':
            statvals = df.groupby("zbin")["ez"].apply(sigma)
        elif stat=='OUTLIER':
            statvals = df.groupby("zbin")["ez"].apply(outlier)
        elif stat=='CATASTROPHIC':
            statvals = df.groupby("zbin")["ez"].apply(catastrophic_outlier)
                          
    # get correct label, linestyle                              
    if stat=='BIAS':
        label = 'bias'
    elif stat =='SIGMA':
        label = '$\sigma$'
    elif stat=='OUTLIER':
        label = '$\eta$'
    elif stat=='CATASTROPHIC':
        label = '$\eta_{cat}$'

    ax.plot(zvals, statvals, label=label, **kwargs)


def all_stats_vs_redshift(df, photo_z, true_z,  ax, selection=None, zbins=[0,2,0.2], catastrophic=False, 
                          **kwargs):
    """Plot sigma, bias, AND outlier rate vs redshift
     
       @param df           dataframe containing photo-z data 
       @param photo_z      column name (str) containing photo-z estimates
       @param true_z       column name (str) containing true redshifts
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
       @param catastrophic outlier type, if true use catastrophic outlier definition
       @param **kwargs     matplotlib.pyplot.plot keyword arguments
    """
    
    # set outlier type
    stats = ['SIGMA', 'BIAS']
    if (catastrophic):
        stats.append('CATASTROPHIC')
    else:
        stats.append('OUTLIER')
        
    # set linestyles
    linestyles = ['solid', 'dashed', 'dotted']
        
    # plot each statistic in turn
    for stat, ls in zip(stats, linestyles):
        one_stat_vs_redshift(df, photo_z, true_z,  ax, stat, selection, linestyle=ls, **kwargs)

    
def add_LSST_req(ax, req_type='SIGMA', **kwargs):
    """Add LSST requirements to a plot
       
       @param ax           axes handle to plot on
       @param req_type     requirement type, one of: SIGMA, BIAS, OUTLIER, CATASTROPHIC
       @param **kwargs     matplotlib.pyplot.plot keyword arguments
    
    """
    xlims = ax.get_xlim()


    if req_type=='BIAS':
        req = lsst_bias_max
    elif req_type =='SIGMA':
        req = lsst_sigma_goal
    elif req_type=='OUTLIER' or req_type=='CATASTROPHIC':
        req = lsst_outlier_max
    
    ax.plot(xlims, [req, req], **kwargs)
    if req_type=='BIAS':
        ax.plot(xlims, [-req, -req], **kwargs)
        

def n_of_z(df, ez, true_z, ax, selection=None, zbins=[0,2,0.2], **kwargs):
    """Plot number of galaxies vs redshift
     
       @param df           dataframe containing photo-z data 
       @param ez           column name (str) containing ez = (zp-z)/(1+z)
       @param true_z       column name (str) containing true redshifts
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
       @param **kwargs     matplotlib.pyplot.plot keyword arguments
    """
    # unpack redshift binning
    zmin = zbins[0]
    zmax = zbins[1]
    dbin = zbins[2]
    zvals = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # add column that indicates redshift bin "category" for each object
    df["zbin"] = pd.cut(df[true_z], np.arange(zmin, zmax+dbin, dbin))
    
    # groupby redshift bin and count galaxies
    if (selection is not None):
        nz = df.loc[selection,:].groupby("zbin")["ez"].apply(lambda x: len(x))
    else:
        nz = df.groupby("zbin")["ez"].apply(lambda x: len(x))
        
    ax.plot(zvals, nz, **kwargs)
        
        
def distribution_vs_redshift(df, photo_z, true_z, ax, selection=None, zbins=[0,2,0.2], 
                             color='blue', **kwargs):
    """Plot distribution of photo-z in redshift bins
    
       @param df           dataframe containing photo-z data 
       @param photo_z      column name (str) containing photo-z estimates
       @param true_z       column name (str) containing true redshifts
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
       @param **kwargs     matplotlib.pyplot.hist keyword arguments
           
    """
    # unpack redshift binning
    zmin = zbins[0]
    zmax = zbins[1]
    dbin = zbins[2]
    zvals = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # add column that indicates redshift bin "category" for each object
    df["zbin"] = pd.cut(df[true_z], np.arange(zmin, zmax+dbin, dbin))
    

    # plot the boxes
    boxwidth = dbin/2.
    bp = df.boxplot(column=photo_z, by="zbin", return_type='dict', 
                     ax=ax, positions=zvals, widths=boxwidth, **kwargs)
    
    # set the tick labels along the redshift axes
    xtickNames = plt.setp(ax, xticklabels=zvals)
    
    # these boxes will be outlined by 'color'
    plot_color.color_outlines(bp, photo_z, color)
        
    # remove the auto-labelling
    ax.set_title("")
    ax.set_xlabel("")
    plt.suptitle("")
    
    return bp
    
    
def violin_vs_redshift(df, photo_z, true_z, ax, selection=None, zbins=[0,2,0.2], 
                             color='blue', **kwargs):
    """Plot distribution of photo-z in redshift bins
    
       @param df           dataframe containing photo-z data 
       @param photo_z      column name (str) containing photo-z estimates
       @param true_z       column name (str) containing true redshifts
       @param ax           axes handle to plot on
       @param selection    Series containing bool values
       @param zbins        list of: [z-min, z-max, dz-bin]
       @param **kwargs     matplotlib.pyplot.hist keyword arguments
           
    """
    # unpack redshift binning
    zmin = zbins[0]
    zmax = zbins[1]
    dbin = zbins[2]
    zvals = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # add column that indicates redshift bin "category" for each object
    df["zbin"] = pd.cut(df[true_z], np.arange(zmin, zmax+dbin, dbin))
    

    # plot the boxes
    boxwidth = dbin/2.
    bp = df.violinplot(column=photo_z, by="zbin", return_type='dict', 
                     ax=ax, positions=zvals, widths=boxwidth, **kwargs)
    
    # set the tick labels along the redshift axes
    xtickNames = plt.setp(ax, xticklabels=zvals)
    
    # these boxes will be outlined by 'color'
    plot_color.color_outlines(bp, photo_z, color)
        
    # remove the auto-labelling
    ax.set_title("")
    ax.set_xlabel("")
    plt.suptitle("")
    
    return bp
