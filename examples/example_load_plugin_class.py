"""Example of how to load particular photometric error model

"""
import photErrorModel as pem
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# list of LSST filters
filterList = ["LSST_u", "LSST_g", "LSST_r", "LSST_i", "LSST_z", "LSST_y"]
nFilter = len(filterList)


# Select the class you want to load by supplying it's name as a string to getattr
# However this doesn't create an instance
ep = getattr(pem, "LSSTErrorModel")
print ep

minMag = 14.
maxMag = 40.
nMag = 100
dMag = (maxMag - minMag)/(nMag - 1.)

magErrors = np.zeros([nMag, nFilter+1])
for i in xrange(nMag):
    
    mag = minMag + i*dMag
    magErrors[i,0] = mag
    
    for j in xrange(nFilter):
    
        # need to call the initializer at the same time as the method
        err = ep().getMagError(mag, filterList[j])
        magErrors[i,j+1] = err
        
fig = plt.figure(figsize=(10,10))
            
for j in xrange(nFilter):
    
    ax = fig.add_subplot(3,3,j+1)
    plt.plot(magErrors[:,0], magErrors[:,j+1],  linestyle='solid',color='black')
    ax.set_xlabel('magnitude', fontsize=24)
    ax.set_ylabel('$\sigma_m$', fontsize=24)
    ax.set_title(filterList[j])
    ax.set_yscale('log')
    
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

plt.show(block=True)
