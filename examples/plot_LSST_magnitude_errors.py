"""Example script that plots the LSST magnitude errors as a function of magnitude

"""
import photErrorModel as pem
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# list of LSST filters
filterList = ["LSSTu", "LSSTg", "LSSTr", "LSSTi", "LSSTz", "LSSTy"]
nFilter = len(filterList)


# initialise error model

# dictionary containing error model parameters, if not supplied defaults are assumed
# all values below ARE the defaults (so defining pars is redundant here, just for illustration)
pars = {} 
pars["tvis"] = 30.          # single visit exposure time (2x15s)
pars["sigmaSys"] = 0.0025   # expected irreducible error
pars["nYrObs"] = 1          # number of years of observations
# number of visits per year
pars["nVisYr"] = {'LSSTu':6,'LSSTg':8,'LSSTr':18,'LSSTi':18,'LSSTz':16,'LSSTy':16} 
# band dependent parameter
pars["gamma"] = {'LSSTu':0.037,'LSSTg':0.038,'LSSTr':0.039,'LSSTi':0.039,'LSSTz':0.040,'LSSTy':0.040}
# band dependent parameter
pars["Cm"] = {'LSSTu':23.60,'LSSTg':24.57,'LSSTr':24.57,'LSSTi':24.47,'LSSTz':24.19,'LSSTy':23.74}
# sky brightness
pars["msky"] = {'LSSTu':21.8,'LSSTg':22.0,'LSSTr':21.3,'LSSTi':20.0,'LSSTz':19.1,'LSSTy':17.5}
# seeing
pars["theta"] = {'LSSTu':0.77,'LSSTg':0.73,'LSSTr':0.70,'LSSTi':0.67,'LSSTz':0.65,'LSSTy':0.63}
# extinction coefficient
pars["km"]= {'LSSTu':0.48,'LSSTg':0.21,'LSSTr':0.10,'LSSTi':0.07,'LSSTz':0.06,'LSSTy':0.06}
pars["airMass"] = 1.2       # air mass
pars["extendedSource"] = 0. # extended source model: simple constant added to m5 (makes m5 fainter)
pars["minFlux"] = 2.5e-40   # set minimum allowed flux (equivalent to mag=99)

ep = pem.LSSTErrorModel(pars)

# custom printing of error model parameters
print ep

# magnitude range to plot over
minMag = 14.
maxMag = 40.
nMag = 100
dMag = (maxMag - minMag)/(nMag - 1.)


# store errors at each magnitude and in each filter
magErrors = np.zeros([nMag, nFilter+1])
for i in xrange(nMag):
    
    mag = minMag + i*dMag
    magErrors[i,0] = mag
    
    for j in xrange(nFilter):
    
        err = ep.getMagError(mag, filterList[j])
        magErrors[i,j+1] = err
        
fig = plt.figure(figsize=(10,10))
            
for j in xrange(nFilter):
    
    ax = fig.add_subplot(2,3,j+1)
    plt.plot(magErrors[:,0], magErrors[:,j+1],  linestyle='solid',color='black')
    ax.set_xlabel('magnitude', fontsize=24)
    ax.set_ylabel('$\sigma_m$', fontsize=24)
    ax.set_title(filterList[j])
    ax.set_yscale('log')
    
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

plt.show(block=True)
