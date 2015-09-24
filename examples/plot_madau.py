"""Example script to plot Madau IGM absorption model

   Compare plot to Figure 3 in Madau 1995

"""
import igmModel
import numpy as np
import matplotlib.pyplot as plt

d = {}
d["lineNumStop"] = 5  # max number of Lyman lines to account for
d["doLyCont"] = True  # turn on/off Lyman continuum contribution
madauIGM = igmModel.MadauIGM(d)

# custom instance printing
print madauIGM

# source galaxy redshifts to compute IGM stransmission for 
z1 = 2.5
z2 = 3.5
z3 = 4.5

# observed wavelengths to compute end IGM transission at
# angstroms
lamMin = 2000.
lamMax = 8000.
nLam = 1000
dLam = (lamMax - lamMin)/(nLam - 1.)

wavelengths = []
madau = np.zeros([nLam, 3])
for i in range(nLam):
    
    lam = lamMin + i*dLam
    wavelengths.append(lam)
    
    lam *= 1e-10 # covert angstroms to meters
    madau[i,0] = madauIGM.getObserverFrameTransmission(lam, z1)
    madau[i,1] = madauIGM.getObserverFrameTransmission(lam, z2)
    madau[i,2] = madauIGM.getObserverFrameTransmission(lam, z3)



fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)
ax.plot(wavelengths, madau[:,0], color='black', label=z1)
ax.plot(wavelengths, madau[:,1], color='black', label=z2)
ax.plot(wavelengths, madau[:,2], color='black', label=z3)
ax.set_xlabel('wavelength (angstroms)', fontsize=24)
ax.set_ylabel('transmission', fontsize=24)

plt.show(block=True)
