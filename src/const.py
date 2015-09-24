class _const:
    class ConstError(TypeError): pass
    def __setattr__(self,name,value):
        if self.__dict__.has_key(name):
            raise self.ConstError, "Can't rebind const(%s)"%name
        self.__dict__[name]=value
import sys
sys.modules[__name__]=_const()

# Ok now set all constants wanted  below here:
import const
import cmath

# maths constants
const.pi = 3.141592653589793
const.i = cmath.sqrt(-1.)
const.e = 2.71828182845904523536028747135266249775724709369995

# physics constants
const.cKms = 2.9979e5  # speed of light km/s
const.cms = 2.9979e8   # speed of light m/s
const.GmKgs = 6.67e-11 # gravitational constant m^3/Kg/s^2

# astro conversions
const.Mpc2m = 3.0856e22 # m
const.m2pc = 1./3.08567758e16 # pc
const.m2Mpc = 1./const.Mpc2m  # Mpc
const.deg2rad = const.pi/180. # degrees into radians
const.rad2deg = 180./const.pi # radians into degrees
const.arcsec2rad = 3600.*const.deg2rad # arcsecs into radians
const.arcsec2deg = 1./3600 # arcsecs into degrees
const.arcsec2kpc = 1e3*const.pi/(180.*3600.) # arcsecs into kpc
const.arcsec2Mpc = 1e6*const.pi/(180.*3600.) # arcsecs into Mpc
const.deg2arcsec = 3600. # degrees into arcsecs

# astro constants
const.msolarKg = 1.989e30      # solar mass in Kg
const.skySteRad = 4.*const.pi  # area of full sky in steradians
const.skySqDeg = const.skySteRad*(180./const.pi)*(180./const.pi) # area of full sky in sq. degrees=41252.9612

# cosmology constants
const.H100hs = 100.*1e3/const.Mpc2m # Hubble constant in h/s
const.T_CMB_K = 2.725  # temperature of CMB in K
const.T_NU_K = 1.9     # temperature of neutrinos in K
## critical density
const.rho_crit = 3.*(const.H100hs*const.H100hs)/(8.*const.pi*const.GmKgs)/ \
       const.msolarKg*(const.Mpc2m*const.Mpc2m*const.Mpc2m) # h^2 M_solar/Mpc^-3
const.RHOCRIT_h2KGM3 = 1.879e-26                # in h^2 Kg/m^3
const.RHOCRIT_h2MSOLMPC3 = 2.77536627e11        # in MSol/Mpc^3
const.RHOCRIT_h2GEVCM3 = 1.0539e-5              # in GeV/cm^3
const.FREQ_21CM_HI_IN_GHZ = 1.420405751786  # frequency of 21cm radiation in GHz
const.PHOTON_DENSITY_TODAY_KGM3 = 4.6417e-31    # Photon density today (kg/m^3) 
const.PHOTON_NUMDENSITY_PER_CM3_TODAY = 410.4   # Photon number density today /cm^3 
const.PHOTON_AND_NU_DENSITY_TODAY_KGM3 = 7.8042e-31 # Photon+neutrino density (kg/m^3) 
const.Z_RECOMBINATION = 1088.                   # Value in WMAP ApJ paper
const.DELTA_C = 1.686                           # Overdensity for collapse (spherical model)


