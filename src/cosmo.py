"""

 To contain class that performs standard cosmological calculations:
 
 - Hubble parameter
 - luminosity distance
 - angular diameter distance
 - comoving distance
 - etc

"""

import math

import const

## Class for computing cosmological quantities
#
# This class is a python port from the C++ version by R. Ansari, C. Magneville
#
# important getter methods, call after setting redshift with setEmissionRedShift(z):
#
# DH()                          returns c/H0 in Mpc (Mpc/h if h=1)
# AngularDiameterDistance()     returns unitless angular diameter distance
# AngularDiameterDistanceMpc()  returns angular diameter distance in Mpc (Mpc/h if h=1)
# LuminosityDistance()          returns unitless luminosity distance
# LuminosityDistanceMpc()       returns luminosity distance in Mpc (Mpc/h if h=1)
# LineOfSightComovDistance()    returns unitless line of sight comoving distance 
# LineOfSightComovDistanceMpc() returns line of sight comoving in Mpc (Mpc/h if h=1)
# TransComovDistance()          returns unitless transverse comoving distance 
# TransComovDistanceMpc()       returns transverse comoving distance in Mpc (Mpc/h if h=1)
#
# important setter methods, though really should set through instantiation
# sethubble(h)
# setOmegaMatter(om)
# setOmegaLambda(ol)
# setDarkEnergy(ol, w0, wa)
#
# Last modified: 29 July 2014 by AA
class cosmologyCalculator(object):


    ## The constructor
    def __init__(self, h=1., omegamat=0.3, omegaDE=0.7, wX=-1., wXa=0.):
    
    
        # cosmological parameters
        self.h = h               # H0 in units of 100 km/s/Mpc : 0.5 -> 50 km/s/MPc  
        self.omegamat = omegamat # Total matter density in universe @z=0 (inc baryons) 
        self.omegaDE = omegaDE   # cosmological constant or dark energy density        
        self.wX = wX             # Dark energy density equation of state constant par  
        self.wXa = wXa           # Dark energy density equation of state evolving par  
        
        # can be changed using methods ->
        ## Baryon density (included in omegamat_)
        self.omegabaryon = 0.0441 
        # Total radiation density at z=0 
        self.omegarad = const.PHOTON_AND_NU_DENSITY_TODAY_KGM3 /               \
                        (self.h*self.h*const.RHOCRIT_h2KGM3)
        # Photon energy density, (included in omegarad_)   
        self.omegaphot = const.PHOTON_DENSITY_TODAY_KGM3 /                     \
                        (self.h*self.h*const.RHOCRIT_h2KGM3)   
        
        # True-> Non zero cosmological constant/dark energy
        if (self.omegaDE>1e-79):
            self.hasDE = True
        else:
            self.hasDE = False
            
        # True-> Eqn. of state of DE component is cosmo. constant 
        if ( (math.fabs(self.wX+1.) > 0.) or (math.fabs(self.wXa) > 0.) ):
            self.isLambda = False
        else: 
            self.isLambda = True 
            
        # Set omegaCurv, kcurvature if needed
        self.computeCurvature()        
                   
        # We initilalize self.ze, ...  to today
        self.initialize()
        #print 'in init',self.omegaCurv,self.kcurvature
        
        
    ## Redefine __repr__ and __str__ to be more useful
    
    def __repr__(self):
        print_str = 'cosmologyCalculator( h = {0.h!r}, \n'
        print_str +='                    Om = {0.omegamat!r}, Ol = {0.omegaDE!r}, Ob = {0.omegabaryon!r}, \n'
        print_str +='                    Or = {0.omegarad!r}, Op = {0.omegaphot!r}, \n'
        print_str +='                    w0 = {0.wX!r}, wa = {0.wXa!r}, \n'
        print_str +='                    Ok = {0.omegaCurv!r}, K = {0.kcurvature!r}, \n'
        print_str +='                    current redshift = {0.ze!r}, \n'
        print_str +='                    l.o.s dc = {0.dcl!r}, trans. dc = {0.dct!r}, chi = {0.chie!r}, \n'
        print_str +='                    dL = {0.dl!r}, dA ={0.da!r} )'
        
        return print_str.format(self)
        
        
    def __str__(self):
        print_str = 'cosmologyCalculator( h = {0.h!s}, \n'
        print_str +='                    Om = {0.omegamat!s}, Ol = {0.omegaDE!s}, Ob = {0.omegabaryon!s}, \n'
        print_str +='                    Or = {0.omegarad!s}, Op = {0.omegaphot!s}, \n'
        print_str +='                    w0 = {0.wX!s}, wa = {0.wXa!s}, \n'
        print_str +='                    Ok = {0.omegaCurv!s}, K = {0.kcurvature!s}, \n'
        print_str +='                    current redshift = {0.ze!s}, \n'
        print_str +='                    l.o.s dc = {0.dcl!s}, trans. dc = {0.dct!s}, chi = {0.chie!s}, \n'
        print_str +='                    dL = {0.dl!s}, dA ={0.da!s} )'
        return print_str.format(self)
        
        
    ## Compute curvature component of universe from other components' values
    # @param self    The object pointer
    def computeCurvature(self):
    
        self.omegaCurv = 1. - (self.omegamat + self.omegarad + self.omegaDE)
    
        if (math.fabs(self.omegaCurv) < 1.e-39):
            self.omegaCurv = 0.
            self.kcurvature = 0
        else:
            if (self.omegaCurv < 0):
                self.kcurvature = 1  
            else:
                self.kcurvature = -1
        
        self.initialize()
        #print 'in computecurvature',self.omegaCurv,self.kcurvature
    
    
    ## Initialize "current" values of various quantities to zero
    # @param self    The object pointer
    def initialize(self):

        # We initilalize self.ze, ...  to today
        self.ze = 0.       # Emission redshift
        self.te = 0.       # Emission time (0=today)  
        self.chie = 0.     # Chi coordinate 
        self.da = 0.       # Angular diameter distance (units c/H0)   
        self.dl = 0.       # Luminosity distance (units c/H0) 
        self.dcl = 0.      # Line of sight comoving distance (units c/H0) 
        self.dct = 0.      # Transverse comoving distance (units c/H0) 
        self.EofZe = 1.    # E(ze)    
        self.integGz = 0   # Integral[G(z) dz]_[ze_ ... 0] 
        self.integGTz = 0  # Integral[GT(z) dz]_[ze_ ... 0]
        self.dVc = 0.      # Comoving volume element


    ## Set matter density
    # @param self        The object pointer
    # @param omegamat    matter density
    def setOmegaMatter(self,omegamat):

        self.omegamat = omegamat
        if (self.omegamat < self.omegabaryon):
            print " cosmologyCalculator.setOmegaMatter()/Warning  Setting \
                    OmegaBaryon=OmegaMatter=", omegamat 
            self.omegabaryon = self.omegamat
        
        self.computeCurvature()



    ## Set dark energy density
    # @param self       The object pointer
    # @param omegaDE    dark energy density
    def setOmegaLambda(self, omegaDE):

        self.omegaDE = omegaDE
        if (self.omegaDE > 1.e-39):
            self.hasDE = True
        else:
            self.hasDE = False

        self.computeCurvature()



    ## Set radiation density
    # @param self       The object pointer
    # @param omegarad   radiation density
    def setOmegaRadiation(self, omegarad):

        self.omegarad = omegarad
        if (self.omegarad < self.omegaphot):
            print " cosmologyCalculator.setOmegaRadiation()/Warning  Setting \
                    OmegaPhoton=OmegaRad=", omegarad 
            self.omegaphot = self.omegarad
        
        self.computeCurvature()



    ## Set baryon density
    # @param self       The object pointer
    # @param omegab     baryon density
    def setOmegaBaryon(self, omegabar):

        self.omegabaryon = omegabar
        if (self.omegabaryon > self.omegamat):
            print " cosmologyCalculator.setOmegaBaryon()/Warning  Setting \
                    OmegaMatter=OmegaBaryon=", omegabar
            self.omegamat = self.omegabaryon
            
        self.computeCurvature()
         
         
    ## Set photon density
    # @param self       The object pointer
    # @param omegaphot  Photon density
    def setOmegaPhoton(self, omegaphot):

        self.omegaphot = omegaphot
        if (self.omegaphot > self.omegarad):
            print " cosmologyCalculator.setOmegaPhoton()/Warning  Setting \
                    OmegaRad=OmegaPhoton=", omegaphot 
            self.omegarad = self.omegaphot
            
        self.computeCurvature()
        
  
    ## Set dark enerqy
    # @param self       The object pointer
    # @param omegaDE    Dark energy density
    # @param wX         dark energy equation of state parameter
    # @param wXa        dark energy equation of state parameter
    def setDarkEnergy(self, omegaDE, wX, wXa):


        self.omegaDE = omegaDE
        self.wX = wX
        self.wXa= wXa
        if (self.omegaDE > 1.e-39):
            self.hasDE = True
        else: 
            self.hasDE = False
        
        if ( (math.fabs(self.wX+1.) > 0.) or (math.fabs(self.wXa) > 0.) ):
            self.isLambda = False
        else:
            self.isLambda = True
    
        self.computeCurvature()



    ## Set dark enerqy equation of state
    # @param self       The object pointer
    # @param wX         dark energy equation of state parameter
    # @param wXa        dark energy equation of state parameter
    def setDarkEnergyEoS(self, wX, wXa=0.):

        # First check there is some dark energy
        if ( self.omegaDE < 1.e-9 ):
            emsg = "setDarkEnergyEoS warning: Negligible amount of dark energy"
            print emsg

        self.wX = wX
        self.wXa= wXa
        if ( (math.fabs(self.wX+1.) > 0.) or (math.fabs(self.wXa) > 0.) ):
            self.isLambda = False
        else:
            self.isLambda = True



    ## Set flat universe using matter density
    # @param self       The object pointer
    def setFlatUniverse_OmegaMatter(self):

        self.omegamat = 1. - (self.omegarad + self.omegaDE)
  
        if (self.omegamat < 0.): 
            print " cosmologyCalculator.setFlatUniverse_OmegaMatter()/Error ->  \
                    OmegaMatter=",  self.omegamat , " < 0: SETTING TO ZERO" 
            self.omegamat = 0.
        
        self.computeCurvature()



    ## Set flat universe using dark energy density
    # @param self       The object pointer
    def setFlatUniverse_OmegaLambda(self):

        self.omegaDE = 1. - (self.omegamat + self.omegarad)
        if (self.omegaDE < 0.):
            print " cosmologyCalculator.setFlatUniverse_OmegaLambda()/Error ->  \
                    OmegaL=",  self.omegaDE , " < 0: SETTING TO ZERO" 
            self.omegaDE = 0.

        if (self.omegaDE > 0.):
            self.hasDE = True
        else:
            self.hasDE = False
        
        self.computeCurvature()
       
    ## Set hubble parameter
    # @param self      The object pointer
    # @param h    
    def sethubble(self,h):
        self.h = h
    
    def H0(self):
        return 100*self.h
         
    def DH(self):
        return const.cKms/self.H0()    
            
    def AngularDiameterDistance(self):
        return self.da
        
    def AngularDiameterDistanceMpc(self):
        return self.da*self.DH()
        
    def LuminosityDistance(self):
        return self.dl
        
    def LuminosityDistanceMpc(self):
        return self.dl*self.DH()
        
    def LineOfSightComovDistance(self):
        return self.dcl
        
    def LineOfSightComovDistanceMpc(self):
        return self.dcl*self.DH()
        
    def TransComovDistance(self):
        return self.dct
       
    def TransComovDistanceMpc(self):
        return self.dct*self.DH()
        
    def comovingVolumeElement(self):
        return self.dVc

        
            
    ## Return true if universe is flat
    # @param self       The object pointer
    def isFlat(self):
        return self.kcurvature == 0
        
    
    ## Return true if universe is open
    # @param self       The object pointer
    def isOpen(self):
        return self.kcurvature == -1
        
        
    ## Return true if universe is closed
    # @param self       The object pointer
    def isClosed(self):
        return self.kcurvature == 1


    ## Return E(z) : Notation of Peebles / Principle of Physical Cosmology 
    # (Chap 13)
    # @param self   The object pointer                                                  
    # @param z      redshift  
    def Ez(self,z):

        zp1 = 1. + z
        ez2 = zp1*zp1*( self.omegarad*zp1*zp1 + self.omegamat*zp1 + self.omegaCurv)
        if (self.hasDE and self.isLambda):
            ez2 += self.omegaDE
        elif (self.hasDE and not self.isLambda):
            ez2 += self.omegaDE*pow(zp1,3.*(1.+self.wX+self.wXa))* \
                   math.exp(-3*self.wXa*(1.-1./zp1))
  
        return math.sqrt(ez2)
        
        
            
    ## Perform calculations at redshift
    # @param self       The object pointer
    # @param ze         redshift of photon emission
    # @param prec       precision of calculation
    # @param fginc      do incremental calculation
    def setEmissionRedShift(self, ze, prec=0.001, fginc=True):

  
        if (ze < 0.): # z is unphysical
            print " cosmologyCalculator::setEmissionRedShift/Error ze = ", ze , \
                  " less than 0 ! SETTING TO ZERO"
            ze = 0.


        if (ze < 1.e-39): # if ~today
            self.ze = 0.
            self.te = 0.
            self.chie = 0.
            self.da = 0.
            self.dl = 0.
            self.dcl = 0.
            self.dct = 0.
            self.EofZe = 1.
            self.integGz = 0
            self.integGTz = 0
            return
        
        zelast = self.ze      # set previous emission redshift to be last 
        self.ze = ze          # set current emission redshift
        self.EofZe = self.Ez(ze)

        # Reference: Principles of Physical Cosmology - P.J.E. Peebles  
        #            Princeton University Press - 1993
        #              ( See Chapter 13)
        # We have to integrate Integral(dz / E(z)) from 0 to ze  (cf 13.29)
        #      E(z) = Sqrt(Omega0*(1+z)^3 + OmegaCurv*(1+z)^2 + OmegaL) (cf 13.3)
        # G(z) = 1/E(z) : integration element for the calculation of D_A  
        # GT(z) = 1/[(1+z)*E(z)] : integration element for the calculation of 
        #                          age/time (dt = GTz(z)*dz)

    
        idzoez = 0.
        idzToez = 0.
        if (fginc):  # Incremental calculation
            inGz = 0.
            inGTz = 0.
        
            if (zelast < self.ze):
                # with precision prec:
                # Numerical integration of i) inGz: G(z)dz = 1/E(z) dz 
                #                             between [zelast,ze]
                # Numerical integration of ii) inGTz: GT(z)dz = 1/[(1+z)*E(z)] dz 
                #                              between [zelast,ze]
                inGz, inGTz = self.NumIntegrateGzGTz(zelast, self.ze, prec)
                self.integGz += inGz
                self.integGTz += inGTz
            
            else:
                inGz, inGTz = self.NumIntegrateGzGTz(self.ze, zelast, prec)
                self.integGz -= inGz
                self.integGTz -= inGTz
            
        else: # not incremental
             self.integGz, self.integGTz = self.NumIntegrateGzGTz(0., self.ze, prec) 
  
  
        idzoez = self.integGz
        idzToez = self.integGTz
    
  
        self.dcl = idzoez # line of sight comoving distance (units c/H0)
        # Notation of Hogg 1999:
        # Hubble distance:
        # dh = c/H0
        # Line of sight comoving distance:
        # dcl = dh int_0^z dz/E(z')
        # Transverse comoving distance:
        # i)     OPEN UNI dm = dh 1/sqrt(OmK)sinh[sqrt(OmK) dc/dh]
        # ii)    FLAT UNI dm = dc
        # iii) CLOSED UNI dm = dh 1/sqrt(OmK)sin[sqrt(|OmK|) dc/dh]
    
        if ( self.isFlat() ):
            # here self.chie= int_z1^z2 dz/E(z')
            # therefore it is equal to dc/dh=dm/dh
            self.dct = self.chie = self.da = idzoez  
        else:
            if (self.isOpen()):
                sov = math.sqrt(self.omegaCurv)
            else:
                sov = math.sqrt(-self.omegaCurv)

            # here self.chie= sqrt(OmK) int_z1^z2 dz/E(z') = sqrt(OmK) dc/dh
            # **and is NOT equal to dc/dh or dm/dh**
            self.chie = sov * idzoez
            if (self.isOpen()):
                self.da = math.sinh(self.chie)
            else: 
                self.da = math.sin(self.chie)
            self.da /= sov
            self.dct = self.da
        
        self.da /= (1+self.ze)

        # There is a factor (1+self.ze)^2 between self.da and self.dl
        # > (1+self.ze) to get self.da to the value z=0 (today) (1)
        # > (1+self.ze) to compensate for the redshift of the photons (2)
        # See for example Weinberg, Gravitation and Cosmology, page 423 (14.4.22)
        self.dl = self.da*(1+self.ze)*(1+self.ze)   

        # The look back time
        self.te = idzToez
        #print ze,self.da,self.dl,self.dct

        # Comoving volume element in Mpc^3 (or (Mpc/h)^3 if h=1)
        daMpc = self.AngularDiameterDistanceMpc()
        self.dVc = self.DH()*(1+self.ze)*(1+self.ze)*daMpc*daMpc/self.EofZe
        


    ## Numerical integration of G(z) dz = 1/E(z) dz and numerical integration 
    # of GT(z) dz = 1/[(1+z)*E(z)] dz between [z1,z2]
    # @param self   The object pointer 
    # @param z1     Lower redshift bound
    # @param z2     Upper redshift bound
    # @param prec   The target precision used to compute the adaptative 
    #               integration step
    # returns resG,resGT
    # @param resG   Resulting integration of G(z)
    # @param resGT  Resulting integration of GT(z)                               
    def NumIntegrateGzGTz(self, z1, z2, prec):
        resG, resGT = self.NumIntegGzGTz(z1, z2, prec);
        return resG, resGT
        
        
    ## Numerical integration of G(z) dz = 1/E(z) dz and numerical integration 
    # of GT(z) dz = 1/[(1+z)*E(z)] dz between [z1,z2]
    # @param self   The object pointer 
    # @param z1     Lower redshift bound
    # @param z2     Upper redshift bound
    # @param prec   The target precision used to compute the adaptative 
    #               integration step
    # returns resG,resGT
    # @param resG   Resulting integration of G(z)
    # @param resGT  Resulting integration of GT(z)  
    def NumIntegGzGTz(self, z1, z2, prec):
   
        if ((z1 < 0.) or (z2 < 0.) or (z2 < z1)):
            print " SimpleUniverse::NumIntegGzGTz()/Error invalid values for z1,z2: ", \
                  z1 , "," , z2 ,
        
        # How the precision is used to calculate the step - 
        # We use the relation :
        #  | Integ_{a,b} f(x) dx - (b-a)f(a) | <= (b-a)^2/2 Sup_{a,b} |f'(x)| 
        # The function to integrate G(z) is a decreasing positive function
        # We approximate Sup |f'(x)| par |f'(a)| = -f'(a) 
        # I = Integ_{a,b} f(x) dx  >= (b-a)f(b)    , delta I / I <= prec
        #  ===> delta I/I <= (b-a)/2 | f'(a) / f(b) |  < prec 
        # step = (b-a) ===> b-a <= 2*prec* | f(b)/f'(a) | 
        # We will use : step = prec* | f(b)/f'(a) |
   
        resG = resGT = 0.
        
        z1orig = z1
        z2orig = z2

        # Integration with a fixed step up to z=zbrk (=10.)
        zbrk = 10.
        # for z < zbrk, we take a=0 , G'(0) ~ 1 G(1) ~ 1 ===> step~prec 
        cfac = math.fabs(self.DerivateGz(0.,1.e-3));
        if (cfac < 1.): 
            cfac = 1.
        
        dz = prec*self.Gz(1)/cfac;
        if (dz > 0.1):
            dz = 0.1
        dzini = dz
    
        #DBG  cout << " **DBG*A**  prec=" << prec << " zbrk=" << zbrk
        #DBG << " G(1.)= " << Gz(1.) << " dG(0)=" << DerivateGz(0.,dz) 
        #    << " ==>dz " << dz << endl;
   
        zc = z1
        accG = accGT = 0.
        if (z1<zbrk): # Integration up to zbrk 
            
            if (z2 > zbrk):
                z2 = zbrk
     
            while (zc < z2): # Calculate time and distance
                accG += self.Gz(zc)
                accGT += self.GTz(zc)
                zc += dz;   
            
            zc -= dz
            accG -= 0.5*(self.Gz(z1)+self.Gz(zc))
            accGT -= 0.5*(self.GTz(z1)+self.GTz(zc))
            accG *= dz 
            accGT *= dz
            dz = z2-zc
            
            if (dz > 1.e-69): # Correction for the last step 
                accG += 0.5*dz*(self.Gz(zc)+self.Gz(z2))
                accGT += 0.5*dz*(self.GTz(zc)+self.GTz(z2))
        
   
        resG = accG
        resGT = accGT

        if (z2orig <= zbrk):
            return resG,resGT
   
        z1 = zbrk
        # if integration after z1 > zbrk 
        if (z1orig > zbrk):
            z1 = z1orig
            dz = prec*math.fabs(Gz(z1*4.)/self.DerivateGz(z1,0.1))
            #DBG    cout << " **DBG*C**  dz=" << dz << endl; 
        else:
            dz = dzini

        z2 = z2orig;


        kkMAX = 100
        fgencore = True
        MAXRAPDZ = 10.
   
        zc = z1;
        while ( (zc < z2) and fgencore ):

            Delz = float(kkMAX)*dz
            if (Delz < 1):
                Delz = 1.

            # Calculation of step :
            lastdz = dz
            dz = prec*math.fabs(self.Gz(zc+Delz)/self.DerivateGz(zc,dz))
            # dzs = dz;

            if (dz/lastdz > MAXRAPDZ):
                dz = lastdz*MAXRAPDZ
                
            kkmx = long(((z2-zc)/dz)-1)
            if (kkmx > kkMAX):
                kkmx = kkMAX

        #DBG cout << " **DBG*B**  prec=" << prec << " zc=" << zc << " Delz=" << Delz 
        #DBG  << " G(zc)= " << Gz(zc+Delz) << " dG=" << DerivateGz(zc,dz) << " ==>dz " 
        #DBG  << dz << " dzs=" << dzs << " kkmx=" << kkmx << endl; 

            za = zc
            accG = accGT = 0.
            for kkp in range(kkmx):
                accG += self.Gz(zc)
                accGT += self.GTz(zc)
                zc += dz
            
            zc -= dz

            accG -= 0.5*(self.Gz(za)+self.Gz(zc))
            accGT -= 0.5*(self.GTz(za)+self.GTz(zc))
            accG *= dz
            accGT *= dz
            resG += accG
            resGT += accGT
             
            if (kkmx < kkMAX):
                fgencore = False

   
        dz = z2-zc
        if (dz > 1.e-69): # Correction for the last step 
            resG += 0.5*dz*(self.Gz(zc)+self.Gz(z2))
            resGT += 0.5*dz*(self.GTz(zc)+self.GTz(z2))

        return resG,resGT
        
    
        

    ## Value of d G(z) / dz (no dark energy)
    # G(z) = ( sum_components )^(-0.5)
    # dG/dz = -1/2*(sum_components)^(3/2) * d(sum_components)/dz
    # sum_components is equal to E^2
    # @param self   The object pointer 
    # @param z    redshift                                                  
    def DerivateGz_NoLambda(self, z):
    
        return ( -(2.*self.omegaCurv*(1+z) + 3.*self.omegamat*(1+z)*(1+z)  \
                    + 4.*self.omegarad*(1+z)*(1+z)*(1+z)) /                \
                  (2.*pow(self.omegaCurv*(1+z)*(1+z)                       \
                    + self.omegamat*(1+z)*(1+z)*(1+z)                      \
                    + self.omegarad*(1+z)*(1+z)*(1+z)*(1+z),  1.5)) ) 


    ## Value of d G(z) / dz (non-zero dark energy but is cosmological constant)
    # @param self   The object pointer 
    # @param z    redshift      
    def DerivateGz_Lambda(self, z):

        return ( -(2.*self.omegaCurv*(1+z) + 3.*self.omegamat*(1+z)*(1+z)      \
                    + 4.*self.omegarad*(1+z)*(1+z)*(1+z)) /                    \
                  (2.*pow(self.omegaDE + self.omegaCurv*(1+z)*(1+z)            \
                    + self.omegamat*(1+z)*(1+z)*(1+z)                          \
                    + self.omegarad*(1+z)*(1+z)*(1+z)*(1+z),  1.5)) )



    ## Value of d G(z) / dz (non-zero dark energy but is NOT cosmological const)
    # @param self   The object pointer 
    # @param z    redshift  
    def DerivateGz_Lambda_X(self, z):


        Esq =  self.omegaCurv*(1.+z)*(1.+z) 
        Esq += self.omegamat*(1.+z)*(1.+z)*(1.+z) 
        Esq += self.omegarad*(1.+z)*(1.+z)*(1.+z)*(1.+z) 
        # + omegaX_*pow(1+z, 3.*(1.+self.wX)),1.5) ) 
        Esq += self.omegaDE*pow(1.+z,3.*(1.+self.wX+self.wXa))                 \
               *math.exp(-3*self.wXa*(1-1./(1.+z)))

        #print self.omegaCurv,self.omegamat,self.omegarad,self.omegaDE,
        #print self.wX ,self.wXa,": ",Esq

        part1 = -( 2*self.omegaCurv*(1.+z) + 3*self.omegamat*(1.+z)*(1.+z)     \
                  + 4*self.omegarad*(1.+z)*(1.+z)*(1.+z)                       \
                 + 3.*(1.+self.wX+self.wXa)*self.omegaDE ) 

        part2 = math.exp(-3*self.wXa*(1-1./(1.+z)))                            \
                 *(   pow(1.+z,3.*(2./3.+self.wX+self.wXa))                    \
                    - 3*self.wXa*pow(1.+z,1.+3.*(self.wX+self.wXa)) )          \
                    / (2.*pow(Esq,1.5))
        
        #print Esq,part1,part2
        return part1*part2

        """return ( -(    2*self.omegaCurv*(1.+z) + 3*self.omegamat*(1.+z)*(1.+z) \
                 + 4*self.omegarad*(1.+z)*(1.+z)*(1.+z)                        \
        #+ 3*omegaX_*(1+self.wX)*pow(1.+z, -1.+3.*(1.+self.wX)) ) 
                 + 3.*(1.+self.wX+self.wXa)*self.omegaDE   )                   \
                 * math.exp(-3*self.wXa*(1-1./(1.+z)))                         \
                 *(   pow(1.+z,3.*(2./3.+self.wX+self.wXa))                    \
                    - 3*self.wXa*pow(1.+z,1.+3.*(self.wX+self.wXa)) )          \
                    / (2.*pow(Esq,1.5)) )"""
                    
                      

    ## Return G(z) = 1/E(z) : integration element for the calculation of D_A
    # @param self   The object pointer
    # @param z    redshift                                                  
    def Gz(self, z):
        return ( 1./self.Ez(z) )
        
    
    ## Return GT(z) = 1/[(1+z)*E(z)] : integration element for the calculation 
    # of age/time (dt = GTz(z)*dz)
    # @param self   The object pointer 
    # @param z    redshift                                                  
    def GTz(self, z):
        return ( 1./((1+z)*self.Ez(z)) )
    

    ## Value of d G(z) / dz
    # @param self   The object pointer 
    # @param z    redshift                                                  
    # @param dz   derivation step?                                          
    def DerivateGz(self, z, dz):

        if (dz > 1.e-19): # Calculate a numerical derivative
            return ((self.Gz(z+dz)-self.Gz(z))/dz)
            
        else: # Using the analytic form, problem ??
            if (self.hasDE and self.isLambda):
                return self.DerivateGz_Lambda(z)
            elif (self.hasDE and not self.isLambda):
                return self.DerivateGz_Lambda_X(z)
            else:
                return self.DerivateGz_NoLambda(z)
                
        
    ## Print parameter values
    # @param self   The object pointer        
    def print_pars(self): 

        print " OmegaMatter = ", self.omegamat,"OmegaRad = ", self.omegarad
        print " OmegaBaryon = ", self.omegabaryon," OmegaCDM= ",self.omegamat-self.omegabaryon
        print " OmegaPhoton = ", self.omegaphot
        

        if (self.hasDE):
            print " Non-zero dark energy: OmegaDE = ", self.omegaDE
            if (self.isLambda):
                print " Dark energy is cosmlogical constant:",
            else:
                print " Dark energy is quintessence:",
            print " w0 = ",self.wX,'wa = ',self.wXa
        else:
            print " No dark energy: OmegaDE = ", self.omegaDE
            
        print " OmegaCurv = ", self.omegaCurv, " KCurv=",self.kcurvature,

        if (self.kcurvature == 0):
            print " Flat "
        elif (KCurvature() < 0):
            print " Open "
        else:
            print " Closed "

        print " H0 = ",self.H0(),"\n\n"
        
        
        
