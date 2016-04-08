

class Reddening(object):

    def __init__():
        #Ebmv_=Ebmv;

         
    def calzetti(lam):
        """Returns the Calzetti reddening law \f$ k(\lambda) \f$ for starburst 
           galaxies. Reddening can then be applied to a galaxy spectrum, \f$ S(\lambda) \f$,
           via: \f$ S^{red}(\lambda) = S(\lambda)10^{-0.4k(\lambda)E(B-v)} \f$
        
           @param lam   wavelength in angstroms (in rest-frame of galaxy)   
        """
        # The Calzetti law is given by :
        # k(λ) = 2.659(−1.857 + 1.040/λ) + Rv,	if 0.63 μm ≤ λ ≤ 2.20 μm
        # k(λ) = 2.659(−2.156 + 1.509/λ − 0.198/λ^2+ 0.011/λ^3 ) + Rv if 0.12 μm ≤ λ ≤ 0.63 μm
        # with Rv =4.05 and E(B−V) varies from 0 to 0.3.

        #if (lambda>2.20e-6||lambda<0.12e-6):
        #     print "WARNING! Calzetti law is undefind for lambda =", lambda
		
		# convert wavelength to m
		lam/=1e10
		
		Rv = 4.05
        lambda_micron = lam*1e6
        
        if (lam>=0.63e-6):
            k = 2.659*(-1.857 + 1.040/lambda_micron) + Rv
        elif (lam<0.63e-6):
            k = 2.659*(-2.156 + 1.509/lambda_micron - 0.198/pow(lambda_micron,2) + 0.011/pow(lambda_micron,3) ) + Rv
        else:
            raise ValueError("ERROR! don't understand wavelength value");

        return k


    def cardelli(lam, Rv):
        """Returns the Cardelli reddening law \f$ k(\lambda) \f$ for elliptical and spiral galaxies. Reddening can then be applied to a 
	         galaxy spectrum, \f$ S(\lambda) \f$,
           via: \f$ S^{red}(\lambda) = S(\lambda)10^{-0.4k(\lambda)E(B-v)} \f$
        
           @param lam   wavelength in angstroms (in rest-frame of galaxy)        
        """
        
        # Ebmv is the B-V color excess
        # The Cardelli law:
        # k(λ) = a(x) + b(x)/Rv
        # where x = 1/λ, and a and b take different forms depending on λ :
        # a(x) = 0.574x^1.61 
        # b(x) = −0.527x^1.61
        # if 0.3μm−1 ≤ x ≤ 1.1μm−1
        # y=x-1.82
        # a(x) = 1 + 0.17699y − 0.50447y^2 − 0.02427y^3 + 0.72085y^4 + 0.01979y^5 − 0.77530y^6 + 0.32999y^7
        # b(x) = 1.41338y + 2.28305y^2 + 1.07233y^3 − 5.38434y^4 − 0.62251y^5 + 5.30260y^6 − 2.09002y^7 
        # if 1.1μm−1 ≤ x ≤ 3.3μm−1
        # a(x) = 1.752 − 0.31x − 0.104/[(x − 4.67)^2 + 0.341] + Fa(x) 
        # b(x) = −3.090 + 1.825x + 1.206/[(x − 4.62)^2 + 0.263] + Fb(x)
        # Fa(x) = −0.04473(x − 5.9)^2 − 0.009779(x − 5.9)^3	8μm−1 ≥ x ≥ 5.9μm−1
        # 0 x ≤ 5.9μm−1
        # Fb(x) = 0.2130(x − 5.9)2 + 0.1207(x − 5.9)^3 8μm−1 ≥ x ≥ 5.9μm−1
        # 0 x ≤ 5.9μm−1
        # The value of Rv is usually taken to be equal to 3.1 and 
        # E(B − V ) varies from 0 to 0.1 for early type galaxies and from 0 to 0.3 for late types.


        # if (lambda>8e-6||lambda<0.3e-6):
            # print "WARNING! Cardelli law is undefind for lambda =", lambda
	
        # convert wavelength to m
		lam/=1e10
		
        lambda_micron = lam*1e6
        x = 1/lambda_micron
        y = x - 1.82
        check = 0
	
        # INFRARED
        if (x<1.11):
	
            if (check!=0):
			    print "ERROR!"
            check = 1
            a = 0.574*pow(x, 1.61) 
            b = -0.527*pow(x, 1.61)

        # OPTICAL
        elif ( (x>=1.11) and (x<3.33) ):
	
            if (check!=0):
                print "ERROR!"
            check = 1
            
            a = 1 + 0.17699*y - 0.50447*pow(y,2) - 0.02427*pow(y,3) + 0.72085*pow(y,4) + 
		            0.01979*pow(y,5) - 0.77530*pow(y,6) + 0.32999*pow(y,7)
            b = 1.41338*y + 2.28305*pow(y,2) + 1.07233*pow(y,3) - 5.38434*pow(y,4) - 
		            0.62251*pow(y,5) + 5.30260*pow(y,6) - 2.09002*pow(y,7)
		# UV 
		elif ( (x>=3.33) and (x<=8) ):
		    
		    if(check!=0):
		        print "ERROR!"
            check = 1
		
		    if ( (x>=5.9) and (x<=8) ):
                Fa = -0.04473*pow(x - 5.9,2) - 0.009779*pow(x - 5.9,3)
                Fb = 0.2130*pow(x - 5.9,2) + 0.1207*pow(x - 5.9,3) # CHECK THIS LINE!
            else:
                Fa = 0
                Fb = 0
                
            a = 1.752 - 0.31*x - 0.104/(pow(x - 4.67,2) + 0.341) + Fa
            b = -3.090 + 1.825*x + 1.206/(pow(x - 4.62,2) + 0.263) + Fb
		
        # FAR-UV
        elif (x>8):
            
            if (check!=0):
                print "ERROR!"
                
            check = 1
            
            a = -1.073 - 0.628*(x-8) + 0.137*pow(x-8,2) - 0.070*pow(x-8,3)
            b = 13.670 + 4.257*(x-8) - 0.420*pow(x-8,2) + 0.374*pow(x-8,3)
            
        else:
            raise ValueError("ERROR! don't understand wavelength value")
            
        k = a + b/Rv
        return k


# The Calzetti law is given by :
# k(λ) = 2.659(−1.857 + 1.040/λ) + Rv,	if 0.63 μm ≤ λ ≤ 2.20 μm
# k(λ) = 2.659(−2.156 + 1.509/λ − 0.198/λ^2+ 0.011/λ^3 ) + Rv if 0.12 μm ≤ λ ≤ 0.63 μm
# with Rv =4.05 and E(B−V) varies from 0 to 0.3.

# The Cardelli law:
# k(λ) = a(x) + b(x)/Rv
# where x = 1/λ, and a and b take different forms depending on λ :
# a(x) = 0.574x^1.61 
# b(x) = −0.527x^1.61
# if 0.3μm−1 ≤ x ≤ 1.1μm−1
# a(x) = 1 + 0.17699y − 0.50447y^2 − 0.02427y^3 + 0.72085y^4 + 0.01979y^5 − 0.77530y^6 + 0.32999y^7
# b(x) = 1.41338y + 2.28305y^2 + 1.07233y^3 − 5.38434y^4 − 0.62251y^5 + 5.30260y^6 − 2.09002y^7 
# if 1.1μm−1 ≤ x ≤ 3.3μm−1
# a(x) = 1.752 − 0.31x − 0.104/[(x − 4.67)^2 + 0.341] + Fa(x) 
# b(x) = −3.090 + 1.825x + 1.206/[(x − 4.62)^2 + 0.263] + Fb(x)
# Fa(x) = −0.04473(x − 5.9)^2 − 0.009779(x − 5.9)^3	8μm−1 ≥ x ≥ 5.9μm−1
# 0 x ≤ 5.9μm−1
# Fb(x) = 􏰞0.2130(x − 5.9)2 + 0.1207(x − 5.9)3 8μm−1 ≥ x ≥ 5.9μm−1
# 0 x ≤ 5.9μm−1
# The value of Rv is usually taken to be equal to 3.1 and 
# E(B − V ) varies from 0 to 0.1 for early type galaxies and from 0 to 0.3 for late types.
