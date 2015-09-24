"""
  Higher level objects

  Will probably have the following classes:

  CosmoSimsGalaxyObj

  GalaxyObj

  ObservedGalaxyObj


"""
# will import these once implemented
#import sedMapper
#import sedGenerator

class CosmoSimsGalaxyObj(object):
    """ Imagine this object as holding one row from CosmoSims table
    
    """

    def __init__(self, galaxy_properties):
        pass
        
        
    def compute_sim_color(self):
        """ returns color of CosmoSims galaxy
        """
        pass
        
    
class GalaxyObject(object):

    def __init__(self, galaxy_properties, filterDict, cosmoModel, igmModel, pars):
        """Initialise all the galaxy data:
           galaxy_properties = z, phot_sim, phys_pars, size (from CosmoSims table). Maybe CosmoSimsGalaxyObj
                               instead?
           filterDict = list of filters to calculate photometry in
           cosmoModel = cosmological parameters and basic calculations
           igmModel = model for attenuation by IGM
           pars = i.e. which SEDMapper, SEDGenerator, EmLine module to use (actually maybe IGM too)
        """
        pass
        
    def genSED(self, z, phot_sim, phys_pars):
        """Use SEDMapper model to map z, phys_pars to galaxy SED via SEDGenerator
           [probably won't need these args as they will be attributes] 
           
        """
        # SEDMapper generates "SED parameters" TBD given input galaxy data (z, phot_sim, phys_pars) 
        
        # SEDGenerator generates SED given "SED parameters" TBD
        
        # once generated SED, store SED object as class attribute
        
        # add IGM model to SED as well
        pass
        
        
    def genEmLine(self, z, phot_sim, phys_pars):
        """Use EmLine model to generate emission lines to add to this galaxy
           [probably won't need these args as they will be attributes]   
        """
        
        # generate emission line(s) using EmLine model
        
        # add emission lines directly to SED object stored in this class as attribute
        pass
        
        
    def calcPhot(self):
        """Calculate true magnitudes in filterDict given SED
        
        """
        pass


class ObservedGalaxyObject(GalaxyObject):
    """Add here observational effects, Galactic reddening, photometric errors, filter observation effects"""

    def __init__(self):
        pass
        
        
        
