"""
to add AtomicCalcs and VoigtProfile and similar generic astro-based calculations

"""

import const

class AtomicCalcs(object):

    def __init__(self):
        """Set relevant constants """
    
        self.SIGMA_LYMANLIM_CMSQ = 6.30e-18    # In cm^2
        self.WAVE_LYMANLIM_METERS =  91.175e-9 # in m
        self.FREQ_LYMANLIM_INVSEC = const.cms/self.WAVE_LYMANLIM_METERS
    
        self.LINE_NUM_LYALPHA = 2
        self.LINE_NUM_LYMAN_MAX = 32
        self.NGAMMA_MAX = 24
        
        
    def getWavelengthLymanSeries(self, lineNum):
        """Get the wavelength of the Lyman line number in meters"""
        
        if (lineNum<self.LINE_NUM_LYALPHA):
            raise ValueError("ERROR! Lyman line number cannot be less than 2")
            
        if (lineNum>self.LINE_NUM_LYMAN_MAX):
            raise ValueError("ERROR! Lyman line number can't be larger than " + str(self.LINE_NUM_LYMAN_MAX) )
            
        return self.WAVE_LYMANLIM_METERS/(1. - 1./(lineNum*lineNum))
