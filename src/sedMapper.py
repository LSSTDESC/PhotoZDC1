"""

  Classes that return SED "parameters" given input galaxy properties (e.g. z, phot_sim, phys_pars)

  To use them this will be the process:

  import sedMapper

  sedMap = getattr(sedMapper, "name_of_model")
  sedMap(pos_args).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the SED mapper model to be used
  
  The  model classes must all follow the same template (will need
  to make an abstract base class to ensure this is followed). They 
  must all take the same number of arguments and have the same methods defined.
"""
