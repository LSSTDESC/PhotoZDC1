"""
  
  Classes that return SED given SED "parameters" that were input

  To use them this will be the process:

  import sedGenerator

  sedGen = getattr(sedGenerator, "name_of_model")
  sedGen(sed_pars).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the SED generator model to be used
  
  The  model classes must all follow the same template (will need
  to make an abstract base class to ensure this is followed). They 
  must all take the same number of arguments and have the same methods defined.



"""
