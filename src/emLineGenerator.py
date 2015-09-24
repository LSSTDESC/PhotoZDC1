"""

  Classes that generate emission line(s)

  These classes will be instantiated inside another class that generates galaxy 
  observations. They will take as an argument the bulk of all the galaxy
  data contained within that class. Each particular class will access which
  attributes of the bulk of data it needs to do its emission line generation.

  To use them this will be the process:

  import emLineGenerator

  eLG = getattr(emLineGenerator, "name_of_model")
  elG(gal_data).method_that_does_stuff()

  where "name_of_model" will correspond to the name of the
  class that describes the emission line model to be used

  So within method_that_does_stuff() the data will get accessed like:
  
  gal_data["mstellar"]
  gal_data["z"]

  etc

  Therefore in gal_data there can be data the emission line class
  *doesn't* use, but it must have the data to be accessed, all with the
  same names

  
  The emission line classes must all follow the same template (will need
  to make an abstract base class to ensure this is followed). They 
  must all take the same number of arguments and have the same methods defined.


"""
