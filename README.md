# Photo-z Task Force Simulations Generator


### Project aim

The photo-z Task Force simulations will be the starting point of simulating systematic effects in the LSST
photometry that impact the quality of the photo-z. This includes (eventually) simulating galaxy properties
realistically enough that we can mimic the kinds of spectroscopic training and calibration samples that will
be obtainable.

We will use baseline properties of galaxies from mock galaxy catalogs produced by the DESC simulation working
group, and these will contain LSST photometry with no systematics. Using these simulations we will be able to 
perturb the spectral energy distributions (SEDs) of galaxies according to a variety of different systematics,
e.g. dust, IGM absorption, blended galaxies, etc --- to determine how LSST photometry is affected. These
simulations should be considered a secondary "layer" on top of the cosmological galaxy simulations where we
can generate many observed galaxy catalogs with differing properties very quickly. We will have the option to 
turn on and off these various sources of systematic perturbations, and can determine which are the most 
important, and characterize their impact.

The resulting galaxy catalogs will not only be used to test the quality of the photo-z due to a variety of
systematic effects, but also to test the robustness of photo-z algorithms in producing robust probability of
redshifts *p(z)*.



### Brief Project Summary

Using mock catalogs developed in the Simulations group, the Task Force will determine a way of populating the 
sims with galaxy SEDs (e.g. the Brown templates or a continuous parameterization of them), connecting physical
parameters such as stellar mass to the SEDs.  The goal is to construct more realistic obtainable training and
calibration sets from these catalogs (e.g. via a stellar mass complete sample). As well as generating large 
photo-z catalogs, this will enable us to study training set incompleteness and its effect on the cosmological
analysis.


### Tasks 

**To be implemented!**

* define a continuous parametrization of realistic SEDs, see `sedGenerator.py`
* construct the mapping from SEDs onto simulated photometry and physical galaxy parameters, see `sedMapper.py`
* add prescriptions for incorporating the impact of emission lines, see `emLineGenerator.py`
* add prescriptions for including the variation of the line of sight IGM absorption for high redshift galaxies,
see `igmModel.py`


### Deliverables

* Catalogs of *p(z)*


### Code installation

1. python 2.X distribution with all the usual (numpy, scipy, matplotlib, etc) and the following packages: cython (optional), scikit-learn, pandas, astropy. [Anaconda](https://store.continuum.io/cshop/anaconda/) is a good choice.
2. clone repository to local directory
3. add `src/` directory to `PYTHONPATH` environment variable
4. try running some examples!


### Code organisation/directory structure

* `src/` contains all the important code, hopefully docstrings make everything clear-ish
* `examples/` demonstrates some current calculations the code can do and demonstrates how to construct programs
* `sed_data/` contains the SED data, see README file for more information
* `filter_data/` contains filter transmission functions, see README file for more information
* `igm_data/` currently empty, will contain line of sight transmission data for the IGM
* `tests/` contains unit tests


### Warnings

This code is currently under *extreme* development! 

If using python < 2.7.6 and/or scipy < 0.14 you may get a lot of printing of:
`IntegrationWarning: The maximum number of subdivisions (50) has been achieved.'
from `scipy.integrate.quad`

The integration method is under development because `scipy.quad` is too slow.
