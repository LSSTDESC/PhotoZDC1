# Filter transmission functions

This directory will contain the filter transmission data. Currently it is in the form of text files (no 
special header format) with two columns: first column is wavelengths in **ANGSTROMS** and second column is the
filter transmission (no units, value between 0 and 1). 

Each filter is contained within its own text file. The name of the text file is important as it is used to 
refer to the filter throughout the code. The wavelengths of the filters in the file should begin and end at 
the filter edges (to avoid unnecessary integration).

The `*.filter` files are the list of all filters belonging to a particular filter set and are used by the code
to read in the desired set of filters.

The LSST filter transmissions were obtained from the [LSST github repo](https://github.com/lsst/throughputs)
