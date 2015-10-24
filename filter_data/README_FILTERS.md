# Filter transmission functions

This directory will contain the filter transmission data. Currently it is in the form of text files (no 
special header format) with two columns: first column is wavelengths in **ANGSTROMS** and second column is the
filter transmission (no units, value between 0 and 1). 

ALL FILTER NAMES MUST FOLLOW THIS FORMAT:
[FILTERSETNAME]_[FILTERBAND].[EXTENSION]
e.g.
SDSS_u.dati
TOPHAT_UV2800.dati
LSST_u.dat

Each filter is contained within its own text file. The name of the text file is important as it is used to 
refer to the filter throughout the code. The wavelengths of the filters in the file should begin and end at 
the filter edges (to avoid unnecessary integration).

The `*.filter` files are the list of all filters belonging to a particular filter set and are used by the code
to read in the desired set of filters.

The LSST filter transmissions were obtained from the [LSST github repo](https://github.com/lsst/throughputs)

The Johnson filter transmissions are from BPZ

Other filter transmissions are from the Theoretical Astrophysical Observatory (TAO)
https://tao.asvo.org.au/tao/

UKIRT_Y_filter.dati
UKIRT_K_filter.dati
UKIRT_J_filter.dati
UKIRT_H_filter.dati
TOPHAT_UV2800top.dati
TOPHAT_UV2300top.dati
TOPHAT_UV1500top.dati
SDSS_sdss_z.dati
SDSS_sdss_u.dati
SDSS_sdss_r.dati
SDSS_sdss_i.dati
SDSS_sdss_g.dati
PACS_pacs70.dati
PACS_pacs160.dati
PACS_pacs100.dati
IRAC_irac_ch4.dati
IRAC_irac_ch3.dati
IRAC_irac_ch2.dati
IRAC_irac_ch1.dati
GALEX_galex_NUV.dati
GALEX_galex_FUV.dati
2MASS_Ksband_2mass.dati
2MASS_Jband_2mass.dati
2MASS_Hband_2mass.dati
