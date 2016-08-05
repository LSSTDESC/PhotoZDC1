"""Functions concerned with reading in photo-z catalogs 
   in FITS format and returning the catalogs in format 
   required for certain photo-z algorithms

"""


import interface


def convert_fits_to_skynet(training_file, test_file, input_cols, output_cols, outroot):
    """Convert data stored in a FITS file to the format required to 
       run in Skynet
       
       @param training_file    filename that contains training data
       @param test_file        filename that contains test data
       @param input_cols       columns names to be the "inputs" to the model (list)
       @param output_cols      column names to be the "outputs" to the model (list)
       @parma outroot          output root name for training and test files
    
    """
    n_inputs = len(input_cols)
    n_outputs = len(output_cols)
    
    # open files for writing
    ftrain = open(outroot + "_train.txt", 'w')
    ftest = open(outroot + "_test.txt", 'w')
    
    # SkyNet requirement:
    # First line should have nin (no. of inputs), 
    # Second line should have nout (no. of outputs),
    
    # write first two lines nin, nout
    ftrain.write(str(n_inputs) + ',\n')
    ftrain.write(str(n_outputs) + ',\n')
    
    ftest.write(str(n_inputs) + ',\n')
    ftest.write(str(n_outputs) + ',\n')
    
    # Read FITS files into pandas dataframes
    training_set = interface.ReadCosmoSim(training_file)
    training_df = training_set._data
    
    test_set = interface.ReadCosmoSim(training_file)
    test_df = test_set._data
    
    # write data
    write_skynet_data(ftrain, training_df, input_cols, output_cols)
    write_skynet_data(ftest, test_df, input_cols, output_cols)      
            
            
def write_skynet_data(f, df, input_cols, output_cols):
    """Writes input and output data (after first two lines) 
       from the dataframe into the format required by SkyNet:
       
       Data files should have comma separated entries.
       
       After the first two lines, rest of the file consists 
       of data entries with inputs for a given entry on one 
       line & outputs on the line below the inputs. 
    """

    # Loop down file writing each line
    for index, row in df.iterrows():
    
        for col in input_cols:
            f.write(str(row[col]) + ',')
        f.write('\n')
        
        for col in output_cols:
            f.write(str(row[col]) + ',')
        f.write('\n')

    f.close()
    
    
def convert_fits_to_bpz(training_file, test_file, outroot, filters, m0, zs, 
                        zp_errors=0.01, zp_offsets=0., testHasZ=True):
    """Convert data stored in a FITS file to the format required to 
       run in BPZ. Creates both a data text file and the .columns file
       
       @param training_file   filename that contains training data
       @param test_file       filename that contains test data
       @parma outroot         output root name for training and test files
       @param filters         list of filter info tuple (filtname(str), mag(str), errmag(str))
                              filtname is the name of the filter in BPZ, 
                              mag, errmag are the column names of the corresponding mag, error on mag in file
       @param m0              column name of magnitude in prior
       @param zs              column name of true/spec redshift
       @param zp_errors       value(s) of zero point error in each filter (float or list/array)
       @param zp_offsets      value(s) of zero point offset in each filter (float or list/array)
       @param testHasZ        test file has redshifts
       
       BPZ data file will be written like:
       mag_cols, err_mag_cols, zs, other[id]
       
       or 
       
       mag_cols, err_mag_cols, zother[id]
       for the test file if it has no redshifts
    """
    
    # Read FITS files into pandas dataframes
    training_set = interface.ReadCosmoSim(training_file)
    training_df = training_set._data
    
    test_set = interface.ReadCosmoSim(test_file)
    test_df = test_set._data
    
    
    # extract filter info
    filt_names = [f[0] for f in filters]
    mag_cols = [f[1] for f in filters]
    emag_cols = [f[2] for f in filters]
    
    
    # open files to write to
    ftrain = open(outroot + "_train.txt", 'w')
    ftest = open(outroot + "_test.txt", 'w')
    
    
    # create .columns file
    create_columns_file(outroot, filters, m0, True, zp_errors, zp_offsets)
    
    
    # write the BPZ data for the training set
    write_bpz_data(ftrain, training_df, mag_cols, emag_cols, zs)
    
    # write the BPZ data for the test set
    if (testHasZ):
        write_bpz_data(ftest, test_df, mag_cols, emag_cols, zs)
    else:
        create_columns_file(outroot + '_testset', filters, m0, testHasZ, 
                            zp_errors, zp_offsets)
                            
        write_bpz_data(ftest, test_df, mag_cols, emag_cols, zs=None)
    
    
def create_columns_file(outroot, filters, m0, hasZ=True, zp_errors=0.01, zp_offsets=0.):
    """Write BPZ .columns file (assumes AB mags). 
       No facility yet to write multiple OTHER cols
    
       @param outroot     output root name [outroot.columns]
       @param filters     list of tuples of (filtname(str), mag_col#(int), magerr_col#(int))
       @param m0          column name of magnitude in prior
       @param hasZ        if file will have true/spec redshifts
       @param zp_errors   value(s) of zero point error in each filter (float or list/array)
       @param zp_offsets  value(s) of zero point offset in each filter (float or list/array)
    """
    
    
    # convert zp_errors, zp_offsets into list if needed
    num_filters = len(filters)
    zp_errors = parse_zp_pars(zp_errors, num_filters)
    zp_offsets = parse_zp_pars(zp_offsets, num_filters)
    
    # open file to write to
    f_columns = open(outroot + '.columns', 'w')
    
    # write the header
    header = "# Filter    columns    AB/Vega    zp_error    zp_offset\n"
    f_columns.write(header)
    
    im0 = -1
    for i, filt in enumerate(filters):
        
        filt_name = filt[0]
        mag_col = filt[1]
        magerr_col = filt[2]
        
        if mag_col == m0:
            im0 = i+1
        
        f_columns.write(filt_name + '    ')
        f_columns.write(str(i+1) + ',' + str(i+num_filters+1) + '    AB    ')
        f_columns.write(str(zp_errors[i]) + '    ' + str(zp_offsets[i]) + '\n')
        
    if (im0<0):
        raise ValueError("Magnitude prior column not found!")
    
    # prior magnitude column
    f_columns.write('M_0    ' + str(im0) + '\n')
    
    # redshift column, other[id] column
    if (hasZ):
        f_columns.write('Z_S    ' + str(num_filters*2+1) + '\n')
        f_columns.write('OTHER    ' + str(num_filters*2+2) )
    else:
        f_columns.write('OTHER    ' + str(num_filters*2+1) )


def parse_zp_pars(zp_pars, num_filters):
    """Parse the zero point parameters, either zp_error
       or zp_offset to be set in BPZ .columns file
       
       @param zp_pars       either value(s) to set for zp_error or zp_offset
       @param num_filters   number of filters
       
    """
    
    # check if list/array or not
    isList = False
    try: 
        np = len(zp_pars)
        isList = True
    except:
        # if not list, check is float
        if isinstance(zp_pars, float):
    
            # if float, turn into list
            zp_par = zp_pars
            zp_pars = []
            for i in range(num_filters):
                zp_pars.append(zp_par)  
        else:
            raise ValueError("data type of zp_pars not understood")
    
    
    if isList and len(zp_pars) != num_filters:
        raise ValueError("zp_pars length doesn't match filter length")
    
    
    return zp_pars
            
        
            
def write_bpz_data(f, df, mag_cols, emag_cols, zs=None):
    """Writes data from the dataframe into the format required by BPZ
    
       @param f           open file
       @param df          dataframe containing data
       @param mag_cols    column names in df corresponding to magnitudes
       @param emag_cols   column names in df corresponding to magnitude errors
       @param zs          column name in df corresponding to true/spec redshift
    
       BPZ data file will be written like:
       mag_cols, err_mag_cols, zs, other[id]
       
       or 
       
       mag_cols, err_mag_cols, zother[id]
       if it has no redshifts
    """
    
    if zs is not None:
        hasZ = True
    else:
        hasZ = False

    # Loop down file writing each line
    for index, row in df.iterrows():
    
        for col in mag_cols:
            f.write(str(row[col]) + ' ')
        
        for col in emag_cols:
            f.write(str(row[col]) + ' ')
            
        if hasZ:
            f.write(str(row[zs]) + ' ')
        
        f.write(str(index) + ' ')
            
        f.write('\n')

    f.close()      
