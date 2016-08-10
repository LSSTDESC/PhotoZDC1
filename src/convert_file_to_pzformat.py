"""Functions concerned with reading in photo-z catalogs 
   and returning the catalogs in file format required 
   for certain photo-z algorithms

"""


import interface


#### SKYNET ####

def convert_fits_to_skynet(training_file, test_file, input_cols, output_cols, outroot):
    """Convert data stored in a FITS file to the format required to 
       run in Skynet
       
       @param training_file    filename that contains training data
       @param test_file        filename that contains test data
       @param input_cols       columns names to be the "inputs" to the model (list)
       @param output_cols      column names to be the "outputs" to the model (list)
       @parma outroot          output root name for training and test files
    
    """
    
    # file names for writing
    ftrain = outroot + "_train.txt"
    ftest = outroot + "_test.txt"
    
    
    # Read FITS files into pandas dataframes
    training_set = interface.ReadCosmoSim(training_file)
    training_df = training_set._data
    
    test_set = interface.ReadCosmoSim(test_file)
    test_df = test_set._data
    
    
    # write data
    write_skynet_data(ftrain, training_df, input_cols, output_cols)
    write_skynet_data(ftest, test_df, input_cols, output_cols)      
            
            
def write_skynet_data(outfile, df, input_cols, output_cols):
    """Writes input and output data (after first two lines) 
       from the dataframe into the format required by SkyNet:
       
       Data files should have comma separated entries.
       
       SkyNet requirement:
       First line should have nin (no. of inputs), 
       Second line should have nout (no. of outputs),
       
       After the first two lines, rest of the file consists 
       of data entries with inputs for a given entry on one 
       line & outputs on the line below the inputs. 
    """
    # open file to write to
    f = open(outfile, 'w')
    
    n_inputs = len(input_cols)
    n_outputs = len(output_cols)
    
    # write first two lines nin, nout
    f.write(str(n_inputs) + ',\n')
    f.write(str(n_outputs) + ',\n')

    # Loop down file writing each line
    for index, row in df.iterrows():
    
        for col in input_cols:
            f.write(str(row[col]) + ',')
        f.write('\n')
        
        for col in output_cols:
            f.write(str(row[col]) + ',')
        f.write('\n')

    f.close()
    
    
def create_inp_file(inp_file_name,
                    input_root,
                    output_root,
                    hidden_layers_nodes,
                    #output_layer_nodes,
                    classification=False,
                    hidden_layers_func=1,
                    output_layer_func=0,
                    mini_batch_fraction=1,
                    validation_data=True,
                    whitenin=1,
                    whitenout=1,
                    sigma=0.05,
                    confidence_rate=0.5,
                    confidence_rate_min=0.1,
                    max_iter=2000,
                    convergence_function=1,
                    verbose=3,
                    iteration_print_freq=50,
                    fix_seed=False,
                    fixed_seed=0,
                    prior=True,
                    noise_scaling=False,
                    set_whitened_noise=False,
                    historic_maxent=False,
                    resume=False,
                    reset_alpha=False,
                    reset_sigma=False,
                    randomise_weights=0.01,
                    line_search=0,
                    calculate_evidence=False,
                    pretrain=False,
                    nepoch=10,
                    recurrent=False,
                    norbias=False):
                                  
                    
    """Write SkyNet parameters file
    
       @param inp_file_name            filename for parameter file
       @param input_root               /path/to/[filestem][train|test].txt
       @param output_root              /path/to/[filestem] (where output files go)
       @param classification           False=do regression, True=do classification
       @param hidden_layers_nodes      list of number of nodes in each hidden layer
       @param output_layer_nodes       number of nodes in output layer
       @param hidden_layers_func       set activation function of hidden layer connections
                                       0=linear, 1=sigmoid, 2=tanh,
                                       3=rectified linear, 4=softsign
                                       if int all hidden layers have same function
                                       if list, list must have len=len(hidden_layers_nodes)
       @param output_layer_func        as above, but for output layer therefore cannot be list
       @param mini_batch_fraction      fraction of training data to be used in each batch
       @param validation_data          is there validation data to test against?
       @param whitenin                 input whitening transform to use
                                       0=none, 1=[0,1] range, 2=standard normal dist
       @param whitenout                as above for output transform  
       @param sigma                    initial noise level, set on (un-)whitened data
       @param confidence_rate          initial learning rate step size factor, 
                                       higher values are more aggressive
       @param confidence_rate_min      minimum confidence rate allowed
       @param max_iter                 max no. of iterations allowed
       @param convergence_function     function to use for convergence testing:
                                       1=log-posterior, 2=log-likelihood, 
                                       3=correlation, 4=error squared
       @param verbose                  verbosity level of feedback sent to stdout (0=min, 3=max)
       @param iteration_print_freq     stdout feedback frequency
       @param fix_seed                 use a fixed seed?
       @param fixed_seed               seed to use
       
    
       PARAMETERS BELOW LESS COMMON TO CHANGE
       @param prior                    use L2 weight regularization. Strongly advised
       @param noise_scaling            if noise level (std of outputs) is to be estimated
       @param set_whitened_noise       if the noise is to be set on whitened data
       @param historic_maxent          experimental implementation of MemSys's historic maxent option
       @param resume                   resume from a previous job
       @param reset_alpha              reset hyperparameter upon resume
       @param reset_sigma              reset hyperparameters upon resume
       @param randomise_weights        random factor to add to saved weights upon resume
       @param line_search	           perform line search for optimal distance
                                       0=none, 1=golden section, 2=linbcg lnsrch
       @param calculate_evidence       whether to calculate the evidence at the convergence
       @param pretrain                 perform pre-training using restricted BM?
       @param nepoch                   number of epochs to use in pre-training
       @param recurrent                use a RNN?
       @param norbias                  use a bias for the recurrent hidden layer connections?
       
    """

    f = open(inp_file_name, 'w')
    
    # where to read files from, where to write files to
    f.write('#input_root\n')
    f.write(input_root + '\n')
    f.write('#output_root\n')
    f.write(output_root + '\n')

    # set hidden layer parameters
    nhidden = len(hidden_layers_nodes)
    if isinstance(hidden_layers_func, list):
        if  len(hidden_layers_func) != nhidden:
            raise ValueError("len(hidden_layers_nodes) != len(hidden_layers_func)")
    else:
        hidden_layers_func = [hidden_layers_func for i in range(nhidden)]
    
    activation = ''
    for i,nhid in enumerate(hidden_layers_nodes):
        f.write('#nhid\n')
        f.write(str(nhid) + '\n')
        
        activation += str(hidden_layers_func[i])
        
    # activation function for each hidden layer and output layer
    activation += str(output_layer_func)
    f.write('#activation\n')
    f.write(activation + '\n')
    
    # set model type
    f.write('#classification_network\n')
    if (classification):
        f.write('1\n')
    else:
        f.write('0\n')
        
    # mini-batch fraction of data
    if mini_batch_fraction>1 or mini_batch_fraction<0:
        raise ValueError("mini_batch_fraction<0 or >1")
    f.write('#mini-batch_fraction\n')
    f.write(str(mini_batch_fraction) + '\n')
    
    # is there validation data
    f.write('#validation_data\n')
    if (validation_data):
        f.write('1\n')
    else:
        f.write('0\n')
        
    # data whitening
    if whitenin not in [0,1,2]:
        raise ValueError("whiten input choice not understood")
    f.write('#whitenin\n')
    f.write(str(whitenin) + '\n')
    if whitenout not in [0,1,2]:
        raise ValueError("whiten output choice not understood")
    f.write('#whitenout\n')
    f.write(str(whitenout) + '\n')
    
    # initial noise level
    f.write('#sigma\n')
    f.write(str(sigma) + '\n')
    
    # confidence rate/learning rate
    f.write('#confidence_rate\n')
    f.write(str(confidence_rate) + '\n')
    f.write('#confidence_rate_minimum\n')
    f.write(str(confidence_rate_min) + '\n')
    
    # maximum number of iterations
    f.write('#max_iter\n')
    f.write(str(max_iter) + '\n')
    
    # convergence_function 
    f.write('#convergence_function\n')
    f.write(str(convergence_function) + '\n')
    
    # printing to stdout
    f.write('#verbose\n')
    f.write(str(verbose) + '\n')
    f.write('#iteration_print_frequency\n')
    f.write(str(iteration_print_freq) + '\n')
    
    # fix seed?
    if (fix_seed):
        f.write('#fix_seed\n')
        f.write('1\n')
        f.write('#fixed_seed\n')
        f.write(str(fixed_seed) + '\n')
    else:
        f.write('#fix_seed\n')
        f.write('0\n')
        
    # prior
    f.write('#prior\n')
    if (prior):
        f.write('1\n')
    else:
        f.write('0\n')
             
    # estimate noise level of outputs
    f.write('#noise_scaling\n')
    if (noise_scaling):
        f.write('1\n')
    else:
        f.write('0\n')
        
    # if set noise on whitened data
    f.write('#set_whitened_noise\n')
    if (set_whitened_noise):
        f.write('1\n')
    else:
        f.write('0\n')  
    
    # ?    
    f.write('#historic_maxent\n')
    if (historic_maxent):
        f.write('1\n')
    else:
        f.write('0\n') 
        
    # if resuming from a previous job
    f.write('#resume\n')
    if (resume):
        f.write('1\n')
        
        f.write('#reset_alpha\n')
        if (reset_alpha):
            f.write('1\n')
        else:
            f.write('0\n')
            
        f.write('#reset_sigma\n')
        if (reset_sigma):
            f.write('1\n')
        else:
            f.write('0\n')
            
        f.write('#randomise_weights\n')
        if (randomise_weights):
            f.write('1\n')
        else:
            f.write('0\n')
            
    else:
        f.write('0\n') 
        
    # line search for optimal distance
    if line_search not in [0,1,2]:
        raise ValueError("line_search choice not understood")
    f.write('#line_search\n')
    f.write(str(line_search) + '\n')
    
    # calculate evidence  
    f.write('#calculate_evidence\n')
    if (calculate_evidence):
        f.write('1\n')
    else:
        f.write('0\n') 
        
    # pre-training
    f.write('#pretrain\n')
    if (pretrain):
        f.write('1\n')
        
        f.write('#nepoch\n')
        f.write(str(nepoch) + '\n')
    else:
        f.write('0\n') 

    # RNN
    f.write('#recurrent\n')
    if (resume):
        f.write('1\n')
        
        f.write('#norbias\n')
        if (norbias):
            f.write('1\n')
        else:
            f.write('0\n')
    else:
        f.write('0\n')
     
#### BPZ ####
    
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
    

