import interface


def convert_fits_to_skynet(training_file, test_file, input_cols, output_cols, outroot):
    """Convert data stored in a FITS file to the format required to 
       run in Skynet
       
       @param training_file    filename that contains training data
       @param test_file        filename that contains test data
       @param input_cols       columns names to be the "inputs" to the model (list)
       @param output_cols      column names to be the "outputs" to the model (list)
    
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
