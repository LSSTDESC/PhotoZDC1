"""Classes and functions to interface with data sets, our "Butler"

  

"""

import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import pandas as pd
import operator

# operator look up dictionary
ops = {"==": operator.eq, "!=": operator.ne, ">": operator.gt, ">=": operator.ge, "<": operator.lt,
       "<=": operator.le}


class ReadCosmoSim(object):
    """Read cosmological simulation from a FITS or text file
    
       Text file must contain column names that can be parsed by numpy.genfromtxt().
       FITS file must have extension .fit or .fits
    
    """

    def __init__(self, file_to_read, nhdu=1, delimiter=' ', index=False):
        """Read in data from file
        
           @param file_to_read   full path and name of file containing simulation
           @param nhdu           (FITS ONLY) HDU number containing data (indexed from zero)
           @param delimiter      (TEXT ONLY) character that delimits text file
           @param index          (TEXT ONLY) column number to use as row index (default is none)
        """
        ext = file_to_read.split(".")[-1]
        
        if (ext=="fits" or ext=="fit"):
            self._filetype = "FITS"
            self._read_fits(file_to_read, nhdu)
        else:
            self._filetype = "TEXT"
            self._read_text(file_to_read)


    def _read_fits(self, file_to_read, nhdu=1):
        """Read data and column names from a FITS file

           @param file_to_read   full path and name of file containing simulation
           @param nhdu           HDU number containing data (indexed from zero)
        """
        
        table = Table.read(file_to_read)
        data = table._data # will be deprecated but Table.as_array doesn't work????
        # Fix byte order.
        # See https://github.com/astropy/astropy/issues/1156
        data = data.byteswap().newbyteorder()
        self._data = pd.DataFrame.from_records(data)
        self._col_names = self._data.columns
        
        # old attempt
        #hdulist = fits.open(file_to_read)
        #col_names = [x.name for x in hdulist[nhdu].data.columns]
        #self._data = pd.DataFrame(np.asarray(hdulist[nhdu].data), columns=col_names)
        #self._data = pd.DataFrame.from_records(np.asarray(hdulist[nhdu].data))
        #self._data.columns = col_names
        #self._col_names = self._data.columns
        # ValueError: Big-endian buffer not supported on little-endian compiler
        
        
    def _read_text(self, file_to_read, delimiter=' ', index=False):
        """Read data and column names from a text file

           @param file_to_read   full path and name of file containing simulation
           @param delimiter      delimiter in text file
           @param index          column to use as row index (False=none)
        """
        self._data = pd.read_csv(file_to_read, delimiter=delimiter, index_col=index) 
        #np.genfromtxt(file_to_read, dtype='float', names=True) # this is WAY slower
        self._col_names = self._data.columns

        
    def get_number_of_cols(self):
        """Return number of columns in file
        """
        return len(self._col_names)
        
        
    def get_number_of_rows(self):
        """Return number of entries in file
        """
        return len(self._data)
        
    def get_dimensions(self):
        """Return the size of each dimension in the file
        """
        return self._data.shape
        
    
    def get_column_names(self):
        """Return column names
        """
        return self._col_names
        
        
    def get_column(self, column, selection=None):
        """Return column via its index number or name
        
           @param column       column name (str) or index (int)
           @param selection    perform selection on column tuple describing condition on column
                               e.g. (">", 20.)
        """
        
        self._check_column_valid(column)
        
        if (selection==None):
            return self._data[column]
        else:
            condition = [(column, selection[0], selection[1])]
            return self.select_given_all_true(condition, cols_to_select=[column])
        
        
    # could make data selection method
    # args: cols_to_select, conditions (tuple of: (col, cond, val) e.g. (colname, '>', 10)
    #    df[ (df.a > 100) & (df.b == "test") ]
    def select_given_all_true(self, conditions, cols_to_select='all'):
        """Select all the listed columns and return all entries where all the conditions are true. Returns
           dataframe
        
           @param conditions       list of tuples describing conditions upon each column
                                   (colname/index, 'operator', value)
           @param cols_to_select   list of column names OR list of column numbers OR select all columns
          
           For example:
           select_given_all_true([("uband", "<", 21), ("redshift", ">=", 0.1)])
           select_given_all_true([("uband", "<", 21), ("redshift", ">=", 0.1)], ["ra", "dec"])
           select_given_all_true([(10, "<", 21), (1, ">=", 0.1)], [2, 3])
        """
        
        # check all the column names to return
        if isinstance(cols_to_select, list):
            for column in cols_to_select:
                self._check_column_valid(column)
        
        # make condition, start with everything true
        final_condition = pd.Series(np.ones(self._data.shape[0], dtype=bool))
        for condition in conditions:
                        
            condition_col = condition[0]
            op = condition[1]
            val = condition[2]
            
            print "Adding condition that column: ", condition_col , op , val
            final_condition = operator.and_( final_condition, ops[op](self._data[condition_col],val) )

        if isinstance(cols_to_select, list):
            return self._data[final_condition][cols_to_select]
        else:
            return self._data[final_condition]
            
            
    def group_by_column(self, column, bins):
        """Return a dataframe grouped into bins based upon a column
           
           @param column     column name (str) or index (int)
           @param bins       bin edges to be used to group the data in column
        """ 
        self._check_column_valid(column)
        
        return self._data.groupby(pd.cut(self._data[column], bins=bins))
        
        
    def group_by_color(self, column1, column2, bins):
        """Return a dataframe grouped into bins based upon a "color" (difference between two columns)
           
           @param column     column name (str) or index (int)
           @param bins       bin edges to be used to group the data in column
        """ 
        self._check_column_valid(column1)
        self._check_column_valid(column2)
        
        return self._data.groupby(pd.cut(self._data[column1]-self._data[column2], bins=bins))
        
        
    def get_random_subsample(self, nsample, cols_to_select='all'):
        """Return a random sub-sample of the catalog of size nsample. Returns dataframe
        
           @param nsample           number of entries of catalog to select [int]
           @param cols_to_select    list of column names OR list of column numbers OR select all columns
        """
        
        # check all the column names to return
        if isinstance(cols_to_select, list):
            for column in cols_to_select:
                self._check_column_valid(column)
                
        ids = np.random.permutation(self.get_number_of_rows())[:nsample]
        if isinstance(cols_to_select, list):
            return self._data.iloc[ids][cols_to_select]
        else:
            return self._data.iloc[ids]
        
        
    def get_array_of_colors(self, mag_columns, nsample=None):
        """Return numpy array of colors by subtracting magnitude columns in order supplied. Also returns
           the standard deviation of each sample of colors
        
           @param mag_columns     names of columns containing magnitudes. Must be supplied in correct order
                                  so color1 = mag_col1 - mag_col2
                                     color2 = mag_col2 - mag_col3 etc
           @param nsample         number of entries from table to select, if less than table size, randomly
                                  selects the entries
                                  
           Returns: 
           - dataframe of colors for each galaxy in sample, column names are colors
           
        """
        
        if (nsample==None):
            nsample = self.get_number_of_rows()
        
        if (nsample<self.get_number_of_rows()):
            data = self.get_random_subsample(nsample)
        else:
            data = self._data
            
        
        # check all the column names to return
        for column in mag_columns:
            self._check_column_valid(column)
        
        color_array = np.zeros((nsample, len(mag_columns)-1))
        color_names = []
        for i in range(len(mag_columns)-1):
        
            color_array[:,i] = data[mag_columns[i]] - data[mag_columns[i+1]]
            color_names.append(str(mag_columns[i]) + "-" + str(mag_columns[i+1]) )
        
        # convert to dataframe and return
        return pd.DataFrame(color_array, columns=color_names)

        
        
    def _check_column_valid(self, column):
        """Check the column name or index is valid
           
           @param column    column name (str) or index (int)
        """
        if (isinstance(column, (int, long) )):
            if (column<0 and column>=self.get_number_of_cols()):
                raise ValueError("ERROR! column number (" + str(column) + ") not valid!")
                
        if (isinstance(column, str )):
            if (column not in self._col_names):
                raise ValueError("ERROR! column name (" + column + ") not valid!")
                
                
                
