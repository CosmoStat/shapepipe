"""! 
   @package scatalog.scatalog Simple catalog management
   @author Marc Gentile
   @file scatalog.py
   Simple catalog management for astronomical data
"""

# -- Python imports
import os, sys
import string
import re
import numpy as np
import asciidata as asc    # text catalogs
import pyfits              # FTS catalogs

# -- Internal imports
from scatalog_helper import *    # Helper methods

# --------------------------------------------------------------------------------------------------
class BaseCatalog(object):

   """! 
      Base catalog management class
      @note This class is not meant to be used directly 
   """

   def __init__(self, fullpath):
      """!  BaseCatalog class constructor """

      # Public members
      self._directory, self._filename = os.path.split(fullpath)   # catalog filepath
      self._format = BaseCatalog.InputFormat.Undefined            # input/output format
      self._helper = Helper()

      # Private members
      self._cat_data = None       # catalog internal data (e.g. AsciiData class)
      #self.__table_data = None    # internal numpy matrix    


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def fullpath(self):
      """! 
         get the full path of the catalog file 
         @return full path of the catalog file
      """
      return os.path.join(self._directory, self._filename)

   @property
   def directory(self):
      """! 
         Get the directory of the catalog file 
         @return directory of the catalog file
      """
      return self._directory

   @property
   def filename(self):
      """! 
         Get the name of the catalog file 
         @return name of the catalog file
      """
      return self._filename

   @property
   def format(self):
      """!
         Get the default input/output format of the catalog (e.g. Text, SExtractor, FITS) 
         @return the input/output format of the catalog
         @note the returned format will always be 
               scatalog.scatalog.BaseCatalog.InputFormat.Undefined for this class
      """
      return self._format

   @property
   def helper(self):
      """! 
         Get the helper class 
         @return helper class
      """
      return self._helper


   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   def get_nb_rows(self):
      """! 
         Get the number of rows in the catalog 
         @return number of rows
      """
      raise BaseCatalog.FeatureNotImplemented("get_nb_rows()")

   def get_nb_cols(self):
      """! 
         Get the number of columns in the catalog 
         @return number of columns
      """
      raise BaseCatalog.FeatureNotImplemented("get_nb_cols()")

   def get_col_names(self):
      """! 
         Get the list of column names in the catalog 
         @return list of column names
      """
      raise BaseCatalog.FeatureNotImplemented("get_col_names()")

   def get_col_formats(self):
      """!  Get the list of column formats in the order of columns """
      if not self._cat_data is None:
         return [c.get_format() for c in self._cat_data]   
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def append_col(self, column):
      """!
         Append a Column object after index @c index
         @param column an object derrived from BaseCatalog.Column
      """
      raise BaseCatalog.FeatureNotImplemented("append_col()")
      
   
   # ~~~~~~~~~~~~~
   # Inner Classes 
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class InputFormat:
      """! Supported input catalog formats"""
      (Undefined, TabulatedText, SExtractor, FITS) = (0,1,2,4)


   # -----------------------------------------------------------------------------------------------
   class Column(object):
      """! Represents a column in the catalog """

      def __init__(self, name, format=None, comment=None, data=None):

         self._name = name          # column name
         self._format = format      # column format (Python)
         self._comment = comment    # comment string
         self._data = data
         self._helper = Helper()
         #self._data = self.helper.fix_data_type(self._data)

      # ~~~~~~~~~~~
      # Properties 
      # ~~~~~~~~~~~

      # --- Getters

      @property
      def name(self):
         """! 
            Get the name of the column 
            @return name of the column
         """
         return self._name

      @property
      def format(self):
         """! 
            Get the format of the column 
            @return format of the column
         """
         return self._format

      @property
      def comment(self):
         """! 
            Get the column comment string 
            @return comment comment string
         """
         return self._comment

      @property
      def data(self):
         """! 
            Get the data associated with the column 
            @return data data associated with the column
         """
         return self._data


      @property
      def helper(self):
         """! 
            Get the helper class 
            @return helper class
         """
         return self._helper

      # --- Setters

      @name.setter
      def name(self, name):
         """! 
            Set the name of the column
            @param name of the column
         """
         self._name = name

      @format.setter
      def format(self, format):
         """! 
            Set the format of the column
            @param format new format to set
         """
         self._format = format

      @comment.setter
      def comment(self, comment):
         """! 
            Set the comment string of the column
            @param comment comment string
         """
         self._comment = comment

      @data.setter
      def data(self, data):
         """! 
            Set the data string of the column
            @param data data associated with the column
         """
         self._data = data


      # ~~~~~~~~~~~~~~
      # Public Methods 
      # ~~~~~~~~~~~~~~

      # --- Getters

      def get_nb_rows(self):
         """! Retrieve the number of rows of the column """
         raise BaseCatalog.FeatureNotImplemented("get_nb_rows()")

      def get_info(self):
         """! Retrieve information about the column """
         raise BaseCatalog.FeatureNotImplemented("get_info()")

      def get_type(self):
         """! 
            Get the data type of the column 
            @return Python data type of the column
          """
         raise BaseCatalog.FeatureNotImplemented("get_type()")


   # -----------------------------------------------------------------------------------------------
   class FeatureNotImplemented(NotImplementedError):
      """ Feature Not Implemented """

      def __init__(self, msg):
         self._msg = msg

      def __str__(self):
         return "SCatalog *** ERROR ***: Feature: {0} is not implemented in this class".format(
                                                                                         self._msg)
   # -----------------------------------------------------------------------------------------------
   class CatalogNotOpen(Exception):
      """ Catalog has not been open yet """

      def __init__(self, filepath):
         """!
            Exception constructor
            @param filepath file path of the catalog file 
         """
         self._filepath = filepath

      def __str__(self):
         return "SCatalog *** ERROR ***: Catalog: {0} is not open".format(self._filepath)

   # ------------------------------------------------------------------------------
   class CatalogFileNotFound(Exception):
      """!
         Exception thrown when a catalog file is not found on disk 
      """
      def __init__(self, filepath):
         """!
            Exception constructor
            @param filepath file path of the catalog file 
         """
         self._filepath = filepath

      def __str__(self):
         """! String representation of the exception object """
         return "SCatalog *** ERROR ***: file {0} no found".format(self._filepath)

   # ------------------------------------------------------------------------------
   class ColumnNotFound(Exception):
      """!
         Exception thrown when a named catalog column could not be found
      """
      def __init__(self, col_name):
         """!
            Exception constructor
            @param col_name name of the column 
         """
         self._col_name = col_name

      def __str__(self):
         """! String representation of the exception object """
         return "SCatalog *** ERROR ***: column {0} no found".format(self._col_name)


# ~~~~~~~~~~~
# Subclasses 
# ~~~~~~~~~~~

# --------------------------------------------------------------------------------------------------
class TextCatalog(BaseCatalog):
   """! 
      Catalog management in tabulated text (.ascii) format.
      @note The current implementation is based on the AstroAsciiData module from 
            http://www.stecf.org/software/PYTHONtools/astroasciidata/, but a reimplementation
            using numpy is envisaged.
   """

   def __init__(self, fullpath, sep_char='', comment_char='#', null_char='*'):
      """! 
         Create a catalog object and open/create the corresponding catalog file 
         @param fullpath full path of the catalog file to open or create
         @param sep_char character to use as column separator
         @param comment_char character to use to represent comments
         @param null_char character to use to represent null entries
      """
      BaseCatalog.__init__(self, fullpath)

      self._format = BaseCatalog.InputFormat.TabulatedText      # default input/output format

      self._sep_char = sep_char           # character to use as column delimiter
      self._comment_char = comment_char   # character used to represent comments
      self._null_char = null_char         # character to prepresent null entries

      self._cat_data = None
      #self.__table_data = None

   def __str__(self):
      info = self._cat_data.info()
      return info

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def sep_char(self):
      """ Return the character used as column delimiter """
      return self._sep_char

   @property
   def comment_char(self):
      """ Return the character used to represent comments"""
      return self._comment_char

   @property
   def null_char(self):
      """ Return the character to use to represent null entries """
      return self._null_char
   
   def get_nb_rows(self):
      """! 
         Get the number of rows in the catalog 
         @return number of rows
      """
      if not self._cat_data is None:
         return self._cat_data.nrows
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 

   def get_nb_cols(self):
      """! 
         Get the number of columns in the catalog 
         @return number of columns
      """
      if not self._cat_data is None:
         return self._cat_data.ncols
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 

   def get_col_names(self):
      """! 
         Get the list of column names in the catalog 
         @return list of column names
      """
      if not self._cat_data is None:
         return [c.colname for c in self._cat_data]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 


   # --- Setters  

   @sep_char.setter
   def sep_char(self, sep_char):  
      """! 
         Set the character used to as column separator 
         @param sep_char character used as column separator
      """
      self._sep_char = sep_char

   @comment_char.setter
   def comment_char(self, comment_char):  
      """! 
         Set the character used to represent comments 
         @param comment_char character for comments
      """
      self._comment_char = comment_char 

   @null_char.setter 
   def null_char(self, null_char):  
      """! 
         Set the character to prepresent null entries 
         @param null_char character for null entries   
      """
      self._null_char = null_char 



   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   # --- Opening, closing, creating

   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):
         # --- Open catalog file
         self._cat_data = asc.open(self.fullpath, 
                                   #delimiter=self.sep_char, 
                                   comment_char=self.comment_char, 
                                   null=self.null_char)

         # --- If a header is present, rename the columns based on that header
         self._rename_cols_from_header(self.get_header())   
            
      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)


   def create(self, nb_rows=0, nb_cols=0, header=None):
      """! 
         Create a new catalog with a text-tabulated format
         @param nb_rows initial number of rows
         @param nb_cols initial number of columns
         @param header string to use as header, if required   
      """
      if nb_rows == 0 or nb_cols == 0:
         # --- Create an empty catalog
         self._cat_data = asc.create(1, 1, null=self.null_char)
         del self._cat_data[0] 
      else:
         # --- Create a catalog with an initial number of rows and columns
         self._cat_data = asc.create(nb_cols, nb_rows, 
                                     null=self.null_char)
      if not header is None:
         self._cat_data.header.append(header)
      
      self._cat_data.writeto(self.fullpath, True, True)      


#   def create_from_numpy(self, matrix, col_names):
#      """!
#         Create a new catalog from a two-dimensional numpy array
#         @param matrix two-dimensional numpy array
#         @param col_names list of column names to use as header
#      """
#      [nbrows, nbcols] = matrix.shape
#      new_cat = self.create()
#      for icol in xrange(nbcols):
#         col_obj = TextCatalog.Column("Column{0}".format(icol), list(matrix[:,icol]))         
#         new_cat.append_col(col_obj)       
#      for col_name in col_names:
#         new_cat[col_names.index(col_name)].rename(col_name)

#      new_cat.writeto(self.fullpath, True, True)      


   def create_from_numpy(self, matrix, col_names):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
      """
      # --- Create catalog file
      np.savetxt(self.fullpath, matrix)
      self._cat_data = asc.open(self.fullpath, 
                                #delimiter=self.sep_char, 
                                comment_char=self.comment_char,
                                null=self.null_char)

      # --- Create header
      self._cat_data.header.reset()
      self._create_header_from_col_names(col_names, sep_char=" ")

      # --- Rename columns
      self._rename_cols_from_header(self._cat_data.header)

      # --- Update file on disk      
      self._cat_data.writeto(self.fullpath, False, True)      


   def close(self):
      if not self._cat_data is None:
         # Note: adding flush() here causes problems!
         del self._cat_data
         self._cat_data = None
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

   def save(self, write_header=True):
      if not self._cat_data is None:
         self._cat_data.writeto(self.fullpath, False, write_header)      
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

#   def flush(self):
#      """ Flush the internal catalog buffers """
#      if not self._cat_data is None:
#         self._cat_data.flush()
#      else:
#         raise BaseCatalog.CatalogNotOpen(self.fullpath)       



   # --- Getters

   def get_header(self, index=0):
      """! 
         Return the catalog header possibly (at a specific index) defined in the catalog 
         @param index zero-based index of the header to retrieve 
      """
      if not self._cat_data is None:
         if len(self._cat_data.header) > index:
            return self._cat_data.header[index].strip()
         else:
            return ""
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_headers(self):
      """! get the list of catalog headers possibly defined in the catalog """
      if not self._cat_data is None:
         if len(self._cat_data.header) > 0:
            return [h.strip() for h in self._cat_data.header]
         else:
            return []
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_info(self):
      """! Retrieve some information about the catalog """
      if not self._cat_data is None:
         return self._cat_data.info()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_col_index(self, col_name):
      """! 
         Get the column index from a column name 
         @param col_name column name
         @return column index (zero-based)
         @note Checks disabled for performance reasons
      """
      return self._cat_data.find(col_name)

   def get_col(self, col_index):
      """! Get a Column object from its index """
      if not self._cat_data is None:
         cat_col = self._cat_data[col_index]
         col_obj = TextCatalog.Column(cat_col.colname, 
                                      cat_col.get_format(), cat_col.get_colcomment(),
                                      cat_col.tonumpy())
         return col_obj
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_named_col(self, col_name):
      """! Get a Column object from its name """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            return self.get_col(col_index)
         else:
            raise BaseCatalog.ColumnNotFound(col_name)             
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_col_data(self, col_index):
      """! Return the data of a column from its index """
      if not self._cat_data is None:
         return self._cat_data[col_index].tonumpy()
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_named_col_data(self, col_name):
      """! Return the data of a column from its index """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            return self.get_col_data(col_index)
         else:
            raise BaseCatalog.ColumnNotFound(col_name)             
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_data(self):

      """! Convert the catalog to a two-dimensional, transposed numpy array """         
      if not self._cat_data is None:
         col_list = []
         for icol in xrange(self._cat_data.ncols):
            col_list.append(np.asarray([]))
      
         # Concatenate column data
         for icol in xrange(self._cat_data.ncols):
            col_data = np.asarray(list(self._cat_data[icol]))
            col_list[icol] = np.concatenate( (col_list[icol], col_data) )

#         for icol in xrange(self._cat_data.ncols):
#            col_list[icol] = np.asarray(list(self._cat_data[icol]))

         # Create a catalog Numpy matrix from the SE columns
         return np.asmatrix(col_list).transpose()
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_item(self, row, col):
      """! Retrieve an item at row @c row and column @c col
         @note Checks disabled for performance reasons
      """
      return self._cat_data[col][row]


   # --- Setters

   def set_item(self, row, col, value):
      """! 
         Set the value of an item at row @c row and column @c col 
         @param row row of the item to set
         @param col column of the item to set
         @param value value to assign to the item at (@c row, @c col)
         @note row and column indice checking disabled for performance reasons
      """
      self._cat_data[col][row] = value

   def set_header(self, header, index=0):  
      """! 
         Set the header at a specific index 
         @param header header string
         @param index zero-based header index   
      """
      if not self._cat_data is None:
         header_length = len(self._cat_data.header)
         if header_length == 0:
            # Append new header
            self._cat_data.header.append(header)
         else:   
            # Header already created
            if header_length > index:
               self._cat_data.header[index] = header
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def set_headers(self, headers):  
      """! 
         Specify a list of headers  
         @param headers list of header strings    
      """
      if not self._cat_data is None:
         del self._cat_data.header
         self._cat_data.header = asc.Header()
         for h in headers:
            self._cat_data.header.append(h)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def set_col_data(self, col_index, data):
      """! 
         Set the data of a column from its index 
         @param col_index index of the column whose data is to be retrieved
         @param data one-dimensional numpy array
      """
      if not self._cat_data is None:
         old_col = self._cat_data[col_index]
         #data = self.helper.fix_data_type(data)
         new_col = asc.AsciiColumn(list(data))    
         new_col.rename(old_col.colname)
         new_col.reformat(old_col.get_format())
         self._cat_data[col_index] = new_col
         self._cat_data[col_index].set_colcomment(old_col.get_colcomment())
         del old_col
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def set_named_col_data(self, col_name, data):
      """! 
         Set the data of a column from its index 
         @param col_name name of the column whose data is to be retrieved
         @param data one-dimensional numpy array
      """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            self.set_col_data(col_index, data)
         else:
            raise BaseCatalog.ColumnNotFound(col_name)
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def set_data(self, data):
      """! 
         Set the catalog data from a numpy two-dimensional matrix 
         @param data two-dimensional numpy array with the data to set
         @note This method is quite slow: maybe better to use create_from_numpy() 
      """
      if not self._cat_data is None:
         #data = self.helper.fix_data_type(data)
         nb_rows, nb_cols = data.shape
         for col_index in xrange(nb_cols):
#            for row_index in xrange(nb_rows):
#               self._cat_data[col_index][row_index] = 0.0
#               print col_index, type(data[row_index, col_index]), self._cat_data[col_index].get_type()
#               self._cat_data[col_index][row_index] = float(data[row_index, col_index])   

#            print "***", col_index, data[:,col_index]
            self.set_col_data(col_index, data[:,col_index])
     
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)
                  


   # --- Manipulating rows

   def insert_rows(self, start=0, nb_rows=1):
      """! 
         Insert @c nb_rows starting from index @c start
         @param start index where to insert
         @param nb_rows numer of rows to insert
         @note Checks disabled for performance reasons
      """
      self._cat_data.insert(nb_rows, start)
      
   def delete_rows(self, start=0, end=1):
      """! 
         Insert @c nb_rows starting from index @c start
         @param start index where to start deleting
         @param end the first row index not to be deleted
         @note Checks disabled for performance reasons
      """
      self._cat_data.delete(start, end)


   # --- Changing catalog structure

#   def append_col(self, column):
#      """!
#         Append a Column object after index @c index
#         @column an object derrived from BaseCatalog.Column
#         @index index zero-based index of the column to append 
#      """
#      if not self._cat_data is None:
#         nb_rows = self._cat_data.nrows
#         self._cat_data.append(column.name)  
#         self._cat_data.insert(len(column.data))
#         for irow in xrange(len(column.data)):
#            self._cat_data[column.name][irow] = column.data[irow]
#         print "==> APPEND COL:", type(column.data[0])

#         self._cat_data[column.name].rename(column.name)
#         self._create_header_from_col_names(self.get_col_names())
#         self._cat_data[column.name].reformat(column.format)
#         self._cat_data[column.name].set_colcomment(column.comment)
#      else:
#         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def append_col(self, column):
      """!
         Append a Column object after index @c index
         @param column an object derrived from BaseCatalog.Column
      """
      if not self._cat_data is None:
         nb_rows = self._cat_data.nrows
         self._cat_data.append(column.name)  
         self._cat_data.insert(len(column.data))
         self._cat_data.delete(0, nb_rows)

         col_obj = asc.AsciiColumn(list(column.data))
         col_obj.rename(column.name)
         self._cat_data[column.name] = col_obj
         self._create_header_from_col_names(self.get_col_names())
         self._cat_data[column.name].reformat(column.format)
         self._cat_data[column.name].set_colcomment(column.comment)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def remove_col(self, col_index):
      """!
         Delete a column from its index
         @param col_index index of the column to delete
      """   
      if not self._cat_data is None:
         del self._cat_data[col_index]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)  

   def remove_named_col(self, col_name):
      """!
         Delete a named column
         @param col_name name of the column to delete
      """   
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            self.remove_col(col_index)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def col_exists(self, col_name):
      """! 
         Tell whether a named column exists 
         @return True if column exists, False otherwise
      """
      if not self._cat_data is None:
         return self._cat_data.find(col_name) != -1      
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             
      
   # --- Transforming

   def to_plain(self):
      """! Convert the catalog to plain format """
      if not self._cat_data is None:
         self._cat_data.toplain()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             
      
   def to_SExtractor(self):
      """! Convert the catalog to SExtractor format """
      if not self._cat_data is None:
         self._cat_data.toSExtractor()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def to_FITS(self):
      """! Convert the catalog to FITS format """
      if not self._cat_data is None:
         self._cat_data.tofits()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def sort(self, col_index, descending=False):
      """! Sort the catalog data about column  @c col_index  
         @param col_index column index
         @param descending if True, sort in descending order, otherwise, in ascending order
      """
      if not self._cat_data is None:
         self._cat_data.sort(col_index, ordered=True, descending=descending)

#         data = self.get_data()
#         sorted_indice = data[:, col_index].argsort(kind="mergesort")               
#         sorted_data = np.take(data, sorted_indice, 0)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   def sort_by_col_name(self, col_name, descending=False):
      """! Sort the catalog data about column  @c col_name  
         @param col_name column name
         @param descending if True, sort in descending order, otherwise, in ascending order
      """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            self.sort(col_index, descending)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # --- Displaying
   def print_content(self):
      if not self._cat_data is None:
         print self._cat_data
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

            


   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------   
   def _rename_cols_from_header(self, header):
      """! Set column names based on header information if any """
      if len(header) > 0:
         head_string = self.get_header().rstrip('\r')
         head_names = re.sub(r'\s+', ' ', head_string).strip().split()
         if len(head_names) <= len(self.get_col_names()):
            for h in head_names:
               self._cat_data[head_names.index(h)].rename(h)

   # -----------------------------------------------------------------------------------------------
   def _create_header_from_col_names(self, col_names, sep_char=' '):
      """ 
         ! Create a header string from a list of comlumn names 
         @param col_names list of clumn names   
      """
      header = " "
      for col_name in col_names:
         header += (col_name + sep_char)
      self.set_header(header.rstrip(sep_char))


   # ~~~~~~~~~~~~~~
   # Nested Classes 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class Column(BaseCatalog.Column):
      """! Represents a column in the catalog """

      def __init__(self, name, format=None, comment=None, data=None):
         """!
            Create a Column object
            @param name name of the column 
            @param format Python format of the column 
            @param comment comment string
            @param data associated column data 
         """
         BaseCatalog.Column.__init__(self, name, format, comment, data)

         if data is None:
            data = []
         else:
            if not isinstance(data, list):
               data = list(data)
         self._cat_col = asc.AsciiColumn(data)  # column object (AsciiColumn class) 
         self._data = self._cat_col.tonumpy()   # column data (numpy array) 

         self._cat_col.rename(name)
         if not format is None:
            self._cat_col.reformat(format)

      def __str__(self):
         info = self._cat_col.info()
         info += "Nrows: {0}".format(len(self._data)) 
         return info

      # ~~~~~~~~~~~
      # Properties 
      # ~~~~~~~~~~~

      # --- Getters

#      @property
#      def format(self):
#         """! 
#            Get the format of the column 
#            @return format of the column
#         """
#         if self._format is None:
#            return self._cat_col.get_format()
#         

#      @property
#      def comment(self):
#         """! 
#            Get the comment string associated with the column 
#            @return comment comment string
#         """
#         if self._comment is None:
#            return self._cat_col.get_colcomment()

      # --- Setters




      # ~~~~~~~~~~~~~~
      # Public Methods 
      # ~~~~~~~~~~~~~~

      # --- Getters

      def get_nb_rows(self):
         """! Retrieve the number of rows of the column """
         return self._cat_col.get_nrows()   

      def get_info(self):
         """! Retrieve information about the column """
         return self._cat_col.info()   

      def get_type(self):
         """! 
            Get the data type of the column 
            @return Python data type of the column
          """
         return self._cat_col.get_type()

#      def get_format(self):
#         """! 
#            Get the expected format of the column
#            @return the Python format of the column
#         """
#         return self._cat_col._get_format()

#      def get_comment(self):
#         return self._cat_col.get_colcomment()


#      # --- Setters

#      def set_format(self, format):
#         """! 
#            Set the format of the column
#            @param format new format to set
#         """
#         self._cat_col.reformat(format)

#      def set_comment(self, comment):
#         """! 
#            Set the format of the column
#            @param commment new comment to set
#         """
#         self._cat_col.set_colcomment(comment)



# --------------------------------------------------------------------------------------------------
class SExCatalog(TextCatalog):

   """! 
      Catalogs management in SExtractor format 
      @note The current implementation is based on the AstroAsciiData module from 
            http://www.stecf.org/software/PYTHONtools/astroasciidata/, but a reimplementation
            using numpy is envisaged.
   """
   def __init__(self, fullpath, sep_char="", comment_char='#', null_char='*'):
      TextCatalog.__init__(self, fullpath, sep_char, comment_char, null_char)

      self._format = BaseCatalog.InputFormat.SExtractor      # default input/output format

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters



   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):
         # --- Open catalog file
         self._cat_data = asc.open(self.fullpath, 
                                    #delimiter=self.sep_char, 
                                    comment_char=self.comment_char, 
                                    null=self.null_char)

         # --- If a header is present, rename the columns based on that header
         self._rename_cols_from_headers(self.get_headers())   
            
      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   def create(nb_rows, nb_cols, header):
      """! 
         Create a new catalog with a SEXtractor format
         @param nb_rows initial number of rows
         @param nb_cols initial number of columns
         @param header string to use as header, if required   
      """
      if nb_rows == 0 or nb_cols == 0:
         # --- Create an empty catalog
         self._cat_data = asc.createSEx(1, 1, 
                                        null=self.null_char)
         del self._cat_data[0] 
      else:
         # --- Create a catalog with an initial number of rows and columns
         self._cat_data = asc.createSEx(nb_cols, nb_rows, null=self.null_char)

      if not header is None:
         self._cat_data.header.append(header)
      self._cat_data.writeto(self.fullpath, write_col_info=True, write_header=True)


   def create_from_numpy(self, matrix, col_names, col_comments=[], col_formats=[]):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
         @param col_comments list of column comment strings
         @param col_formats list of column formats
      """
      # --- Create ctalog file
      if len(matrix) > 0:
         np.savetxt(self.fullpath, matrix)
         self._cat_data = asc.open(self.fullpath, 
                                   #delimiter=self.sep_char, 
                                   comment_char=self.comment_char, 
                                   null=self.null_char)
      else:
         self._cat_data = asc.createSEx(len(col_names), 1, null=self.null_char)         

      # --- Rename columns
      for col_name in col_names:
         self._cat_data[col_names.index(col_name)].rename(col_name)

      # --- Set column comments
      if len(col_comments) >= len(col_names):   
         for col_idx in xrange(self._cat_data.ncols):
            self._cat_data[col_idx].set_colcomment(col_comments[col_idx]) 
         
      # --- Set column comments
      if len(col_formats) >= len(col_names):   
         for col_idx in xrange(self._cat_data.ncols):
            self._cat_data[col_idx].reformat(col_formats[col_idx])

      # --- Update file on disk      
      self._cat_data.writeto(self.fullpath, True, True)      


   def save(self, write_col_info=True, write_header=True):
      if not self._cat_data is None:
         if write_col_info:
            self._cat_data.header.reset()
         self._cat_data.writeto(self.fullpath, write_col_info, write_header)      
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

   def get_col_comments(self):
      if not self._cat_data is None:
         return [c.get_colcomment() for c in self._cat_data]   
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       


   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~

   def _rename_cols_from_headers(self, headers):
      """! Set column names based on header information if any """
      headers = self.get_headers()
      if len(headers) <= len(self.get_col_names()):
         for h in headers:
            self._cat_data[headers.index(h)].rename(h)

   def _create_header_from_col_names(self, col_names, sep_char=' '):
      """ 
         ! Create a header string from a list of comlumn names 
         @param col_names list of clumn names   
      """
      self._cat_data.header.reset()
      for col_name in col_names:
         self._cat_data.header.append(sep_char + col_name)


# --------------------------------------------------------------------------------------------------
class FITSCatalog(BaseCatalog):
   """! 
      Catalogs management in .FITS format 
      @note The implementation of this class is still work in progress.
      @note At present, .FITS table can only be read, not updated. Also features such like managing
            headers or sorting are not yet implenented.  
   """

   def __init__(self, fullpath, hdu_no=0, mem_map=False, open_mode=0):
      BaseCatalog.__init__(self, fullpath)

      self._format = BaseCatalog.InputFormat.FITS     # default input/output format
      self._hdu_no = hdu_no                           # HDU number of the underlying .FITS table
      self._mem_map = mem_map                         # True to use memory I/O mapping
      self._open_mode = FITSCatalog.OpenMode.ReadOnly # opening mode (see FITSCatalog.OpenMode)

   def __str__(self):

      if self._cat_data is not None:
         info = str(self._cat_data.fileinfo(self.hdu_no))
         info += str(self._cat_data.info())
      else:            
         info = str(pyfits.info(self.fullpath))
      return info

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   # --- Getters

   @property
   def hdu_no(self):
      """! 
         Get the HDU index of the table
         @return HDU index of the table
      """
      return self._hdu_no

   @property
   def mem_map(self):
      """! 
         Tell whether memory mapped I/O are enabled
         @return True if memory mapped I/O are enabled
      """
      return self._mem_map

   @property
   def open_mode(self):
      """! 
         Return the catalog opening mode 
         @return a constant of type FITSCatalog.OpenMode
         @see FITSCatalog.OpenMode
      """
      return self._open_mode

   # --- Setters

   @hdu_no.setter
   def hdu_no(self, hdu_no):
      """! 
         Set the HDU index of the table
         @param hdu_no HDU index of the table
      """
      self._hdu_no = hdu_no

   @mem_map.setter
   def mem_map(self, mem_map):
      """! 
         Enable/Disable memory mapped I/O 
         @param mem_map True to enable memory mapped I/O, False otherwise
      """
      self._mem_map = mem_map

   @open_mode.setter
   def mem_map(self, open_mode):
      """! 
         Set the catalog opening mode
         @param open_mode a constant of type FITSCatalog.OpenMode
      """
      self._open_mode = open_mode

   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_nb_rows(self):
      """! 
         Get the number of rows in the catalog 
         @return number of rows or 0 if None
      """
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 
      return self._cat_data[self.hdu_no].header.get("NAXIS2", 0)

   # -----------------------------------------------------------------------------------------------
   def get_nb_cols(self):
      """! 
         Get the number of columns in the catalog 
         @return number of columns or 0 if None
      """
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 
      return self._cat_data[self.hdu_no].header.get("TFIELDS", 0)
   

   def get_col_names(self):
      """! 
         Get the list of column names in the catalog 
         @return list of column names
      """
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath) 
      else:
         return [attr[1] for attr in self._cat_data[self.hdu_no].header.cards \
                         if attr.keyword.startswith("TTYPE")]

   # --- Opening, closing, creating

   #------------------------------------------------------------------------------------------------
   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):
         # --- Open catalog file
         if self.open_mode == FITSCatalog.OpenMode.ReadOnly:  
            self._cat_data = pyfits.open(self.fullpath, "readonly", memmap=self.mem_map)
         else:
            raise BaseCatalog.FeatureNotImplemented("open() for update")
            #self._cat_data = pyfits.open(self.fullpath, "update", memmap=self.mem_map)

      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   #------------------------------------------------------------------------------------------------
   def create(self, nb_rows=0, nb_cols=0, header=None):
      """! 
         Create a new catalog with a text-tabulated format
         @param nb_rows initial number of rows
         @param nb_cols initial number of columns
         @param header string to use as header, if required   
      """
      raise BaseCatalog.FeatureNotImplemented("create()")      

   #------------------------------------------------------------------------------------------------
   def create_from_numpy(self, matrix, col_names):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
      """
      raise BaseCatalog.FeatureNotImplemented("create()")      

#      # --- Create ctalog file
#      np.savetxt(self.fullpath, matrix)
#      self._cat_data = asc.open(self.fullpath, 
#                                #delimiter=self.sep_char, 
#                                comment_char=self.comment_char,
#                                null=self.null_char)

#      # --- Create header
#      self._cat_data.header.reset()
#      self._create_header_from_col_names(col_names, sep_char=" ")

#      # --- Rename columns
#      self._rename_cols_from_header(self._cat_data.header)

#      # --- Update file on disk      
#      self._cat_data.writeto(self.fullpath, False, True)      


   def close(self):
      if not self._cat_data is None:
         if self.open_mode == FITSCatalog.OpenMode.ReadWrite:
            self._cat_data.flush()
         del self._cat_data
         self._cat_data = None
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

   def save(self, write_header=True):
      raise BaseCatalog.FeatureNotImplemented("save()")      
#      if not self._cat_data is None:
#         self._cat_data.writeto(self.fullpath, False, write_header)      
#      else:
#         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

   def flush(self):
      """ Flush the internal catalog buffers """
      if not self._cat_data is None: 
         if self.open_mode == FITSCatalog.OpenMode.ReadWrite:
            self._cat_data.flush()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       


# TODO:
#   def get_header(self, index=0):
#      """! 
#         Return the catalog header possibly (at a specific index) defined in the catalog 
#         @param index zero-based index of the header to retrieve 
#      """
#      if not self._cat_data is None:
#         if len(self._cat_data.header) > index:
#            return self._cat_data.header[index].strip()
#         else:
#            return ""
#      else:
#         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

#   def get_headers(self):
#      """! get the list of catalog headers possibly defined in the catalog """
#      if not self._cat_data is None:
#         if len(self._cat_data.header) > 0:
#            return [h.strip() for h in self._cat_data.header]
#         else:
#            return []
#      else:
#         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_info(self):
      """! 
         Retrieve some information about the catalog 
         @return a dictionary with detailed information
         @note See the pyfits documentation of the info() function for the details
      """
      if not self._cat_data is None:
         return self._cat_data.fileinfo(self.hdu_no)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)       

   def col_exists(self, col_name):
      """! 
         Tell whether a named column exists 
         @return True if column exists, False otherwise
      """
      return col_name in self.get_col_names()            

   def get_col_index(self, col_name):
      """! 
         Get the column index from a column name 
         @param col_name column name
         @return column index (zero-based)
         @note Checks disabled for performance reasons
      """
      col_names = self.get_col_names()
      if not col_name in col_names:
         raise BaseCatalog.ColumnNotFound(col_name)             
      return col_names.index(col_name)


   def get_col_data(self, col_index):
      """! Return the data of a column from its index """
      if not self._cat_data is None:
         return self._cat_data[self.hdu_no].data.field(col_index)
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_named_col_data(self, col_name):
      """! Return the data of a column from its index """
      if not self._cat_data is None:
         return self._cat_data[self.hdu_no].data[col_name]
      else:   
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             

   def get_data(self):
      """! Convert the catalog into a two-dimensional, data using numpy convention """         
      if not self._cat_data is None:
         cat_col_list = []
         ncols = self.get_nb_cols()
         for icol in xrange(ncols):
            cat_col_list.append(np.asarray([]))
         for icol in xrange(ncols):
            cat_col_list[icol] = np.concatenate( (cat_col_list[icol], 
                                                 self._cat_data[self.hdu_no].data.field(icol)) ) 
         return np.asmatrix(cat_col_list).transpose().squeeze()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)             


   # -----------------------------------------------------------------------------------------------
   class OpenMode:
      """! Supported input catalog open mode """
      (ReadOnly, ReadWrite) = (0,1)



# --------------------------------------------------------------------------------------------------
def main():
   
   catalog = TextCatalog("./input/text_cat.txt")
   catalog.skip_rows = 2

   catalog.read()   


if __name__ == "__main__":
    main()


# -- EOF scatalog.py
