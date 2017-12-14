"""!
   @package scatalog.scatalog Simple catalog management
   @author Marc Gentile
   @file scatalog.py
   Simple catalog management for astronomical data
"""

# -- Python imports
import os
import string
import re
import numpy as np
import asciidata as asc             # for text catalogs
from astropy.io import fits         # for FITS catalogs (New)
import astromatic_wrapper.utils as awu     # for FITS_LDAC catalogs

# -- Internal imports
from scatalog_helper import *       # Helper methods
# from numpy import ndarray
# from astropy.io.fits.fitsrec import FITS_record


# --------------------------------------------------------------------------------------------------
class BaseCatalog(object):

   """!
      Base catalog management class
      @note This class is not meant to be used directly
   """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, fullpath):
      """!  BaseCatalog class constructor """

      # --- Public members
      self._directory, self._filename = os.path.split(fullpath)  # catalog file path
      self._format = BaseCatalog.InputFormat.Undefined  # input/output format
      self._helper = Helper()

      # --- Private members
      self._cat_data = None  # catalog internal data (e.g. AsciiData class)


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

   # -----------------------------------------------------------------------------------------------
   def get_nb_rows(self):
      """!
         Get the number of rows in the catalog
         @return number of rows
      """
      raise BaseCatalog.FeatureNotImplemented("get_nb_rows()")

   # -----------------------------------------------------------------------------------------------
   def get_nb_cols(self):
      """!
         Get the number of columns in the catalog
         @return number of columns
      """
      raise BaseCatalog.FeatureNotImplemented("get_nb_cols()")

   # -----------------------------------------------------------------------------------------------
   def get_col_names(self):
      """!
         Get the list of column names in the catalog
         @return list of column names
      """
      raise BaseCatalog.FeatureNotImplemented("get_col_names()")

   # -----------------------------------------------------------------------------------------------
   def get_col_formats(self):
      """!
          Get the list of column formats in the order of columns
      """
      raise BaseCatalog.FeatureNotImplemented("get_col_names()")

#    # -----------------------------------------------------------------------------------------------
#    def append_col(self, column):
#       """!
#          Append a Column object
#          @param column an object derived from BaseCatalog.Column
#       """
#       raise BaseCatalog.FeatureNotImplemented("append_col()")


   # -----------------------------------------------------------------------------------------------
   def add_col(self, col_name, col_format=None, col_comment=None, col_data=None):
      """!
         Add a Column to the catalog
         @param col_name column name
         @param col_format column python format
         @param col_comment column comment (ignored in for class)
         @param col_data column data as a numpy array
      """
      raise BaseCatalog.FeatureNotImplemented("add_col()")


   # ~~~~~~~~~~~~~
   # Inner Classes
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class InputFormat:
      """! Supported input catalog formats"""
      (Undefined, TabulatedText, SExtractor, FITS, FITS_LDAC) = (0, 1, 2, 4, 5)


   # -----------------------------------------------------------------------------------------------
   class OpenMode:
      """! Supported input catalog open mode """
      (ReadOnly, ReadWrite, Append) = ("readonly", "update", "append")


   # -----------------------------------------------------------------------------------------------
   class Column(object):
      """! Represents a column in the catalog """

      def __init__(self):

         self._cat_col = None
         self._helper = Helper()

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
         raise BaseCatalog.FeatureNotImplemented("column.name")

      @property
      def format(self):
         """!
            Get the format of the column
            @return format of the column
         """
         raise BaseCatalog.FeatureNotImplemented("column.format")

      @property
      def data(self):
         """!
            Get the data associated with the column
            @return data data associated with the column
         """
         raise BaseCatalog.FeatureNotImplemented("column.data")

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

      # --------------------------------------------------------------------------------------------
      def get_nb_rows(self):
         """! Retrieve the number of rows of the column """
         raise BaseCatalog.FeatureNotImplemented("get_nb_rows()")

      # --------------------------------------------------------------------------------------------
      def get_info(self):
         """! Retrieve information about the column """
         raise BaseCatalog.FeatureNotImplemented("get_info()")

      # --------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   class CatalogNotCreated(Exception):
      """ Catalog could not be created """

      def __init__(self, filepath):
         """!
            Exception constructor
            @param filepath file path of the catalog file
         """
         self._filepath = filepath

      def __str__(self):
         return "SCatalog *** ERROR ***: Catalog: {0} could not be created".format(self._filepath)

   # -----------------------------------------------------------------------------------------------
   class OpenModeNotSupported(Exception):
      """ Opening mode not supported by this version """

      def __init__(self, filepath, open_mode):
         """!
            Exception constructor
            @param filepath file path of the catalog file
         """
         self._filepath = filepath
         self._open_mode = open_mode

      def __str__(self):
         return "SCatalog *** ERROR ***: Catalog: {0} Open Mode {1} not supported".format(
                                                               self._filepath, self._open_mode)

   #------------------------------------------------------------------------------------------------
   # ADDED
   class OpenModeConflict(Exception):
       """ Opening mode is in conflict with an action """

       def __init__(self, open_mode, open_mode_needed):
           """!
              Exception constructor
              @param open_mode_used
              @param open_mode_needed
           """
           self._open_mode = open_mode
           self._open_mode_needed = open_mode_needed

       def __str__(self):
           return "SCatalog *** ERROR ***: Catalog has to be open as : {0} , Mode used : {1}".format(self._open_mode_needed, self._open_mode)

# ~~~~~~~~~~~
# Subclasses
# ~~~~~~~~~~~

# --------------------------------------------------------------------------------------------------
class TextCatalog(BaseCatalog):
   """!
      Catalog management in tabulated text (ASCII) format.
      @note The current implementation is based on the AstroAsciiData module from
            http://www.stecf.org/software/PYTHONtools/astroasciidata/
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

      self._format = BaseCatalog.InputFormat.TabulatedText  # default input/output format

      self._sep_char = sep_char  # character to use as column delimiter
      self._comment_char = comment_char  # character used to represent comments
      self._null_char = null_char  # character to represent null entries

   def __str__(self):
      if not self._cat_data is None:
         info = self._cat_data.info()
      else:
         info = "No information"
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
         Set the character to represent null entries
         @param null_char character for null entries
      """
      self._null_char = null_char


   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_nb_rows(self):
      """!
         Get the number of rows in the catalog
         @return number of rows
         @note overridden from scatalog.BaseCatalog
      """
      if not self._cat_data is None:
         return self._cat_data.nrows
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_nb_cols(self):
      """!
         Get the number of columns in the catalog
         @return number of columns
         @note overridden from scatalog.BaseCatalog
      """
      if not self._cat_data is None:
         return self._cat_data.ncols
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_names(self):
      """!
         Get the list of column names in the catalog
         @return list of column names
         @note overridden from scatalog.BaseCatalog
      """
      if not self._cat_data is None:
         return [c.colname for c in self._cat_data]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_formats(self):
      """!
         Get the list of column formats in the order of columns
         @return list of column formats
         @note overridden from scatalog.BaseCatalog
      """
      if not self._cat_data is None:
         return [c.get_format() for c in self._cat_data]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


   # -----------------------------------------------------------------------------------------------
   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):
         # --- Open catalog file
         self._cat_data = asc.open(self.fullpath,
                                   comment_char=self.comment_char,
                                   null=self.null_char)

         # --- If a header is present, rename the columns based on that header
         self._rename_cols_from_header(self.get_header())

      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   def create_from_numpy(self, matrix, col_names):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
      """
      # --- Create catalog file
      np.savetxt(self.fullpath, matrix)
      self._cat_data = asc.open(self.fullpath,
                                # delimiter=self.sep_char,
                                comment_char=self.comment_char,
                                null=self.null_char)

      # --- Create header
      self._cat_data.header.reset()
      self._create_header_from_col_names(col_names, sep_char=" ")

      # --- Rename columns
      self._rename_cols_from_header(self._cat_data.header)

      # --- Update file on disk
      self._cat_data.writeto(self.fullpath, False, True)


   # -----------------------------------------------------------------------------------------------
   def close(self):
      """! Close the catalog file """
      if not self._cat_data is None:
         # Note: adding flush() here causes problems!
         del self._cat_data
         self._cat_data = None
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def save(self, write_header=True):
      """! Write the catalog data to disk """
      if not self._cat_data is None:
         self._cat_data.writeto(self.fullpath, False, write_header)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_header(self, index=0):
      """!
         Return the catalog header possibly (at a specific index) defined in the catalog
         @param index zero-based index of the header to retrieve
         @return header
      """
      if not self._cat_data is None:
         if len(self._cat_data.header) > index:
            return self._cat_data.header[index].strip()
         else:
            return ""
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_headers(self):
      """!
         get the list of catalog headers possibly defined in the catalog
         @return list of headers
      """
      if not self._cat_data is None:
         if len(self._cat_data.header) > 0:
            return [h.strip() for h in self._cat_data.header]
         else:
            return []
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_info(self):
      """!
         Retrieve some information about the catalog
         @return catalog information
      """
      if not self._cat_data is None:
         return self._cat_data.info()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_index(self, col_name):
      """!
         Get the column index from a column name
         @param col_name column name
         @return column index (zero-based)
         @note Checks disabled for performance reasons
      """
      return self._cat_data.find(col_name)

   # -----------------------------------------------------------------------------------------------
   def get_col(self, col_index):
      """!
         Get a Column object from its index
         @param col_index index of column (zero-based)
         @return column object for this class
      """
      if not self._cat_data is None:
         cat_col = self._cat_data[col_index]
         col_obj = TextCatalog.Column(cat_col.colname,
                                      cat_col.get_format(), cat_col.get_colcomment(),
                                      cat_col.tonumpy())
         return col_obj
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_named_col(self, col_name):
      """!
         Get a Column object from its name
         @param col_name column name
         @return column object for this class
      """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            return self.get_col(col_index)
         else:
            raise BaseCatalog.ColumnNotFound(col_name)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_data(self, col_index):
      """!
         Return the data of a column from its index
         @param col_index index of column (zero-based)
         @return column data as a numpy array
      """
      if not self._cat_data is None:
         return self._cat_data[col_index].tonumpy()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_named_col_data(self, col_name):
      """!
         Return the data of a column from its name
         @param col_name column name
         @return column data as a numpy array
      """
      if not self._cat_data is None:
         col_index = self._cat_data.find(col_name)
         if col_index != -1:
            return self.get_col_data(col_index)
         else:
            raise BaseCatalog.ColumnNotFound(col_name)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_data(self):
      """!
         Convert the catalog to a two-dimensional, transposed numpy array
         @return transposed numpy array
      """
      if not self._cat_data is None:
         col_list = []
         for icol in xrange(self._cat_data.ncols):
            col_list.append(np.asarray([]))

         # Concatenate column data
         for icol in xrange(self._cat_data.ncols):
            col_data = np.asarray(list(self._cat_data[icol]))
            col_list[icol] = np.concatenate((col_list[icol], col_data))

#         for icol in xrange(self._cat_data.ncols):
#            col_list[icol] = np.asarray(list(self._cat_data[icol]))

         # Create a catalog Numpy matrix from the SE columns
         return np.asmatrix(col_list).transpose()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_item(self, row, col):
      """! Retrieve an item at row @c row and column @c col
         @param row a row index in the catalog
         @param col a column index in the catalog
         @return item value
         @note Checks disabled for performance reasons
      """
      return self._cat_data[col][row]

   # -----------------------------------------------------------------------------------------------
   def set_item(self, row, col, value):
      """!
         Set the value of an item at row @c row and column @c col
         @param row row of the item to set
         @param col column of the item to set
         @param value value to assign to the item at (@c row, @c col)
         @note row and column indice checking disabled for performance reasons
      """
      self._cat_data[col][row] = value

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   def set_col_data(self, col_index, data):
      """!
         Set the data of a column from its index
         @param col_index index of the column whose data is to be retrieved
         @param data one-dimensional numpy array
      """
      if not self._cat_data is None:
         old_col = self._cat_data[col_index]
         # data = self.helper.fix_data_type(data)
         new_col = asc.AsciiColumn(list(data))
         new_col.rename(old_col.colname)
         new_col.reformat(old_col.get_format())
         self._cat_data[col_index] = new_col
         self._cat_data[col_index].set_colcomment(old_col.get_colcomment())
         del old_col
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   def set_data(self, data):
      """!
         Set the catalog data from a numpy two-dimensional matrix
         @param data two-dimensional numpy array with the data to set
         @note This method is quite slow: maybe better to use create_from_numpy()
      """
      if not self._cat_data is None:
         # data = self.helper.fix_data_type(data)
         _, nb_cols = data.shape
         for col_index in xrange(nb_cols):
#            for row_index in xrange(nb_rows):
#               self._cat_data[col_index][row_index] = 0.0
#               print col_index, type(data[row_index, col_index]), self._cat_data[col_index].get_type()
#               self._cat_data[col_index][row_index] = float(data[row_index, col_index])

#            print "***", col_index, data[:,col_index]
            self.set_col_data(col_index, data[:, col_index])

      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def insert_rows(self, start=0, nb_rows=1):
      """!
         Insert @c nb_rows starting from index @c start
         @param start index where to insert
         @param nb_rows numer of rows to insert
         @note Checks disabled for performance reasons
      """
      self._cat_data.insert(nb_rows, start)

   # -----------------------------------------------------------------------------------------------
   def delete_rows(self, start=0, end=1):
      """!
         Insert @c nb_rows starting from index @c start
         @param start index where to start deleting
         @param end the first row index not to be deleted
         @note Checks disabled for performance reasons
      """
      self._cat_data.delete(start, end)

   # -----------------------------------------------------------------------------------------------
   def add_col(self, col_name, col_format=None, col_comment=None, col_data=None):
      """!
         Add a Column to the catalog
         @param col_name column name
         @param col_format column python format
         @param col_comment column comment (ignored in for class)
         @param col_data column data as a numpy array
      """

      if not self.col_exists(col_name):
         column = SExCatalog.Column(name=col_name, format=col_format, comment=col_comment,
                                                                      data=col_data)
         self._append_col(column)


   # -----------------------------------------------------------------------------------------------
   def remove_col(self, col_index):
      """!
         Delete a column from its index
         @param col_index index of the column to delete
      """
      if not self._cat_data is None:
         del self._cat_data[col_index]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   def col_exists(self, col_name):
      """!
         Tell whether a named column exists
         @return True if column exists, False otherwise
      """
      if not self._cat_data is None:
         return self._cat_data.find(col_name) != -1
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def to_plain(self):
      """! Convert the catalog to plain format """
      if not self._cat_data is None:
         self._cat_data.toplain()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def to_SExtractor(self):
      """! Convert the catalog to SExtractor format """
      if not self._cat_data is None:
         self._cat_data.toSExtractor()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def to_FITS(self):
      """! Convert the catalog to FITS format """
      if not self._cat_data is None:
         self._cat_data.tofits()
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
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

   # -----------------------------------------------------------------------------------------------
   def print_content(self):
      if not self._cat_data is None:
         print self._cat_data
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _append_col(self, column):
      """!
         Append a Column object
         @param column an object derived from BaseCatalog.Column
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

      # --------------------------------------------------------------------------------------------
      def __init__(self, name, format=None, comment=None, data=None):
         """!
            Create a Column object
            @param name name of the column
            @param format Python format of the column
            @param comment comment string
            @param data associated column data
         """
         BaseCatalog.Column.__init__(self)

         if data is None:
            data = []
         else:
            if not isinstance(data, list):
               data = list(data)
         self._cat_col = asc.AsciiColumn(data)  # column object (AsciiColumn class)
         self._comment = comment  # comment string

         self._cat_col.rename(name)
         if not format is None:
            self._cat_col.reformat(format)

      # --------------------------------------------------------------------------------------------
      def __str__(self):
         info = self._cat_col.info()
         info += "Nb rows: {0}".format(self._cat_col.get_nrows())
         return info

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
         return self._cat_col.colname

      @property
      def format(self):
         """!
            Get the format of the column
            @return format of the column
            @note overridden from scatalog.BaseColumn
         """
         return self._cat_col.get_format()

      @property
      def data(self):
         """!
            Get the data associated with the column
            @return data data associated with the column
         """
         return self._cat_col.tonumpy()

      @property
      def comment(self):
         """!
            Get the comment string associated with the column
            @return comment comment string
         """
         return self._cat_col.get_colcomment()


      # --- Setters

      @name.setter
      def name(self, name):
         """!
            Set the name of the column
            @param name string
         """
         self._cat_col.rename(name)

      @format.setter
      def format(self, format):
         """!
            Set the format of the column
            @param format string
         """
         self._cat_col.reformat(format)

      @comment.setter
      def comment(self, comment):
         """!
            Set the comment string of the column
            @param comment comment string
         """
         self._cat_col.set_colcomment(comment)


      # ~~~~~~~~~~~~~~
      # Public Methods
      # ~~~~~~~~~~~~~~

      # --- Getters

      # --------------------------------------------------------------------------------------------
      def get_nb_rows(self):
         """! Retrieve the number of rows of the column """
         return self._cat_col.get_nrows()

      # --------------------------------------------------------------------------------------------
      def get_info(self):
         """! Retrieve information about the column """
         return self._cat_col.info()

      # --------------------------------------------------------------------------------------------
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
   # -----------------------------------------------------------------------------------------------
   def __init__(self, fullpath, sep_char="", comment_char='#', null_char='*'):
      TextCatalog.__init__(self, fullpath, sep_char, comment_char, null_char)

      self._format = BaseCatalog.InputFormat.SExtractor  # default input/output format

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters



   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):
         # --- Open catalog file
         self._cat_data = asc.open(self.fullpath,
                                    # delimiter=self.sep_char,
                                    comment_char=self.comment_char,
                                    null=self.null_char)

         # --- If a header is present, rename the columns based on that header
         self._rename_cols_from_headers(self.get_headers())

      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def create(self, nb_rows, nb_cols, header):
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


   # -----------------------------------------------------------------------------------------------
   def create_from_numpy(self, matrix, col_names, col_comments=[], col_formats=[]):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
         @param col_comments list of column comment strings
         @param col_formats list of column formats
      """
      # --- Create catalog file
      if len(matrix) > 0:
         np.savetxt(self.fullpath, matrix)
         self._cat_data = asc.open(self.fullpath,
                                   # delimiter=self.sep_char,
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
         try:
            for col_idx in xrange(self._cat_data.ncols):
               self._cat_data[col_idx].set_type(np.float)
               self._cat_data[col_idx].reformat(col_formats[col_idx])
         except Exception:
            pass

      # --- Update file on disk
      self._cat_data.writeto(self.fullpath, True, True)


   # -----------------------------------------------------------------------------------------------
   def save(self, write_col_info=True, write_header=True):
      if not self._cat_data is None:
         if write_col_info:
            self._cat_data.header.reset()
         self._cat_data.writeto(self.fullpath, write_col_info, write_header)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_comments(self):
      if not self._cat_data is None:
         return [c.get_colcomment() for c in self._cat_data]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _rename_cols_from_headers(self, headers):
      """! Set column names based on header information if any """
      headers = self.get_headers()
      if len(headers) <= len(self.get_col_names()):
         for h in headers:
            self._cat_data[headers.index(h)].rename(h)

   # -----------------------------------------------------------------------------------------------
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
            headers or sorting are not implemented.
   """

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def __init__(self, fullpath, hdu_no=None,
                      open_mode=BaseCatalog.OpenMode.ReadOnly, memmap=False, SEx_catalog=False):
      BaseCatalog.__init__(self, fullpath)

      self._format = BaseCatalog.InputFormat.FITS  # default input/output format

      self._open_mode  = open_mode      # opening mode (see FITSCatalog.OpenMode)
      self._SEx_catalog = SEx_catalog   # Work with SExtractor fits format or not
      if hdu_no is None:                # HDU number of the underlying .FITS table
          if SEx_catalog:               # Default is 1 (or 2 if you are using )
              self._hdu_no = 2
          else:
              self._hdu_no = 1
      else:
          self._hdu_no = hdu_no
      self._use_memmap = memmap         # use memory mapping or not

   # -----------------------------------------------------------------------------------------------
   def __str__(self):

      if self._cat_data is not None:
         info = "{0}".format(self.get_info())
      else:
         info = "No information"
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
   def open_mode(self):
      """!
         Return the catalog opening mode
         @return a constant of type FITSCatalog.OpenMode
         @see FITSCatalog.OpenMode
      """
      return self._open_mode

   @property
   def use_memmap(self):
      """!
         If True, use memory mapping
         @return whether to use memory mapping of not
      """
      return self._use_memmap

   # --- Setters

   @hdu_no.setter
   def hdu_no(self, hdu_no):
      """!
         Set the HDU index of the table
         @param hdu_no HDU index of the table
      """
      self._hdu_no = hdu_no

   @open_mode.setter
   def open_mode(self, open_mode):
      """!
         Set the catalog opening mode
         @param open_mode a constant of type FITSCatalog.OpenMode
      """
      self._open_mode = open_mode

   @use_memmap.setter
   def use_memmap(self, use_memmap):
      """!
         Use memory mapping if requested
         @param use_memmap wehether to use memory mapping of not
      """
      self._use_memmap = use_memmap

   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   #------------------------------------------------------------------------------------------------
   def open(self):

      """! Open an existing catalog """
      if self.helper.file_exists(self.fullpath):

         # --- Open catalog file
         self._cat_data = fits.open(self.fullpath, mode=self.open_mode, memmap=self.use_memmap)
      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   #------------------------------------------------------------------------------------------------
   # MODIFIED
   def create(self, ext_name=None, overwrite=True, s_hdu=True, sex_cat_path=None):
      """!
         Create an empty catalog with a FITS format
         @param ext_name extension name or number
         @param overwrite if the file already exist overwrite it
         @param secondary_hdu if True add a secondary hdu
      """

      # --- Create catalog
      primary_hdu = fits.PrimaryHDU()
      if self._SEx_catalog:
          if sex_cat_path is not None:
              if self.helper.file_exists(sex_cat_path):
                  sex_cat=FITSCatalog(sex_cat_path, hdu_no=1)
                  sex_cat.open()
                  secondary_hdu=sex_cat._cat_data[1]
                  self._cat_data = fits.HDUList([primary_hdu, secondary_hdu])
                  self._cat_data.writeto(self.fullpath, overwrite=overwrite)
                  sex_cat.close()
                  del(sex_cat)
              else:
                  raise BaseCatalog.CatalogFileNotFound(sex_cat_path)
          else:
              raise ValueError('sex_cat_path need to be provided to create SEXtractor catalog like')
      elif s_hdu:
          secondary_hdu = fits.BinTableHDU(data=None, header=None, name=ext_name)
          self._cat_data = fits.HDUList([primary_hdu, secondary_hdu])
          self._cat_data.writeto(self.fullpath, overwrite=overwrite)
      else:
          self._cat_data = fits.HDUList([primary_hdu])
          self._cat_data.writeto(self.fullpath, overwrite=overwrite)

   #------------------------------------------------------------------------------------------------
   # ADDED
   def copy_hdu(self, fits_file=None, hdu_no=None, hdu_name=None):
       """!
           Copy an HDU from fits
           @param fits_file scatalog where the hdu is
           @param hdu_no hdu to copy default=value use when fits_file was open

           NOTE : the scatalog to copy need to be open
       """

       if fits_file is None:
           fits_file=self

       if self._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(self.fullpath)
       if fits_file._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(fits_file.fullpath)

       if hdu_no is None:
           hdu_no=fits_file.hdu_no

       if hdu_name is None:
           hdu_name=fits_file._cat_data[hdu_no].name

       self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no=hdu_no),name=hdu_name))

   #------------------------------------------------------------------------------------------------
   #    ADDED
   def apply_mask(self, fits_file=None, hdu_no=None, mask=None , hdu_name=None):
       """!
           Apply a mask on data on data of a specified HDU
           @param fits file where the data to mask are
           @param hdu_no on which hdu are data to mask
           @param mask array of boolean or array of index
           @param hdu_name name to give to the new extension with data mask
                   (hdu_name=None -> same as the hdu_no.name)

           NOTE : this will create a new hdu with data[mask] in it
       """
       if fits_file is None:
           fits_file=self
       if self._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(self.fullpath)
       if fits_file._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(fits_file.fullpath)
       if maks is None:
           raise ValueError('Mask not provided')
       if type(mask) is not np.ndarray:
           raise TypeError('Mask need to be a numpy.ndarray')
       if hdu_no is None:
           hdu_no=fits_file.hdu_no

       if hdu_name is None:
           hdu_name=fits_file._cat_data[hdu_no].name

       if mask.dtype==bool:
           mask=np.where(mask==True)
           self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no)[:][mask],name=hdu_name))
       elif mask.dtype==int:
           self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no)[:][mask],name=hdu_name))
       else:
           raise TypeError('Mask type has to be int or bool')

   #------------------------------------------------------------------------------------------------
   # ADDED
   def save_as_fits(self, data=None, names=None, ext_name=None, sex_cat_path=None, image=False, overwrite=False):
       """!
            Save data into an already existing fits or a new one.
            Save data from dict, list, numpy.ndarray, numpy.recarray or astropy.io.fits.fitsrec.FITS_rec (data format in an astropy fits file)
            When creating a new fits to store BinTable data it create a PrimaryHDU.
            When creating a new fits to store Image there is no PrimaryHDU.
            You can create a SExtractor format fits by specifying a SExtractor catalog from where data come from.
            @param data data to store
            @param names names of the diferent column (not needed for dict and rec format)
            @param ext_name name of the HDU where data are stored (DEFAULT = NEW)
            @param sex_cat_path path of the already existing SExtractor catalog to mimic
            @param image if True create a fits image
            @param overwrite only used when creating an image fits

            NOTE : to create a SExtractor like fits you need to specify SEx_catalog=True when declaring the FITSCatalog object.
       """


       if data is None:
           raise ValueError('Data not provided')

       if not image:
           if type(data) is dict:
               r=np.rec.fromarrays([data[i] for i in data.keys()], names=data.keys())
               self._save_from_recarray(r,ext_name, sex_cat_path)

           elif type(data) is np.recarray:
               self._save_from_recarray(data,ext_name, sex_cat_path)

           elif type(data) is fits.fitsrec.FITS_rec:
               self._save_from_recarray(data,ext_name, sex_cat_path)

           elif type(data) is np.ndarray:
               if names is not None:
                   if (data.ndim<=2):
                       if data.ndim==1:
                           data=np.array([data])
                   else:
                       raise ValueError('Data dimension > 2')
                   r=np.rec.fromarrays([data[i] for i in range(data.shape[0])], names=names)
                   self._save_from_recarray(r,ext_name, sex_cat_path)
               else:
                   raise ValueError('Names not provided')

           elif type(data) is list:
               data=np.asarray(data)
               if (data.ndim<=2):
                   if data.ndim==1:
                       data=np.array([data])
               else:
                   raise ValueError('Data dimension > 2')
               r=np.rec.fromarrays([data[i] for i in range(data.shape[0])], names=names)
               self._save_from_recarray(r,ext_name, sex_cat_path)
       else:
           if type(data) is np.ndarray:
               self._save_image(data=data, overwrite=overwrite)
           else:
               raise TypeError('Data need to be a numpy.ndarray')



   #------------------------------------------------------------------------------------------------
   def create_from_numpy(self, matrix, col_names, ext_name=None, ext_ver=None, header=None):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
         @param ext_name extension name or number
         @param ext_ver extension version
         @param header list of dictionaries with keys:
                       'card', name', 'value', 'value_orig', 'comment'
      """

      # --- Data
      col_list = []
      for col_name in col_names:
         icol = col_names.index(col_name)
         col_type = self._get_fits_col_type(matrix[:, icol])
         col_data = fits.Column(name=col_name, format=col_type, array=np.ravel(matrix[:, icol]))
         col_list.append(col_data)

      # --- Header
      fits_header = None
      if header is not None:
         fits_header = fits.Header()
         for (k, v) in header.items():
            fits_header[k] = v

      # --- HDUs
      primary_hdu   = fits.PrimaryHDU()
      secondary_hdu = fits.BinTableHDU.from_columns(col_list, header=fits_header)
      if not ext_name is None:
         secondary_hdu.name = ext_name

      # --- Save to disk
      self._cat_data = fits.HDUList(hdus=[primary_hdu, secondary_hdu])
      self._cat_data.writeto(self.fullpath, overwrite=True)

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def close(self):
      if not self._cat_data is None:
          if self.open_mode == FITSCatalog.OpenMode.ReadWrite:
              self.save()
          self._cat_data.close()
          del self._cat_data
          self._cat_data = None
      else:
          raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def save(self):
      if self.open_mode == FITSCatalog.OpenMode.ReadWrite:
          self._cat_data.flush()
      else:
          raise BaseCatalog.OpenModeConflict(open_mode=self.open_mode, open_mode_needed=FITSCatalog.OpenMode.ReadWrite)

   # -----------------------------------------------------------------------------------------------
   def flush(self):
      """ Flush the internal catalog buffers """
#       if self.open_mode == FITSCatalog.OpenMode.ReadWrite:
#          self._cat_data.flush()
      pass

   # -----------------------------------------------------------------------------------------------
   def get_nb_rows(self, hdu_no=None):
      """!
         Get the number of rows in the catalog
         @param hdu_no HDU index (default is 1)
         @return number of rows or 0 if None
      """
      nb_rows = 0
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

      if len(self._cat_data) > 0:
         if hdu_no is None:
            hdu_no = self.hdu_no

         if self._cat_data[hdu_no].size > 0:
            nb_rows = len(self._cat_data[hdu_no].data)

      return nb_rows

   # -----------------------------------------------------------------------------------------------
   def get_nb_cols(self, hdu_no=None):
      """!
         Get the number of columns in the catalog
         @param hdu_no HDU index (default is 1)
         @return number of columns or 0 if None
      """
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)
      if hdu_no is None:
         hdu_no = self.hdu_no
      return len(self.get_col_names(hdu_no))

   # -----------------------------------------------------------------------------------------------
   def get_col_names(self, hdu_no=None):
      """!
         Get the list of column names in the catalog
         @param hdu_no HDU index (default is 1)
         @return list of column names
      """
      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)
      else:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return self._cat_data[hdu_no].columns.names

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def get_info(self):
      """!
         Retrieve some information about the catalog
         @return a dictionary with detailed information
         @note See the fitsio documentation of the info() function for the details
      """
      if not self._cat_data is None:
         return self._cat_data.info()
      else:
         return fits.info(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def col_exists(self, col_name, hdu_no=None):
      """!
         Tell whether a named column exists
         @param hdu_no HDU index (default is 1)
         @return True if column exists, False otherwise
      """
      return col_name in self.get_col_names()

   # -----------------------------------------------------------------------------------------------
   def get_col_index(self, col_name, hdu_no=None):
      """!
         Get the column index from a column name
         @param col_name column name
         @param hdu_no HDU index (default is 1)
         @return column index (zero-based)
         @note Checks disabled for performance reasons
      """
      col_names = self.get_col_names()
      if not col_name in col_names:
         raise BaseCatalog.ColumnNotFound(col_name)
      return col_names.index(col_name)

   # -----------------------------------------------------------------------------------------------
   def get_col_data(self, col_index, hdu_no=None):
      """!
         Return the data of a column from its index
         @param col_name column name
         @param hdu_no HDU index (default is 1)
         @return data associated with the column as a numpy array
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no

         col_data = self._cat_data[hdu_no].columns[col_index].array
         if isinstance(col_data, fits.column.Delayed):
            _ = self._cat_data[hdu_no].data  # trick to force data to be loaded into memory
         return self._cat_data[hdu_no].columns[col_index].array
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_named_col_data(self, col_name, hdu_no=None):
      """!
         Return the data of a column from its index (zero-based)
         @param col_name column name
         @param hdu_no HDU index (default is 1)
         @return data associated with the column as a numpy array
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no

         col_data = self._cat_data[hdu_no].columns[col_name].array
         if isinstance(col_data, fits.column.Delayed):
            _ = self._cat_data[hdu_no].data  # trick to force data to be loaded into memory
         return self._cat_data[hdu_no].columns[col_name].array
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def get_data(self, hdu_no=None):
      """!
         Return data of the specified hdu
         @param hdu_no HDU index
      """
      if not self._cat_data is None:
          if hdu_no is None:
              hdu_no=self.hdu_no

          return self._cat_data[hdu_no].data
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_header(self, hdu_no=None):
      """!
         Return the catalog header as a list of dictionaries
         @param hdu_no HDU index (default is 1)
         @return header
         @see astropy documentation
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return dict(self._cat_data[hdu_no].header.items())
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_headers(self):
      """!
         Return the catalog header as a list of dictionaries
         @return header
         @see fitsio documentation
      """
      headers = []
      try:
         for hdu in self._cat_data:
            headers.append(dict(hdu.header.items()))
      except Exception:
         pass

      return headers

   # -----------------------------------------------------------------------------------------------
   def get_comments(self, hdu_no=None):
      """!
         Return the list catalog comments
         @param hdu_no HDU index (default is 1)
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return [self._cat_data[hdu_no].header.comments[c] \
                                                  for c in self._cat_data[hdu_no].header.keys()]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


   # -----------------------------------------------------------------------------------------------
   def get_col_comments(self, hdu_no=None):
      """!
         Return the list of column comments
         @param hdu_no HDU index (default is 1)
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         hdr_col_types = [tt for tt in self._cat_data[hdu_no].header.keys() if "TTYPE" in tt]
         return [self._cat_data[hdu_no].header.comments[c] for c in hdr_col_types]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_formats(self, hdu_no=None):
      """!
         Get the list of python column formats in the order of columns
         @param hdu_no HDU index (default is 1)
         @return list of column formats in order
      """
      if not self._cat_data is None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return self._cat_data[hdu_no].columns.formats
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


   # -----------------------------------------------------------------------------------------------
   def add_col(self, col_name, col_format=None, col_comment=None, col_data=None, hdu_no=None):
      """!
         Add a Column to the catalog
         @param hdu_no HDU index (default is 1)
         @param col_name column name
         @param col_format column format: (ignored in this class, automatically determined)
         @param col_comment column comment
         @param col_data column data as a numpy array
      """

      # -- Add new column
      col_type = self._get_fits_col_type(col_data)
      new_col = FITSCatalog.Column(name=col_name, format=col_type, comment=col_comment,
                                                  data=col_data)
      if hdu_no is None:
         hdu_no = self.hdu_no
      self._append_col(new_col, hdu_no)

      # --- Update internal data
      self._cat_data.close()
      del self._cat_data
      self._cat_data = fits.open(self.fullpath, mode=self.open_mode, memmap=self.use_memmap)


   # -----------------------------------------------------------------------------------------------
   def remove_col(self, col_index):
      """!
         Delete a column from its index
         @param col_index index of the column to delete
      """

      # TODO: data corruption after calling del_col(). Pyfits bug?
      raise BaseCatalog.FeatureNotImplemented("remove_col()")

   # -----------------------------------------------------------------------------------------------
   def remove_named_col(self, col_name):
      """!
         Delete a named column
         @param col_name name of the column to delete
      """
      if not self._cat_data is None:
         col_index = self.get_col_index(col_name)
         self.remove_col(col_index)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


#
   # ~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _append_col(self, column, hdu_no=None):
      """!
         Append a Column object
         @param column an object derived from BaseCatalog.Column
      """

      if not self._cat_data is None:
         new_col = fits.Column(name=column.name, format=column.format, array=column.data)

         if hdu_no is None:
            hdu_no = self.hdu_no

         orig_table = fits.open(self.fullpath)[hdu_no].data
         orig_cols = orig_table.columns

         new_col = fits.ColDefs([ fits.Column(name=column.name, format=column.format,
                                             array=np.zeros(len(orig_table))) ])
         col_list = orig_cols + new_col
         hdu = fits.BinTableHDU.from_columns(col_list)
         hdu.data[column.name] = column.data
         hdu.writeto(self.fullpath, overwrite=True)

      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def _get_fits_col_type(self, col_data):

      if col_data is None or len(col_data) == 0:
         col_type = 'D'
      elif type(col_data[0]) is int:
         col_type = 'K'
      elif type(col_data[0]) is float:
         col_type = 'D'
      elif type(col_data[0]) is bool:
         col_type = 'L'
      elif type(col_data[0]) is string:
         col_type = 'A'
      else:
         col_type = 'D'

      return col_type

   # -----------------------------------------------------------------------------------------------
   def _get_python_col_type(self, col_type):

      if col_type in ["B", "I", "J", "K"]:
         pcol_type = "%d"
      elif col_type in ["D", "E"]:
         pcol_type = "%f"
      elif col_type in ["A", "C", "M"]:
         pcol_type = "%s"
      elif col_type == "L":
         pcol_type = "%s"
      else:
         pcol_type = "%f"

      return pcol_type

   #------------------------------------------------------------------------------------------------
   # ADDED
   def _save_from_recarray(self, data=None, ext_name=None, sex_cat_path=None):
       """!
            Save a numpy.recarray or astropy.io.fits.fitsrec.FITS_rec into a fits.
            @param data data to store
            @param names names of the diferent column (not needed for dict and rec format)
            @param ext_name name of the HDU where data are stored (DEFAULT = NEW)
            @param sex_cat_path path of the already existing SExtractor catalog to mimic
       """

       if data is None:
           raise ValueError('Data not provided')

       if self.helper.file_exists(self.fullpath):
           if self._cat_data is None:
               self.open()
           if ext_name is None:
               ext_name='new'
           self._cat_data.append(fits.BinTableHDU(data,name=ext_name))
           self.close()
       else:
           if self._SEx_catalog:
              self.create(s_hdu=False, sex_cat_path=sex_cat_path)
              self.open()
              if ext_name is None:
                  ext_name='LDAC_OBJECTS'
              self._cat_data.append(fits.BinTableHDU(data,name=ext_name))
              self.close()
           else:
              self.create(s_hdu=False)
              self.open()
              if ext_name is None:
                  ext_name='new'
              self._cat_data.append(fits.BinTableHDU(data,name=ext_name))
              self.close()

   #------------------------------------------------------------------------------------------------
   # ADDED
   def _save_image(self, data=None, overwrite=False):
       """!
            Save an image into a fits.
            No PrimaryHDU
            @param data image to store
            @param overwrite only used when creating an image fits
       """

       if (data is not None):
           fits.PrimaryHDU(data).writeto(self.fullpath,overwrite=overwrite)
       else:
           raise ValueError('Data or names not provided')


   # -----------------------------------------------------------------------------------------------
   class Column(BaseCatalog.Column):
      """! Represents a column in the catalog """

      def __init__(self, name, format=None, comment=None, data=None):

         """!
            Create a Column object
            @param name name of the column
            @param format Python format of the column (currently not supported)
            @param comment comment string (currently not supported)
            @param data associated column data
#             @param unit column unit, corresponding to TUNIT keyword
#             @param null column null value, corresponding to TNULL keyword
#             @param bscale column bscale value, corresponding to TSCAL keyword
#             @param bzero column bzero value, corresponding to TZERO keyword
#             @param disp column display format, corresponding to TDISP keyword
#             @param start column starting position, corresponding to TBCOL keyword (ASCII table)
#             @param ascii if True, describes a column for an ASCII table
#             @param dim column dimension corresponding to TDIM keyword

         """
         BaseCatalog.Column.__init__(self)

         self._name = name

         if format is None:
            format = "D"
         self._format = format

         if comment is None:
            comment = name
         self._comment = comment

         if data is None:
            data = np.asarray([])
         else:
            if isinstance(data, list):
               self._data = np.asarray(data)
            else:
               self._data = data

      def __str__(self):
         info = "{0}".format(self._cat_col)
         return info

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
            @note overridden from scatalog.BaseColumn im_head, im_object]
         """
         return self._format

      @property
      def comment(self):
         """!
            Get the comment string associated with the column
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


#       @property
#       def unit(self):
#          """!
#             get the column unit, corresponding to TUNIT keyword
#             @return input unit of the column
#          """
#          return self._cat_col.unit
#
#       @property
#       def dim(self):
#          """!
#             get the column dimension corresponding to TDIM keyword
#             @return input dimension of the column
#          """
#          return self._cat_col.dim
#
#       @property
#       def null(self):
#          """!
#             get the column null value, corresponding to TNULL keyword
#             @return input null value of the column
#          """
#          return self._cat_col.null
#
#       @property
#       def bscale(self):
#          """!
#             get the column bscale value, corresponding to TSCAL keyword
#             @return input bscale value of the column
#          """
#          return self._cat_col.bscale
#
#       @property
#       def bzero(self):
#          """!
#             get the column bzero value, corresponding to TZERO keyword
#             @return input bzero value of the column
#          """
#          return self._cat_col.bzero
#
#       @property
#       def disp(self):
#          """!
#             get the column display format, corresponding to TDISP keyword
#             @return display format of the column
#          """
#          return self._cat_col.disp
#
#       @property
#       def start(self):
#          """!
#             get the column starting position (ASCII table only), corresponding to TBCOL keyword
#             @return column starting position
#          """
#          return self._cat_col.start
#
#       @property
#       def ascii(self):
#          """!
#             If True, describes a column for an ASCII table
#             @return True from an ASCII table, False otherwise
#          """
#          return self._cat_col.ascii

      # --- Setters

      @name.setter
      def name(self, name):
         """!
            Set the name of the column
            @param name string
         """
         self._name = name

      @format.setter
      def format(self, format):
         """!
            Set the format of the column
            @param format string
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
            Set the data associated with the column
            @param data column data (list or numpy array)
         """
         self._data = data

#       @unit.setter
#       def unit(self, unit):
#          """!
#             Set the column unit, corresponding to TUNIT keyword
#             @param unit of the column
#          """
#          self._cat_col.unit = unit
#
#       @dim.setter
#       def dim(self, dim):
#          """!
#             Set the column dimension corresponding to TDIM keyword
#             @param dim dimension of the column
#          """
#          self._cat_col.dim = dim
#
#       @null.setter
#       def null(self, null):
#          """!
#             Set the column null value, corresponding to TNULL keyword
#             @param null null value of the column
#          """
#          self._cat_col.null = null
#
#       @bscale.setter
#       def bscale(self, bscale):
#          """!
#             Set the column bscale value, corresponding to TSCAL keyword
#             @param bscale bscale value of the column
#          """
#          self._cat_col.bscale = bscale
#
#       @bzero.setter
#       def bzero(self, bzero):
#          """!
#             Set the column bzero value, corresponding to TZERO keyword
#             @param bzero bzero value of the column
#          """
#          self._cat_col.bzero = bzero
#
#       @disp.setter
#       def disp(self, disp):
#          """!
#             Set the column display format, corresponding to TDISP keyword
#             @param disp display format of the column
#          """
#          self._cat_col.disp = disp
#
#       @start.setter
#       def start(self, start):
#          """!
#             Set the column starting position (ASCII table only), corresponding to TBCOL keyword
#             @param start column starting position
#          """
#          self._cat_col.start = start
#
#       @ascii.setter
#       def ascii(self, ascii):
#          """!
#             If True, describes a column for an ASCII table
#             @param ascii True from an ASCII table
#          """
#          self._cat_col.ascii = ascii


# --------------------------------------------------------------------------------------------------
class LDACFITSCatalog(FITSCatalog):
   """!
      Catalogs management in .FITS LDAC format
      @note The implementation of this class is still work in progress.
      @note At present, .FITS table can only be read, not updated. Also features such like managing
            headers or sorting are not yet implemented.
   """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, fullpath, frame_no=1, open_mode=BaseCatalog.OpenMode.ReadOnly, memmap=False):
      FITSCatalog.__init__(self, fullpath, open_mode, memmap)

      self._format = BaseCatalog.InputFormat.FITS_LDAC  # default input/output format
      self._frame_no = frame_no                         # Frame number of the underlying table
      self._open_mode = open_mode                       # opening mode (see FITSCatalog.OpenMode)
      self._hdu_no = 1                                  # HDU number of the underlying .FITS table


   # -----------------------------------------------------------------------------------------------
   def __str__(self):

      if self._cat_data is not None:
         info = "{0}".format(self.get_info())
      else:
         info = "No information"
      return info


   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~

   # --- Getters

   @property
   def frame_no(self):
      """!
         Get the Frame number of the table
         @return Frame number of the table
      """
      return self._frame_no


   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   #------------------------------------------------------------------------------------------------
   def open(self):

      """!
         Open an existing catalog in FITS LDAC Format
      """
      if self.helper.file_exists(self.fullpath):

         # --- Open catalog file
         hdu_list = fits.open(self.fullpath, mode=self.open_mode, memmap=self.use_memmap)
         [_, _] = awu.ldac.convert_hdu_to_ldac(hdu_list[1])
         #_ = fits.HDUList([fits.PrimaryHDU(), im_head, im_object])

         # --- Update internal data
         self._cat_data   = hdu_list[1:]
         self._cat_imhead = hdu_list[0]

      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   #------------------------------------------------------------------------------------------------
   def create(self, ext_name=None, header=None):
      """!
         Create a new catalog with a FITS LDAC format
         @param ext_name extension name or number
         @param ext_ver extension version
      """
      # --- Create catalog
      col = fits.Column(name='__dummy__', format='D')
      coldefs = fits.ColDefs([col])
      secondary_hdu = fits.BinTableHDU.from_columns(coldefs)
      [im_head, im_object] = awu.ldac.convert_hdu_to_ldac(secondary_hdu)
      hdu_list = fits.HDUList([fits.PrimaryHDU(), im_head, im_object])

      # --- Save to disk
      hdu_list.writeto(self.fullpath, overwrite=True)
      #col = fits.ColDefs([col])

      # --- Update internal data
      self._cat_imhead = hdu_list[0]
      self._cat_data   = hdu_list[1:]


   # -----------------------------------------------------------------------------------------------
   def add_col(self, col_name, col_format=None, col_comment=None, col_data=None, hdu_no=None):
      """!
         Add a Column to the catalog
         @param hdu_no HDU index (default is 1)
         @param col_name column name
         @param col_format column format: (ignored in this class, automatically determined)
         @param col_comment column comment
         @param col_data column data as a numpy array
      """

      # -- Add new column
      if not  self.col_exists(col_name):

         col_type = self._get_fits_col_type(col_data)
         new_col = FITSCatalog.Column(name=col_name, format=col_type, comment=col_comment,
                                                     data=col_data)
         if hdu_no is None:
            hdu_no = self.hdu_no
         self._append_col(new_col, hdu_no)

         # --- Update internal data
         self._cat_data.close()
         del self._cat_data

         ldac_table = awu.ldac.get_table_from_ldac(self.fullpath, frame=self.frame_no)
         hdu_list = awu.ldac.convert_table_to_ldac(ldac_table)

         self._cat_data   = hdu_list[1:]
         self._cat_imhead = hdu_list[0]

   # -----------------------------------------------------------------------------------------------
   def _append_col(self, column, hdu_no=None):
      """!
         Append a Column object
         @param column an object derived from BaseCatalog.Column
      """

      if not self._cat_data is None:
         new_col = fits.Column(name=column.name, format=column.format, array=column.data)

         if hdu_no is None:
            hdu_no = self.hdu_no
         orig_table = fits.open(self.fullpath)[hdu_no].data
         orig_cols = orig_table.columns

         new_col = fits.ColDefs([ fits.Column(name=column.name, format=column.format,
                                              array=np.zeros(len(orig_table))) ])

         col_list = orig_cols + new_col
         hdu = fits.BinTableHDU.from_columns(col_list)
         hdu.data[column.name] = column.data
         [im_head, im_object] = awu.ldac.convert_hdu_to_ldac(hdu)

         new_hdulist = fits.HDUList([fits.PrimaryHDU(), im_head, im_object])
         new_hdulist.writeto(self.fullpath, overwrite=True)

      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)


   #------------------------------------------------------------------------------------------------
   def create_from_numpy(self, matrix, col_names, ext_name=None, ext_ver=None, header=None):
      """!
         Create a new catalog from a two-dimensional numpy array
         @param matrix two-dimensional numpy array
         @param col_names list of column names to use as header
         @param ext_name extension name or number
         @param ext_ver extension version
         @param header list of dictionaries with keys:
                       'card', name', 'value', 'value_orig', 'comment'
      """

      # --- Data
      col_list = []
      for col_name in col_names:
         icol = col_names.index(col_name)

         col_type = self._get_fits_col_type(matrix[:, icol])
         col_data = fits.Column(name=col_name, format=col_type, array=np.ravel(matrix[:, icol]))
         col_list.append(col_data)

      col_defs = fits.ColDefs(col_list)

      # --- Header
      fits_header = None
      if header is not None:
         fits_header = fits.Header()
         for (k, v) in header.items():
            fits_header[k] = v

      # --- HDUs
      primary_hdu   = fits.PrimaryHDU()
      secondary_hdu = fits.BinTableHDU.from_columns(col_defs, header=fits_header)
      if not ext_name is None:
         secondary_hdu.name = ext_name
      [im_head, im_object] = awu.ldac.convert_hdu_to_ldac(secondary_hdu)
      hdu_list = fits.HDUList([primary_hdu, im_head, im_object])

      # --- Save to disk
      hdu_list.writeto(self.fullpath, overwrite=True)

      # --- Update internal data
      self._cat_imhead = hdu_list[0]
      self._cat_data = hdu_list[1:]


# -- EOF scatalog.py
