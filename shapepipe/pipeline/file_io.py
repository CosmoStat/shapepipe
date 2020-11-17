# -*- coding: utf-8 -*-

"""FILE IO

This file contains methods for file I/O handling.

:Author: Marc Gentile and Axel Guinot

:Version: 1.0.1

"""

# -- Python imports
import os
import string
import re
import numpy as np

import operator
import itertools

from astropy.io import fits
from astropy.table import Table

from modopt.math.stats import sigma_mad


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

   def _file_exists(self, filepath):
      """ File Exists

      Check if input file path is a valid file.

      """
      if not os.path.isfile(filepath):
         return False
      else:
         return True

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
         return "SCatalog *** ERROR ***: Feature: {0} is not implemented in this class".format(self._msg)

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

   class DataNotFound(Exception):
      """ No data found (at given hdu) """

      def __init__(self, filepath, hdu):
         """!
            Exception constructor
            @param filepath file path of the catalog file
         """
         self._filepath = filepath
         self._hdu = hdu

      def __str__(self):
         return 'SCatalog *** ERROR ***: File \'{0}\', hdu={1}: '\
                'data not found'\
                ''.format(self._filepath, self._hdu)

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
         return "SCatalog *** ERROR ***: Catalog: {0} Open Mode {1} not supported".format(self._filepath, self._open_mode)

   # ------------------------------------------------------------------------------------------------
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

      self._open_mode = open_mode      # opening mode (see FITSCatalog.OpenMode)
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

   # ------------------------------------------------------------------------------------------------
   def open(self):

      """! Open an existing catalog """
      if self._file_exists(self.fullpath):

         # --- Open catalog file
         self._cat_data = fits.open(self.fullpath, mode=self.open_mode, memmap=self.use_memmap,
                                    ignore_missing_end=True)
      else:
         raise BaseCatalog.CatalogFileNotFound(self.fullpath)

   # ------------------------------------------------------------------------------------------------
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
              if self._file_exists(sex_cat_path):
                  sex_cat = FITSCatalog(sex_cat_path, hdu_no=1)
                  sex_cat.open()
                  secondary_hdu = sex_cat._cat_data[1]
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

   # ------------------------------------------------------------------------------------------------
   # ADDED
   def copy_hdu(self, fits_file=None, hdu_no=None, hdu_name=None):
       """!
           Copy an HDU from fits
           @param fits_file scatalog where the hdu is
           @param hdu_no hdu to copy default=value use when fits_file was open

           NOTE : the scatalog to copy need to be open
       """

       if fits_file is None:
           fits_file = self

       if self._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(self.fullpath)
       if fits_file._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(fits_file.fullpath)

       if hdu_no is None:
           hdu_no = fits_file.hdu_no

       if hdu_name is None:
           hdu_name = fits_file._cat_data[hdu_no].name

       self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no=hdu_no), name=hdu_name))

   # ------------------------------------------------------------------------------------------------
   #    ADDED
   def apply_mask(self, fits_file=None, hdu_no=None, mask=None, hdu_name=None):
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
           fits_file = self
       if self._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(self.fullpath)
       if fits_file._cat_data is None:
          raise BaseCatalog.CatalogNotOpen(fits_file.fullpath)
       if mask is None:
           raise ValueError('Mask not provided')
       if type(mask) is not np.ndarray:
           raise TypeError('Mask need to be a numpy.ndarray')
       if hdu_no is None:
           hdu_no = fits_file.hdu_no

       if hdu_name is None:
           hdu_name = fits_file._cat_data[hdu_no].name

       if mask.dtype == bool:
           mask = np.where(mask is True)
           self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no)[:][mask], name=hdu_name))
       elif mask.dtype == int:
           self._cat_data.append(fits.BinTableHDU(fits_file.get_data(hdu_no)[:][mask], name=hdu_name))
       else:
           raise TypeError('Mask type has to be int or bool')

   # ------------------------------------------------------------------------------------------------
   # ADDED
   def save_as_fits(self, data=None, names=None, ext_name=None, sex_cat_path=None, image=False, image_header=None, overwrite=False):
       """!
            Save data into an already existing fits or a new one.
            Save data from dict, list, numpy.ndarray, numpy.recarray or astropy.io.fits.fitsrec.FITS_rec (data format in an astropy fits file)
            When creating a new fits to store BinTable data it create a PrimaryHDU.
            When creating a new fits to store Image there is no PrimaryHDU.
            You can create a SExtractor format fits by specifying a SExtractor catalog from where data come from.
            @param data data to store
            @param names names of the diferent column as a list (not needed for dict and rec format)
            @param ext_name name of the HDU where data are stored (DEFAULT = NEW)
            @param sex_cat_path path of the already existing SExtractor catalog to mimic
            @param image if True create a fits image
            @param image_header header to use when saving an image, astropy.io.fits.heade format (optional)
            @param overwrite only used when creating an image fits

            NOTE : to create a SExtractor like fits you need to specify SEx_catalog=True when declaring the FITSCatalog object.
       """

       if self.open_mode != FITSCatalog.OpenMode.ReadWrite:
           raise BaseCatalog.OpenModeConflict(open_mode=self.open_mode, open_mode_needed=FITSCatalog.OpenMode.ReadWrite)

       if data is None:
           raise ValueError('Data not provided')

       if not image:
           if type(data) is dict:
               names = list(data.keys())
               it = list(range(len(names)))
               if len(names) == 1:
                  data = np.array(data[names[0]])
               else:
                  data = [np.array(data[i]) for i in names]
               self._save_to_fits(data, names, it, ext_name, sex_cat_path)

           elif type(data) is np.recarray:
               names = list(data.dtype.names)
               it = names
               self._save_to_fits(data, names, it, ext_name, sex_cat_path)

           elif type(data) is fits.fitsrec.FITS_rec:
               self._save_from_recarray(data, ext_name, sex_cat_path)

           elif type(data) is np.ndarray:
               if names is None:
                  if data.dtype.names is not None:
                     names = data.dtype.names
                     it = names
                  else:
                     raise ValueError('Names not provided')
               else:
                  it = range(len(names))
               self._save_to_fits(data, names, it, ext_name, sex_cat_path)

           elif type(data) is list:
               if names is None:
                   raise ValueError('Names not provided')
               it = range(len(names))
               data = np.asarray(data)
               self._save_to_fits(data, names, it, ext_name, sex_cat_path)

           elif type(data) is Table:
               if names is None:
                   raise ValueError('Names not provided')
               it = names
               self._save_to_fits(data, names, it, ext_name, sex_cat_path)

       else:
           if type(data) is np.ndarray:
               self._save_image(data=data, header=image_header, overwrite=overwrite)
           else:
               raise TypeError('Data need to be a numpy.ndarray')

   # ------------------------------------------------------------------------------------------------
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
      primary_hdu = fits.PrimaryHDU()
      secondary_hdu = fits.BinTableHDU.from_columns(col_list, header=fits_header)
      if ext_name is not None:
         secondary_hdu.name = ext_name

      # --- Save to disk
      self._cat_data = fits.HDUList(hdus=[primary_hdu, secondary_hdu])
      self._cat_data.writeto(self.fullpath, overwrite=True)

   # -----------------------------------------------------------------------------------------------
   # MODIFIED
   def close(self):
      if self._cat_data is not None:
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
      if self._cat_data is not None:
         return self._cat_data.info()
      else:
         return fits.info(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   # ADDED
   def get_ext_name(self, hdu_no=None):
      """
         Return the name of an extansion or all of them.
         @param hdu_no index of the hdu to return. If None return all of them
         @return a string of the hdu name or a list of string.
      """

      if hdu_no is None:
         n = [self._cat_data[i].name for i in range(len(self._cat_data))]
      else:
         n = self._cat_data[int(hdu_no)].name

      return n

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
      if col_name not in col_names:
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
      if self._cat_data is not None:
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
      if self._cat_data is not None:
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
      if self._cat_data is not None:
          if hdu_no is None:
              hdu_no = self.hdu_no

          dat = self._cat_data[hdu_no].data
          if dat is None:
              raise BaseCatalog.DataNotFound(self.fullpath, self.hdu_no)
          return dat
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
      if self._cat_data is not None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return dict(self._cat_data[hdu_no].header.items())
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # ------------------------------------------------------------------------------------------------
   def get_header_value(self, request, hdu_no=None):
       """!
          Return the value of a parameters or a linear combination of parameters and/or numbers
          @param request parameter or a linear combination of parameters in string format
          @param hdu_no hdu containing the header
          @return result of the request as float
       """

       if request is None:
           raise ValueError('request not provided')
       if type(request) is not str:
           raise TypeError('request has to be a string')

       if hdu_no is None:
           hdu_no = self._hdu_no

       header = self.get_header(hdu_no=hdu_no)
       if header is None:
           raise ValueError('Empty header in the hdu : {0}'.format(hdu_no))

       return interpreter(string=request, catalog=header).result

   # -----------------------------------------------------------------------------------------------
   def add_header_card(self, key, value=None, comment=None, hdu_no=None):
       """!
          Add a card in the header of the specified hdu
          @param key the key to add
          @param value the value of the key
          @param comment comment for the key
          @param hdu_no hdu where the header to modified is
       """

       if self.open_mode != FITSCatalog.OpenMode.ReadWrite:
           raise BaseCatalog.OpenModeConflict(open_mode=self.open_mode, open_mode_needed=FITSCatalog.OpenMode.ReadWrite)

       if self._cat_data is None:
           raise BaseCatalog.CatalogNotOpen(self.fullpath)

       if hdu_no is None:
           hdu_no = self._hdu_no

       card = []
       if key is None:
           raise ValueError('key not provided')
       else:
           card.append(key)

       if value is not None:
           card.append(value)
       else:
           if comment is not None:
               card.append('')

       if comment is not None:
           card.append(comment)

       card = tuple(card)

       self._cat_data[hdu_no].header.append(card, end=True)

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
      if self._cat_data is not None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return [self._cat_data[hdu_no].header.comments[c]
                 for c in self._cat_data[hdu_no].header.keys()]
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   def get_col_comments(self, hdu_no=None):
      """!
         Return the list of column comments
         @param hdu_no HDU index (default is 1)
      """
      if self._cat_data is not None:
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
      if self._cat_data is not None:
         if hdu_no is None:
            hdu_no = self.hdu_no
         return self._cat_data[hdu_no].columns.formats
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # -----------------------------------------------------------------------------------------------
   # def add_col(self, col_name, col_format=None, col_comment=None, col_data=None, hdu_no=None):
   #    """!
   #       Add a Column to the catalog
   #       @param hdu_no HDU index (default is 1)
   #       @param col_name column name
   #       @param col_format column format: (ignored in this class, automatically determined)
   #       @param col_comment column comment
   #       @param col_data column data as a numpy array
   #    """

   #    # -- Add new column
   #    col_type = self._get_fits_col_type(col_data)
   #    new_col = FITSCatalog.Column(name=col_name, format=col_type,
   #                                 comment=col_comment, data=col_data)
   #    if hdu_no is None:
   #       hdu_no = self.hdu_no
   #    self._append_col(new_col, hdu_no)

   #    # --- Update internal data
   #    self._cat_data.close()
   #    del self._cat_data
   #    self._cat_data = fits.open(self.fullpath, mode=self.open_mode, memmap=self.use_memmap)

   # -----------------------------------------------------------------------------------------------
   # ADDED
   def add_col(self, col_name, col_data, hdu_no=None, ext_name=None, new_cat=False, new_cat_inst=None):
      """
         Add a Column to the catalog
         @param col_name column name
         @param col_data column data as a numpy array
         @param hdu_no HDU index where to add the column
         @param ext_name change the name of the extansion (optional)
         @param new_cat if True will save the change into a new catalog
         @param new_cat_inst io.FITSCatalog object for the new catalog
      """
      if new_cat:
          open_mode = new_cat_inst.open_mode
          output_path = new_cat_inst.fullpath
      else:
          open_mode = self.open_mode
          output_path = self.fullpath

      if self._cat_data is None:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

      if open_mode != FITSCatalog.OpenMode.ReadWrite:
         raise BaseCatalog.OpenModeConflict(open_mode=open_mode, open_mode_needed=FITSCatalog.OpenMode.ReadWrite)

      if type(col_data) != np.ndarray:
         TypeError("col_data must be a numpy.ndarray")

      if hdu_no is None:
         hdu_no = self.hdu_no
      if ext_name is None:
         ext_name = self._cat_data[hdu_no].name

      n_of_hdu = len(self._cat_data)
      old_hdu_prev = []
      for i in range(0, hdu_no):
         old_hdu_prev.append(self._cat_data[i])
      old_hdu_next = []
      for i in range(hdu_no+1, n_of_hdu):
         old_hdu_next.append(self._cat_data[i])

      new_fits = fits.HDUList(old_hdu_prev)

      col_list = self._cat_data[hdu_no].data.columns

      data_type = self._get_fits_col_type(col_data)
      data_shape = col_data.shape[1:]
      dim = str(tuple(data_shape))
      mem_size = 1
      if len(data_shape) != 0:
         for k in data_shape:
            mem_size *= k
         data_format = '{0}{1}'.format(mem_size, data_type)
         new_col = fits.ColDefs([fits.Column(name=col_name, format=data_format,
                                 array=col_data, dim=dim)])
         col_list += new_col
      elif data_type == 'A':
         mem_size *= len(max(col_data, key=len))
         data_format = '{0}{1}'.format(mem_size, data_type)
         new_col = fits.ColDefs([fits.Column(name=col_name, format=data_format,
                                 array=col_data, dim=str((mem_size,)))])
         col_list += new_col
      else:
         data_format = '{0}{1}'.format(mem_size, data_type)
         new_col = fits.ColDefs([fits.Column(name=col_name, format=data_format,
                                 array=col_data)])
         col_list += new_col

      new_fits.append(fits.BinTableHDU.from_columns(col_list, name=ext_name))

      new_fits += fits.HDUList(old_hdu_next)

      new_fits.writeto(output_path, overwrite=True)

      if not new_cat:
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
      if self._cat_data is not None:
         col_index = self.get_col_index(col_name)
         self.remove_col(col_index)
      else:
         raise BaseCatalog.CatalogNotOpen(self.fullpath)

   # ~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _append_col(self, column, hdu_no=None):
      """!
         Append a Column object
         @param column an object derived from BaseCatalog.Column
      """

      if self._cat_data is not None:
         new_col = fits.Column(name=column.name, format=column.format, array=column.data)

         if hdu_no is None:
            hdu_no = self.hdu_no

         orig_table = fits.open(self.fullpath)[hdu_no].data
         orig_cols = orig_table.columns

         new_col = fits.ColDefs([fits.Column(name=column.name, format=column.format,
                                             array=np.zeros(len(orig_table)))])
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
      elif type(col_data[0]) in [np.int16]:
         col_type = 'I'
      elif type(col_data[0]) in [np.int32]:
         col_type = 'J'
      elif type(col_data[0]) in [int, np.int64]:
         col_type = 'K'
      elif type(col_data[0]) in [float, np.float16, np.float32, np.float64]:
         col_type = 'D'
      elif type(col_data[0]) is bool:
         col_type = 'L'
      elif type(col_data[0]) in [str, np.str, np.str_, np.str0]:
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

   # ------------------------------------------------------------------------------------------------
   # ADDED
   def _save_to_fits(self, data, names, it, ext_name=None, sex_cat_path=None):
      """!
            Save array of data as fits with their associated column names.
            @param data array withe the data
            @param names list with the column names
            @param it iterator
            @param ext_name name of the HDU where data are stored (DEFAULT = NEW)
            @param sex_cat_path path of the already existing SExtractor catalog to mimic
      """

      if data is None:
         raise ValueError('Data not provided')

      if self._file_exists(self.fullpath):
          if self._cat_data is None:
              self.open()
          if ext_name is None:
              ext_name = 'new'
      else:
          if self._SEx_catalog:
             self.create(s_hdu=False, sex_cat_path=sex_cat_path)
             self.open()
             if ext_name is None:
                 ext_name = 'LDAC_OBJECTS'
          else:
              self.create(s_hdu=False)
              self.open()
              if ext_name is None:
                  ext_name = 'new'

      if len(names) == 1:
          data = np.array([data])
      col_list = []
      for i in it:
          data_shape = data[i].shape[1:]
          dim = str(tuple(data_shape))
          name = names[it.index(i)]
          data_type = self._get_fits_col_type(data[i])
          mem_size = 1
          if len(data_shape) != 0:
              for k in data_shape:
                  mem_size *= k
              data_format = '{0}{1}'.format(mem_size, data_type)
              col_list.append(fits.Column(name=name, format=data_format,
                                          array=data[i], dim=dim))
          elif data_type == 'A':
            mem_size *= len(max(data[i], key=len))
            data_format = '{0}{1}'.format(mem_size, data_type)
            col_list.append(fits.Column(name=name, format=data_format,
                                        array=data[i], dim=str((mem_size,))))
          else:
              data_format = '{0}{1}'.format(mem_size, data_type)
              col_list.append(fits.Column(name=name, format=data_format,
                                          array=data[i]))

      self._cat_data.append(fits.BinTableHDU.from_columns(col_list, name=ext_name))
      self.close()

   # ------------------------------------------------------------------------------------------------
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

       if self._file_exists(self.fullpath):
           if self._cat_data is None:
               self.open()
           if ext_name is None:
               ext_name = 'new'
           self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
           self.close()
       else:
           if self._SEx_catalog:
              self.create(s_hdu=False, sex_cat_path=sex_cat_path)
              self.open()
              if ext_name is None:
                  ext_name = 'LDAC_OBJECTS'
              self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
              self.close()
           else:
              self.create(s_hdu=False)
              self.open()
              if ext_name is None:
                  ext_name = 'new'
              self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
              self.close()

   # ------------------------------------------------------------------------------------------------
   # ADDED
   def _save_image(self, data=None, header=None, overwrite=False):
       """!
            Save an image into a fits.
            No PrimaryHDU
            @param data image to store
            @param header external image header (optional)
            @param overwrite only used when creating an image fits
       """

       if (data is not None):
           fits.PrimaryHDU(data, header).writeto(self.fullpath, overwrite=overwrite)
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


# --------------------------------------------------------------------------------------------------

class interpreter(object):
    """Interpreter class

    Class to handle operation/comparison in a string

    Parameters
    ----------
    string : str
        String to interpret
    catralog : dict, recarray or str
        If type(catalog) == str : assume a SExtractor fits catalog and read parameter in it
        else : assume the catalog is already open and look into it
    make_compare : bool
        If true assume a comparison in the string
    mask_dict : dict
        dictionary containing mask usable for the operation
    autorun : bool
        If true return directly the result

    """

    def __init__(self, string, catalog, make_compare=False, mask_dict=None):

        if type(string) is not str:
            raise ValueError('string has to be str type')
        else:
            self._string = string

        if catalog is not None:
            if type(catalog) is str:
                f = FITSCatalog(catalog, SEx_catalog=True)
                f.open()
                self._cat = f.get_data()
                f.close()
            else:
                self._cat = catalog
        else:
            raise ValueError('catalog not provided')

        self._make_compare = make_compare

        if mask_dict is not None:
            self._mask = mask_dict

        self._init_stat_function()
        self._comp_dict = {'<': operator.lt,
                           '>': operator.gt,
                           '<=': operator.le,
                           '>=': operator.ge,
                           '==': operator.eq,
                           '!=': operator.ne}

        self.result = self.interpret(self._string, self._make_compare)

    def interpret(self, string, make_compare=False, make_func=True, make_operate=True):
        """Interpret

        This function handles the different possible operations

        Parameters
        ----------
        str:  string
            string to interpret
        make_compare : bool
            If true look for a comparison
        make_func : bool
            If true look for a function
        make_operate : bool
            If true look for an operation

        Returns
        -------
        array or float
            Result of the current operation.

        Notes
        -----
        This is a recursive function.

        """

        if make_compare:
            result = self._compare(string)
        else:
            if make_operate:
                string_split = re.split(r'\*|\/|\-|\+\s*(?![^()]*\))', string)
                result = self._operate(string, string_split)
            else:
                if make_func:
                    result = self._apply_func(string)
                else:
                    result = self._get_value(string)

        return result

    def _compare(self, string):
        """Handle comparison in a string

        This function transform condition in a string to real condition.

        Parameters
        ----------
        string : str
            strind containing the comparison.

        """

        comp = '<|>|<=|>=|==|!='

        if len(re.split(comp, string)) != 2:
            raise Exception('Only one comparison in [<, >, <=, >=, ==, !=] per lines')

        for i in ['<=', '>=', '<', '>', '==', '!=']:
            terms = re.split(i, string)
            if len(terms) == 2:
                self._make_compare = False
                first = self.interpret(terms[0], self._make_compare)
                second = self.interpret(terms[1], self._make_compare)

                return self._comp_dict[i](first, second)

    def _apply_func(self, string):
        """Parse input string for function name and apply function.

        Parameters
        ----------
        str: string
            input string

        Returns
        -------
        float
            result of the function

        """

        s = re.split(r'\(|\)', string)

        if len(s) == 1:
            return self.interpret(s[0], self._make_compare, make_func=False, make_operate=False)
        elif len(s) == 3:
            try:
                ss = re.split(',', s[1])
                if len(ss) > 1:
                    p = [self.interpret(i, self._make_compare, make_func=False, make_operate=True) for i in ss]
                    return self._stat_func[s[0]](*p)
                else:
                    return self._stat_func[s[0]](self.interpret(s[1], self._make_compare, make_func=False, make_operate=True))
            except:
                raise Exception('Unknown function: {0}'.format(s[0]))
        else:
            raise Exception('Only one function can be applied.'
                            'Problem with the term: {0}'.format(string))

    def _init_stat_function(self):
        """Initialise available stat functions

        Create a dictionary containing the functions.

        """

        self._stat_func = {}

        self._stat_func['mean'] = np.mean
        self._stat_func['median'] = np.median
        self._stat_func['mode'] = self._mode
        self._stat_func['sqrt'] = np.sqrt
        self._stat_func['pow'] = pow
        self._stat_func['log'] = np.log
        self._stat_func['log10'] = np.log10
        self._stat_func['exp'] = np.exp
        self._stat_func['std'] = np.std
        self._stat_func['var'] = np.var
        self._stat_func['sigma_mad'] = self._sigma_mad
        self._stat_func['len'] = len
        self._stat_func['min'] = min
        self._stat_func['max'] = max
        self._stat_func['homogen'] = self._test_homogeneity

    def _mean(self, input):
        """Mean

        Compute the mean of a distribution.

        Parameters
        ----------
        input : numpy.ndarray
            Numpy array containing the data.

        Returns
        -------
        m : float
            mean, if input array has size > 0;
            -1, else
        """

        cat_size = len(input)
        if cat_size == 0:
            return -1
        else:
            return np.mean()

    def _mode(self, input, eps=0.001, iter_max=1000):
        """Mode

        Compute the mode, the most frequent value of a continuous distribution.

        Parameters
        ----------
        input : numpy.ndarray
            Numpy array containing the data.
        eps : float, optional
            Accuracy to achieve (default is 0.001)

        Returns
        -------
        m : float
            mode, if input array has 10 or more elements;
            median, if input array has >0 and <10 elements;
            -1, if input array has 0 elements
        """

        cat_size = len(input)
        if cat_size > 100:
            bins = int(float(cat_size)/10.)
        elif cat_size >= 20:
            bins = int(float(cat_size)/5.)
        elif cat_size > 0:
            return np.median(input)
        else:
            return -1

        data = input
        diff = eps+1.

        k = 0
        while diff > eps:
            hist = np.histogram(data, bins)
            if hist[0].max() == 1:
                break

            b_min = hist[1][hist[0].argmax()]
            b_max = hist[1][hist[0].argmax()+1]

            diff = b_max - b_min

            data = data[(data > b_min) & (data < b_max)]

            if k == iter_max:
                break
            k += 1

        if k == iter_max:
            raise ValueError('Mode computation failed')
        else:
            m = (b_min + b_max) / 2.
            return m

    def _sigma_mad(self, input):
        """Mean absolute deviation

        Compute median absolute deviation (MAD).

        Parameters
        ----------
        input : numpy.nparray
            input data

        Returns
        -------
        mad : float
            MAD, if input size > 0;
            -1 if input size is 0

        """

        if len(input) == 0:
            return -1
        else:
            return sigma_mad(input)

    def _test_homogeneity(self, *args):
        """Test homogeneity

        Test homogeneity on 1D or 2D space.

        Parameters
        ----------
        param1 : numpy.ndarray
            Array on which the homogeneity test is performed
        param2 : numpy.ndarray [optional]
            Array on which the homogeneity test is performed
        nb_cells : int
            Number of cells in the space. (note : has to be a square number)

        Returns
        -------
        float
            Percentage of inhomogeneity compared to worse possible case (based on the standard deviation)

        """

        if len(args) == 2:
            n_param = 1
            param = [args[0]]
            n_cells = args[1]
        elif len(args) == 3:
            n_param = 2
            param = [args[0], args[1]]
            n_cells = args[2]
        else:
            raise ValueError('Inputs should be param_1, param_2 [optional], n_cells')

        if n_param == 2:
            if len(param[0]) != len(param[1]):
                raise ValueError('Both param_1 and param_2 must have the same lenght : {0}, {1}'.format(len(param[0]), len(param[1])))

        if np.sqrt(n_cells) % 1 != 0:
            raise ValueError('N_cells must be a square number')

        n_tot = len(param[0])
        homo_ratio = float(n_tot) / float(n_cells)

        param_min = []
        param_max = []
        for i in param:
            step = (np.max(i) - np.min(i)) / pow(n_cells, 1./float(n_param))
            param_min.append([j for j in np.arange(np.min(i), np.max(i), step)])
            param_max.append([j for j in np.arange(np.min(i) + step, np.max(i) + step, step)])

        if n_param == 1:
            n_obj = np.asarray([float(len(np.where((param[0] >= param_min[0][i]) & (param[0] <= param_max[0][i]))[0])) for i in range(int(n_cells))])
        elif n_param == 2:
            it = itertools.product(range(int(np.sqrt(n_cells))), repeat=2)
            n_obj = np.asarray([float(len(np.where((param[0] >= param_min[0][i]) & (param[0] <= param_max[0][i]) & (param[1] >= param_min[1][j]) & (param[1] <= param_max[1][j]))[0])) for i, j in it])

        actual_std = np.std(n_obj/homo_ratio)

        worse_std = np.zeros((int(n_cells), 1))
        worse_std[0] = n_tot/homo_ratio
        worse_std = np.std(worse_std)

        return actual_std/worse_std*100.

    def _operate(self, string, string_split):
        """Handle operation in a string

        Make operation between catalog's parameters and/or numbers

        Parameters
        ----------
        string : str
            Parameter or linear combination of parameters.

        Returns
        -------
        float
            Result of the operation

        Note
        ----
        It's used as a recursive function

        """

        op = r'\*|\/|\-|\+\s*(?![^()]*\))'
        if string is None:
            raise ValueError("Parameter not specified")
        if string_split is None:
            raise ValueError("Parameters splited not specified")

        if len(re.split(op, string)) == 1:
            return self.interpret(string, make_operate=False)

        tmp = self._string_op_func(re.split(r'\+\s*(?![^()]*\))', string), string_split, operator.add, 0)
        if not np.isscalar(tmp) or tmp != 'pass':
            return tmp
        else:
            tmp = self._string_op_func(re.split(r'\-\s*(?![^()]*\))', string), string_split, operator.sub, 'init')
            if not np.isscalar(tmp) or tmp != 'pass':
                return tmp
            else:
                tmp = self._string_op_func(re.split(r'\*\s*(?![^()]*\))', string), string_split, operator.mul, 1)
                if not np.isscalar(tmp) or tmp != 'pass':
                    return tmp
                else:
                    return self._string_op_func(re.split(r'\/\s*(?![^()]*\))', string), string_split, operator.truediv, 'init')

    def _string_op_func(self, string_op, string_split, op, tmp):
        r"""Perform a specified operation

        This function handle the posible operation between parameters.

        Parameters
        ----------
        string_op : list
            List of parameters to operate.
        string_split : list
            The different parameter splitted using '\*|\/|\-|\+\s*(?![^()]*\))' as delimiter.
        op : func
            The kind of operation provided as an operator function
            (Example : operator.sub).
        tmp : str or float
            Temporary result of the global operation or value to
            initiate operation.

        Returns
        -------
        float or 'pass'
            Result of the operation or 'pass' if their is remaining operations.

        """

        if len(string_op) > 2:
            for i in string_op:
                if tmp == 'init':
                    tmp = self._operate(i, string_split)
                else:
                    tmp = op(tmp, self._operate(i, string_split))
            return tmp
        elif len(string_op) == 2:
            if string_op[0] in string_split:
                first = self.interpret(string_op[0], make_operate=False)
            else:
                first = self._operate(string_op[0], string_split)
            if string_op[1] in string_split:
                second = self.interpret(string_op[1], make_operate=False)
            else:
                second = self._operate(string_op[1], string_split)
            return op(first, second)
        else:
            return 'pass'

    def _get_value(self, string):
        """Get Value

        Return the value of the corresponding parameter. Or return a float with a number as parameter.

        Parameters
        ----------
        string : str
            Parameter of the catalog.

        Returns
        -------
        float
            Vvalue of the parameter. Or float

        Note
        ----
        You can't perform operations here!
        """

        if string is None:
            raise ValueError("Parameter not specified")

        try:
            string_value = float(string)
            return string_value
        except:
            s = re.split(r'\{|\}', string)
            if len(s) == 1:
                try:
                    return self._cat[string]
                except:
                    raise ValueError('string has to be a float or a catalog parameter. {0} not found'.format(string))
            if len(s) == 3:

                if s[1] in self._mask.keys():
                    try:
                        return self._cat[s[0]][self._mask[s[1]]]
                    except:
                        raise ValueError('String has to be a catalog parameter. {0} not found'.format(s[0]))
                else:
                    raise ValueError('mask has to be provided. {0} not found in mask'.format(s[1]))


# -- EOF scatalog.py
