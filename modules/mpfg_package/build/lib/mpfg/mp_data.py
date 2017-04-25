"""!
   @package mpfg.mp_data data access
   @author Marc Gentile
   @file mp_data.py
   Dataset management
""" 

# -- Python imports
import os, sys
import glob

# --- Module-specific imports
from mp_helper import *             # utility functions

# -------------------------------------------------------------------------------------------------
class Dataset(object):
   
   """! 
       Represent a dataset with images, catalogs and any other kinds of files.
       This class can be subclassed to match any particular dataset directory structure and 
       specificities
   """

   def __init__(self, master, name, base_dir="."):
      """! Dataset constructor """

      self._name = name             # name of dataset
      self._base_dir = base_dir     # base directory where the data reside
      self._helper = Helper()       # helper utility functions
      

   def __str__(self):
      """! String representation of a Dataset object """

      return "<dataset: {0}, {1}>".format(self._name, self._base_dir)

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def type(self):
      """! @return the type of the dataset """
      return  self.__class__.__name__

   @property
   def name(self):
      """! @return the name of the dataset """
      return self._name

   @property
   def base_dir(self):
      """! @return the base directory where the data reside """
      return self._base_dir

   # --- Setters

   @base_dir.setter 
   def base_dir(self, base_dir):
      """! Specify the base directory where the data reside """
      self._base_dir = base_dir 


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def query(self, master, search_expr_list, recurse=True, sort=True):
      """! 
         Query the dataset for files matching a search expression list like [*.fits, ".txt] ) 
         @param master master object instance
         @param search_expr_list regular expression for the name the files to fetch (Unix style)
         @param sort [optional] tell whether to sort the output file names (default True)
         @param recurse [optional] tell whether to walk down directories (default True)
         
         @return list of file paths matching the matching the search expression       
      """

      return self.locate_files(master, search_expr_list, self.base_dir, recurse=recurse, sort=sort)

   # -----------------------------------------------------------------------------------------------
   def locate_files(self, master, pattern_list, directory, sort=True, recurse=True, err_check=True):
      """! 
         Locate files matching a search pattern list

         @param master master object instance
         @param pattern_list Unix-style file search pattern list (e.g. [*.fits, ".txt] )
         @param directory base directory from where to search for matching files
         @param sort [optional]  tell whether to sort the output file paths (default True)
         @param recurse [optional] tell whether to walk down directories (default True)
         @param err_check [optional] tell whether to the validity of check directories
         
         @return list of absolute paths of the files matching the search criteria 
         @note the search is through the entire directory tree, not only the top nodes.
      """

      if err_check and not os.path.isdir(directory):
         self.helper.print_warning("{0} could not be found or is not a directory".format(directory))
         return [] 

      filepaths = []

      if recurse:
         
         for pattern in pattern_list:
            # --- Recursively search the whole directory tree, from top nodes to leaves
            for filepath in self._walk_directory(master, pattern, directory):
               if not os.path.isdir(filepath):
                  filepaths.append(filepath)
      else:

         # --- Search only the top nodes of the directory
         for pattern in pattern_list:
            filepaths.extend([os.path.join(directory, f) for f in os.listdir(directory) \
                                if not os.path.isdir(os.path.join(directory, f)) and \
                                   self.match_file(master, directory, f, pattern)]) 

      if sort:
         filepaths = sorted(filepaths)

      return list(set(filepaths))

   # -----------------------------------------------------------------------------------------------
   def match_file(self, master, directory, filename, pattern):
      """! 
          File matching predicate method. Must return True in case of matching, False otherwise.   
          May be overriden by subclasses to set additional criteria. 
          @param master master object instance
          @param directory directory of filename
          @param filename file name
          @param pattern Unix-like file pattern 
          @return True of a match is found, False otherwise
      """

      return glob.fnmatch.fnmatch(filename, pattern)


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _walk_directory(self, master, pattern, directory):
      """! 
         Recursively locate files matching pattern in a given directory
         @param master master object instance
         @param pattern Unix-style file search pattern (e.g. *.fits)
         @param directory: base directory from where to search for matching files
         
         @return list of files matching the search criteria 
      """

      for path, dirs, files in os.walk(directory):
         for filename in [os.path.abspath(os.path.join(path, filename)) \
            for filename in files if self.match_file(master, path, filename, pattern)]:
               yield filename  




# -- EOF mp_data.py
