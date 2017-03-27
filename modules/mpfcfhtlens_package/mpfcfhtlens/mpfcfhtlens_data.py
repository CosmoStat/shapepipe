"""! 
   @package mpfcfhtlens.mpfcfhtlens_data Dataset management
   @author Marc Gentile, Martin Kilbinger
   @file mpfcfhtlens_data.py
   Dataset management
""" 

# -- Python imports
import os, sys
import re

# -- External imports
from mpfx.mpfx_data import *


# -------------------------------------------------------------------------------------------------
class MpfcfhtlensDataset(MpfxDataset):
   
   """! 
       Represent a source dataset.
   """

   def __init__(self, master, name, base_dir=".", dir_list=[], dir_recurse=True):

      MpfxDataset.__init__(self, master, name, base_dir, dir_list, dir_recurse)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def query(self, master, file_search_pattern, dir_list=[],   
                   image_list=[], image_range=[-1,-1], 
                   sort=True, recurse=True):
      """! 
         Query files matching a set of directory and file patterns
         @param master master object instance
         @param file_search_pattern Unix-style file search pattern list (like "*.fits")
         @param dir_list [optional] a list of specific directories to search under base directory
         @param image_list a list of specific image numbers to search, like [0, 100]
         @param image_range a range of specific numbers images to search like [50, 100]
         @param sort [optional] tell whether to sort the output file names (default @c True)
         @param recurse [optional] tell whether to walk down directories (default @c True)
               
         @return a list of file paths matching the search criteria      
      """

      return MpfxDataset.query(self, master, file_search_pattern, dir_list, image_list, image_range, 
                                     sort, recurse)

   # -----------------------------------------------------------------------------------------------
   def match_file(self, master, path, filename, pattern):
      """! 
         Reformat a file name so that it complies with the standard MPF filename format: 
         @code <img_name>-<img_no>-<epoch-no>.<ext> @endcode
         @param master master object instance
         @param path (without the filename) or the file whose name we want to format
         @param filename file name to format 
         @param pattern search pattern
         @return formatted file name (not the full path)
      """
      return MpfxDataset.match_file(self, master, path, filename, pattern)

   # -----------------------------------------------------------------------------------------------
   def format_filename(self, master, path, filename):
      """! 
         Reformat the filename
         May be overriden by subclasses to set additional criteria. 
         @param master master object instance
         @param path path of the file (without the file name)
         @param filename name of the file
         @retval True if the file with name @c filename and path @c path matches @c pattern
         @retval False if no match
      """

      # --- Set customized filename here if needed
   
      return filename

   # -----------------------------------------------------------------------------------------------
   def is_star_catalog(self, master, file_path):
      """! Tell whether a file is a star catalog 
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @note by default, it is assumed  that the filename includes both the names 
                "star" and "catalog". The behavior of this method should probably be overriden
                to comply with the naming conventions of a specific source dataset
      """

      return MpfxDataset.is_star_catalog(self, master, file_path)

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_catalog(self, master, file_path):
      """! Tell whether a file is a galaxy catalog 
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @note by default, it is assumed  that the filename includes both the names 
                "star" and "catalog". The behavior of this method should probably be overriden
                to comply with the naming conventions of a specific source dataset
      """

      return MpfxDataset.is_galaxy_catalog(self, master, file_path)

   # -----------------------------------------------------------------------------------------------
   def is_star_image(self, master, file_path):
      """! Tell whether a file is a star image 
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @note by default, it is assumed  that the filename includes both the names 
                "star" and "catalog". The behavior of this method should probably be overriden
                to comply with the naming conventions of a specific source dataset
      """

      return MpfxDataset.is_star_image(self, master, file_path)

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_image(self, master, file_path):
      """! Tell whether a file is a galaxy image 
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @note by default, it is assumed  that the filename includes both the names 
                "star" and "catalog". The behavior of this method should probably be overriden
                to comply with the naming conventions of a specific source dataset
      """

      return MpfxDataset.is_galaxy_image(self, master, file_path)

   # -----------------------------------------------------------------------------------------------
   def is_catalog(self, master, file_path):
      """! Tell whether a file is that of a catalog 
           @param master master object instance
           @param file_path absolute file path of a file
      """

      return MpfxDataset.is_catalog(self, master, file_path)

   # -----------------------------------------------------------------------------------------------
   def is_image(self, master, file_path):
      """! Tell whether a file is that of an image 
           @param master master object instance
           @param file_path absolute file path of a file
      """

      return MpfxDataset.is_image(self, master, file_path)



# -------------------------------------------------------------------------------------------------
class MpfcfhtlensMultiEpochDataset(MpfxMultiEpochDataset):
   
   """! 
       Represent a source dataset that takes care of files with multiple exposures (epochs). 
   """

   def __init__(self, master, name, base_dir=".", dir_list=[], dir_recurse=True):

      MpfxMultiEpochDataset.__init__(self, master, name, base_dir, dir_list, dir_recurse)


   # -----------------------------------------------------------------------------------------------
   def format_filename(self, master, path, filename):
      """! 
         File matching predicate method. Must return True in case of matching, False otherwise.   
         May be overriden by subclasses to set additional criteria. 
         @param master master object instance
         @param path path of the file (without the file name)
         @param filename name of the file
         @retval True if the file with name @c filename and path @c path matches @c pattern
         @retval False if no match
      """

      # --- Set customized filename here if needed

#      filepath = os.path.join(path, filename)
#      print filepath, "Star catalog:", self.is_star_catalog(master, filepath),
#      print "Gal catalog:", self.is_galaxy_catalog(master, filepath),
#      print "Star image:", self.is_star_image(master, filepath),
#      print "Gal image:", self.is_galaxy_image(master, filepath)  

      return filename



# -- EOF mpfcfhtlens_data.py
