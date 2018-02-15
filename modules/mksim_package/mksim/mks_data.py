"""! 
   @package mpfx.mpfx_data Dataset management
   @author Marc Gentile
   @file mpfx_data.py
   Dataset management
""" 

# -- Python imports
import os
import glob
import numpy

# -- External imports
from mpfx.mpfx_data import MpfxDataset


# -------------------------------------------------------------------------------------------------
class MksDataset(MpfxDataset):
   
   """! 
       Represent a source dataset.
   """

   def __init__(self, master, name, base_dir=".", dir_list=[], dir_recurse=True):

      MpfxDataset.__init__(self, master, name, base_dir)


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def is_catalog(self, master, file_path, 
                                dataset_section_name="IMAGE_DATASET"):
      """! Attempt to tell whether a file is a catalog
           @param master master object instance
           @param file_path absolute file path of a file
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note The behavior of this method may be overriden
                  to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_catalog = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
      else:   
         is_catalog = self.is_star_catalog(master, file_path) or \
                      self.is_galaxy_catalog(master, file_path)
         
      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_image(self, master, file_path, dataset_section_name="IMAGE_DATASET"):
      """! Attempt to tell whether a file is an image
           @param master master object instance
           @param file_path absolute file path of a file
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note The behavior of this method may be overriden
                  to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_image = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            is_image = self.is_star_image(master, file_path) or \
                      self.is_galaxy_image(master, file_path)

      return is_image


   # -----------------------------------------------------------------------------------------------
   def is_star_catalog(self, master, file_path,
                             dataset_section_name="CATALOG_DATASET"):
      """! Attempt to determine if a file is a star catalog 
           @param master master object instance
           @param file_path absolute file path of a file
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note  to change this behavior, override this method in a sub-class of MPfx
      """
   
      _, filename = os.path.split(file_path)

      is_catalog = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         
      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_catalog(self, master, file_path,
                               dataset_section_name="CATALOG_DATASET"):
      """! Attempt to determine if a file is a galaxy catalog 
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note  to change this behavior, override this method in a sub-class of MPfx
      """

      _, filename = os.path.split(file_path)

      is_catalog = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list(
                                                      "FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         
      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_star_image(self, master, file_path,
                           dataset_section_name="IMAGE_DATASET"):
      """! Attempt to determine if a file is a star image 
           @param master master object instance
           @param file_path absolute file path of a file
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note  the behavior of this method should probably be overriden in a subclass of Mpfx
                  in order to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_image = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         
      return is_image

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_image(self, master, file_path,
                             dataset_section_name="IMAGE_DATASET"):
      """! Attempt to determine if a file is a star image 
           @param master master object instance
           @param file_path absolute file path of a file
           @param dataset_section_name config section name with a @c FILENAME_PATTERNS key
           @note  the behavior of this method should probably be overriden in a subclass of Mpfx
                  in order to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_image = False
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])

      return is_image



# -- EOF mks_data.py

