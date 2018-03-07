"""!
   @package mpfx.mpfx_data Dataset management
   @author Marc Gentile
   @file mpfx_data.py
   Dataset management
"""

# -- Python imports
import os
import re
from astropy.io.fits import getheader

# -- External imports
from mpfg.mp_data import *


# -------------------------------------------------------------------------------------------------
class MpfxDataset(Dataset):

   """!
       Represent a source dataset.
   """

   def __init__(self, master, name, base_dir=".", dir_list=[], dir_recurse=True):

      Dataset.__init__(self, master, name, base_dir)

      self._dir_list = dir_list        # list of directories to search under base directory
      self._dir_recurse = dir_recurse  # recursive directory search or not
      self._file_matcher = None        # FileMatcher onject to locate files

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def file_matcher(self):
      """! @return the file macthing method """
      return self._file_matcher

   @property
   def dir_list(self):
      """! @return the optional list of directories to search under base directory """
      return self._dir_list

   @property
   def dir_recurse(self):
      """! @return whether the file search is recursive or not """
      return self._dir_recurse

   # --- Setters

   @file_matcher.setter
   def file_matcher(self, file_matcher):
      """!
         Specify the file macthing method
         @param file_matcher the file macthing method
      """
      self._file_matcher = file_matcher

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

      # --- Look for files matching the specified criteria
      self._file_matcher = MpfxFileMatcher(master, self, file_search_pattern,
                                           dir_list, image_list,
                                           image_range, sort, recurse)

      _ = self._file_matcher.match_files(master)

      # --- Here we are interested in having a set of file paths per image, each path corresponding
      #     to a type of image (e.g. image-xxx-00.fits or galaxy_catalog-xxx.fits)
      #     So we use the file path dictionary from the MpfxFileMatcher class.

      return self._file_matcher.filepath_dico

   # -----------------------------------------------------------------------------------------------
   def match_file(self, master, path, filename, pattern):
      """!
         File matching predicate method. Must return True in case of matching, False otherwise.
         May be overriden by subclasses to set additional criteria.
         @param master master object instance
         @param path path of the file (without the file name)
         @param filename name of the file
         @param pattern search pattern
         @retval True if the file with name @c filename and path @c path matches @c pattern
         @retval False if no match
      """

      return self._file_matcher.match_file(master, path, filename, pattern)

   # -----------------------------------------------------------------------------------------------
   def format_filename(self, master, path, filename):
      """!
         Reformat the filename to make sure it has a MPF compliant format
         May be overriden by subclasses to set additional criteria.
         @param master master object instance
         @param path path of the file (without the file name)
         @param filename name of the file
         @retval True if the file with name @c filename and path @c path matches @c pattern
         @retval False if no match
      """

      return filename

   # -----------------------------------------------------------------------------------------------
   def get_filepath_dico(self):
      """!
         Return the dictionary of matched file paths
         @return the dictionary of matched file paths
      """
      return self._file_matcher.filepath_dico

   # -----------------------------------------------------------------------------------------------
   def is_fits(self, master, file_path, hdu_no=0):
      """!
         Tell if a file has a FITS format
         @param master master object instance
         @param file_path full path of the file
         @param hdu_no optional HDU number (default 0)
         @retval True if the file has FITS format
         @retval False if not
      """

      is_fits = True

      try:
         # --- Check file header against a .FITS header
         getheader(file_path, hdu_no)

      except Exception:
         # --- Malformed header or not a fits file
         is_fits = False

      return is_fits

   # -----------------------------------------------------------------------------------------------
   def is_catalog(self, master, file_path):
      """! Attempt to tell whether a file is a catalog
           @param master master object instance
           @param file_path absolute file path of a file
           @note The behavior of this method may be overriden
                  to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_catalog = False
      dataset_section_name = "PRIMARY_DATASET.CATALOG_PROPERTIES"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            # Check if a key exists common to stars and galaxies
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            is_catalog = self.is_star_catalog(master, file_path) or \
                         self.is_galaxy_catalog(master, file_path)
      else:
         self.helper.print_error("is_catalog() - Section [{0}] not found".format(
                                                                            dataset_section_name))


      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_image(self, master, file_path):
      """! Attempt to tell whether a file is an image
           @param master master object instance
           @param file_path absolute file path of a file
           @note The behavior of this method may be overridden
                  to comply with the naming conventions of a specific source dataset
      """
      _, filename = os.path.split(file_path)

      is_image = False
      dataset_section_name = "PRIMARY_DATASET.IMAGE_PROPERTIES"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            # Check if a key exists common to stars and galaxies
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            is_image = self.is_star_image(master, file_path) or \
                       self.is_galaxy_image(master, file_path)
      else:
         self.helper.print_error("is_image() - Section [{0}] not found".format(
                                                                            dataset_section_name))

      return is_image


   # -----------------------------------------------------------------------------------------------
   def is_star_catalog(self, master, file_path,
                             ):
      """! Attempt to determine if a file is a star catalog
           @param master master object instance
           @param file_path absolute file path of a file
           @note  to change this behavior, override this method in a sub-class of MPfx
      """

      _, filename = os.path.split(file_path)

      is_catalog = False
      dataset_section_name = "PRIMARY_DATASET.CATALOG_PROPERTIES.STAR"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
      else:
         self.helper.print_error("is_star_catalog() - Section [{0}] not found".format(
                                                                            dataset_section_name))


      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_catalog(self, master, file_path):
      """! Attempt to determine if a file is a galaxy catalog
           @param master master object instance
           @param file_path absolute file path of a file (image, catalog, etc.)
           @note  to change this behavior, override this method in a sub-class of MPfx
      """

      _, filename = os.path.split(file_path)

      is_catalog = False
      dataset_section_name = "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list(
                                                      "FILENAME_PATTERNS", dataset_section_name)
            is_catalog = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            self.helper.print_error("is_galaxy_catalog() - Key: {0}.{1} not found".format(
                                                      dataset_section_name, "FILENAME_PATTERNS"))
      else:
         self.helper.print_error("is_galaxy_catalog() - Section [{0}] not found".format(
                                                                            dataset_section_name))

      return is_catalog

   # -----------------------------------------------------------------------------------------------
   def is_star_image(self, master, file_path):
      """! Attempt to determine if a file is a star image
           @param master master object instance
           @param file_path absolute file path of a file
           @note  the behavior of this method should probably be overridden in a subclass of Mpfx
                  in order to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_image = False
      dataset_section_name = "PRIMARY_DATASET.IMAGE_PROPERTIES.STAR"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            self.helper.print_error("is_star_image() - Key: {0}.{1} not found".format(
                                                      dataset_section_name, "FILENAME_PATTERNS"))

      else:
         self.helper.print_error("is_star_image() - Section [{0}] not found".format(
                                                                            dataset_section_name))
      return is_image

   # -----------------------------------------------------------------------------------------------
   def is_galaxy_image(self, master, file_path,
                             ):
      """! Attempt to determine if a file is a star image
           @param master master object instance
           @param file_path absolute file path of a file
           @note  the behavior of this method should probably be overridden in a subclass of Mpfx
                  in order to comply with the naming conventions of a specific source dataset
      """

      _, filename = os.path.split(file_path)

      is_image = False
      dataset_section_name = "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY"
      if master.config.has_section(dataset_section_name):
         if master.config.has_key("FILENAME_PATTERNS", dataset_section_name):
            filename_patterns = master.config.get_as_list("FILENAME_PATTERNS", dataset_section_name)
            is_image = numpy.any([glob.fnmatch.fnmatch(filename, p) for p in filename_patterns])
         else:
            self.helper.print_error("is_galaxy_image() - Key: {0}.{1} not found".format(
                                                      dataset_section_name, "FILENAME_PATTERNS"))
      else:
         self.helper.print_error("is_galaxy_image() - Section [{0}] not found".format(
                                                                            dataset_section_name))

      return is_image


   # ~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _is_fits_image(self, file_path, hdu_no=0):
      """! Currently Unused """

      is_fits_image = False

      try:
         # --- Check file header agains a .FITS header
         image_header = getheader(file_path, hdu_no)
         is_fits_image = ("NAXIS" in image_header and image_header["NAXIS"] > 0)
      except Exception:
         # --- Malformed header or not a fits file, assume not a .FITS image
         is_fits_image = False

      return is_fits_image

   # -----------------------------------------------------------------------------------------------
   def _is_fits_table(self, file_path, hdu_no=1):
      """! Currently Unused """

      is_fits_catalog = False

      # --- Check for a .FITS table
      try:
         # --- For for a table extension in the first HDU (not the primary)
         table_header = getheader(file_path, hdu_no)
         if "XTENSION" in table_header:
            is_fits_catalog = ("BINTABLE" in table_header["XTENSION"] or \
                               "TABLE" in table_header["XTENSION"])
         else:
            # --- Possibly an image extension
            is_fits_catalog = False

      except IndexError :
         # --- bad HDU, assume not a .FITS catalog
         is_fits_catalog = False

      return is_fits_catalog


# --------------------------------------------------------------------------------------------------
class MpfxFileMatcher(object):

   """!
        The responsibility of this class is to locate all files in the source dataset
        that match specified criteria (directory list, file pattern, etc.)
        @note: this should be an inner class but this creates problems when pickling the object
   """

   def __init__(self, master, dataset, file_patterns, dir_list,
                      image_list, image_range,
                      sort=True, recurse=True):
      """!
         Construct a FileMatcher class
         @param master master object instance
         @param dataset the Dataset object that instantiated this class
         @param file_patterns a list of Unix-style file matching patterns, like
                 ["*.fts", "*.txt"]
         @param dir_list [optional] a list of specific directories to search under base directory
         @param image_list a list of specific image numbers to search, like [0, 100]
         @param image_range a range of specific numbers images to search like [50, 100]
         @param sort [optional] tell whether to sort the output file names (default @c True)
         @param recurse [optional] tell whether to walk down directories (default @c True)
      """

      self._dataset = dataset             # owner dataset
      self._file_patterns = file_patterns # file search pattern (Unix style)
      self._dir_list = dir_list           # list of directories to search under base dir (opt.)
      self._image_list = image_list       # a list of specific image numbers to search
      self._image_range = image_range     # a range of specific images numbers to search
      self._sort = sort                   # sort results if True
      self._recurse = recurse             # walk down directories if True

      self._filepath_dico = {}            # dictionary with key: img_no, value: G3 files/type

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def dataset(self):
      """! @return the associated dataset object """
      return self._dataset

   @property
   def file_patterns(self):
      """! @return the pattern to match with """
      return self._file_patterns

   @property
   def dir_list(self):
      """! @return the optional directory list to search for matching files """
      return self._dir_list

   @property
   def image_list(self):
      """! @return the image_list to match with """
      return self._image_list

   @property
   def image_range(self):
      """! @return the image_range to match with """
      return self._image_range

   @property
   def sort(self):
      """!
         Tells whether file paths should be sorted or not
         @retval True file path sorted
         @retval False otherwise
      """
      return self._sort

   @property
   def recurse(self):
      """!
         Tells if directories are recursively searched
         @retval True directories recursively searched
         @retval False otherwise
      """
      return self._recurse

   @property
   def filepath_dico(self):
      """!
          @return a dictionary in with keys are inmages and values are a dictionary of
          file paths for each image
      """
      return self._filepath_dico

   # --- Setters

   @file_patterns.setter
   def file_patterns(self, file_patterns):
      """!
         Specify the pattern to match with
         @param file_patterns list of Unix-style file matching patterns
      """
      self._file_patterns = file_patterns

   @dir_list.setter
   def dir_list(self, dir_list):
      """!
         Specify the optional directory list to search for matching files
         @param dir_list a list of directory names
      """
      self._dir_list = dir_list

   @image_list.setter
   def image_list(self, image_list):
      """!
         Specify the image list to match with
         @param image_list a list of specific image numbers to search, like [0, 100]
      """
      self._image_list = image_list

   @image_range.setter
   def image_range(self, image_range):
      """!
         Specify the image range to match with
         @param image_range a range of specific numbers images to search like [50, 100]
      """
      self._image_range = image_range

   @sort.setter
   def sort(self, sort):
      """!
         Sort file paths if set to true
         @param sort tell whether to sort the output file names (@c True) or not (@c False)
      """
      self._sort = sort

   @recurse.setter
   def recurse(self, recurse):
      """!
         Walk down directories if set to True
         @param recurse tell whether to recursively search directories (@c True) or not (@c False)
      """
      self._recurse = recurse

   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   # --------------------------------------------------------------------------------------------
   def match_files(self, master):
      """!
          Return the list of files matching the specified search criteria.
          May be overriden by subclasses to set additional criteria.
          @return list of files matching the search criteria
      """

      # --- Build  a list of paths to search for file with the specified pattern
      search_paths = self._build_search_paths(master)

      # --- Locate the files
      file_paths = []
      for search_path in search_paths:
         file_paths.extend(self.dataset.locate_files(master, self.file_patterns, search_path,
                                                     sort=self.sort, recurse=self.recurse,
                                                     err_check=True) )

      #print "*** MPFX: filepaths:", file_paths

      return file_paths     # set of matching files (without duplicates)

   # --------------------------------------------------------------------------------------------
   def match_file(self, master, path, filename, pattern):
      """!
         File matching predicate method, called for each file found in the source dataset.
         Must return @c True in case of matching, @c False otherwise.

         @param master master object instance
         @param path path of the file (without the file name)
         @param filename name of the file
         @param pattern search pattern
         @retval True if the file with name @c filename and path @c path matches @c pattern
         @retval False if no match
         @note May be overriden by subclasses to set additional criteria.
      """

      #print "*** MPFG: match_file() input path:", path

      # --- Reformat the filename to make sure it has a MPF compliant format
      filename = self.dataset.format_filename(master, path, filename)

      #print "*** MPFX: input path 2:", path, filename

      # --- Match filename with specified Unix pattern
      matched = glob.fnmatch.fnmatch(filename, pattern)
      if matched:

         # --- Check image range and list if specified
         if self.image_range == -1 or self.image_range == [-1]:
            self.image_range = [-1, -1]  # in case -1 was specified instead of [-1, -1]

         if len(self.image_list) >= 0 or self.image_range != [-1, -1]:

            regex_dash = re.compile("\-(.*[0-9])\-(.*[0-9])\.")  # xxx-000-00.* (with epoch index)
            regex_dot  = re.compile("\-(.*[0-9])\.")             # xxx-000.*    (without epoch index)

            (img_no, epoch_no) = self._get_img_no_epoch(filename, regex_dash, regex_dot)

            if img_no is not None:
               # --- First look at the image number list: has priority over any range specified
               if len(self.image_list) > 0:

                  matched = img_no in self.image_list
               else:

                  # --- Then look for a non-empty image range
                  if len(self.image_range) > 0:

                     if len(self.image_range) == 1 or len(self.image_range) > 2:

                        self.dataset.helper.print_error(
                              "Image range {0} must be specified as "\
                              "IMAGE_RANGE=[min img_no, max img_no]".format(self.image_range))

                     else:
                        if self.image_range != [-1, -1]:

                           [min_no, max_no] = self.image_range

                           if min_no == -1:
                              matched = (img_no <= max_no)
                           elif max_no == -1:
                              matched = (img_no >= min_no)
                           else:
                              matched = (img_no >= min_no and img_no <= max_no)

               if matched:

                  #print "*** MFGX: calling build_filepath_dico(() with:", self._filepath_dico, path, filename, img_no

                  # --- Populate filepath dictionary per image (with unformatted file paths)
                  self.build_filepath_dico(self._filepath_dico, path, filename, img_no, epoch_no)

                  #print "==> FINAL filepath_dico:", self.filepath_dico, "\n"

            else:
               self.dataset.helper.print_warning(
                     "File: {0}  does not have the required format "\
                     "<img_name-img_no-epoch.ext> => Ignored".format(filename))

      return matched


   # --------------------------------------------------------------------------------------------
   def build_filepath_dico(self, filepath_dico, path, filename, img_no, epoch_index):
      """!
         Build a dictionary with (key=@c tuple(path_tree, img_no, epoch) and
         value=@c filepath for each tuple (file_type, extension)
      """

      #print "*** MPFX base dir:", self.dataset.base_dir

      base_dir = self.dataset.base_dir
      if base_dir.endswith("/"):
         base_dir = base_dir.rstrip("/")

      if len(self.dataset.dir_list) > 0:
         tree_list = path.split(base_dir+"/")[1].split('/')
      else:
         tree_list = ["."]

      tree_list.extend([img_no])

      #print "*** MPFX tree list:", tree_list

      (file_prefix, ext) = os.path.splitext(filename)
      if '-' in file_prefix:
         # Extracts the type of image e.g. "image in "image-000-0-fits"
         file_components = file_prefix.rsplit('-', 2)


         if not epoch_index is None:
            tree_list.extend([epoch_index])  # include epoch index in the key
            file_type = "".join(file_components[:-2])
         else:
            tree_list.extend([0])   # no epoch index, take 0 => dictionary remains the same
            file_type = "".join(file_components[:-1])

         # Create a file_path_dico dict. with key: (path_tree, img_no, epoch)
         key = tuple(tree_list)
         if not key in filepath_dico:
            filepath_dico[key] = {}

         # Record filepath, keeping track of multiple file extensions for catalogs
         # (e.g. -fits and .txt for the same file type...)
         filepath_dico[key][file_type+ext] = os.path.abspath(os.path.join(path, filename))

      #print "MPFX: filepath_dico:", filepath_dico

      #print "*** MPFGX tree list:", filepath_dico

   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   # --------------------------------------------------------------------------------------------
   def _build_search_paths(self, master):
      """! Build  a list of paths to search for file with the specified pattern """

      search_paths = []

      if len(self.dir_list) > 0:
         for dir_name in self.dir_list:
            search_paths.append(os.path.join(self.dataset.base_dir, dir_name))
      else:
         search_paths.append(self.dataset.base_dir)

      return search_paths

   # --------------------------------------------------------------------------------------------
   def _get_img_no_epoch(self, filename, regex_dash, regex_dot):
      """! Extract the image and epoch numbers from the filename """

      img_no = None
      epoch = None

      res_dash = regex_dash.search(filename)
      res_dot  = regex_dot.search(filename)

      #print "MPFX _get_img_no_epoch() filename:", filename

      if not res_dash is None:
         file_main, file_ext = os.path.splitext(filename)

         img_num = res_dash.group(1)
         epoch_num = res_dash.group(2)

         #print "img_num:", img_num, "epoch_num:", epoch_num

         if img_num.isdigit() and epoch_num.isdigit():
            img_no = int(img_num)   # image number from: xxx-000-00.*
            epoch  = int(epoch_num)    # epoch index from: xxx-000-00.*
         else:
            img_num_list = file_main.rsplit('-', 2)
            #print "img_num_list:", img_num_list

            if len(img_num_list) > 2:
               if img_num_list[-1].isdigit() and img_num_list[-2].isdigit():
                  img_no = int(img_num_list[-2])   # image number from: xxx-000-00.*
                  epoch  = int(img_num_list[-1])   # epoch number from: xxx-000-00.*

            elif len(img_num_list) > 1:
               if img_num_list[-1].isdigit():
                  img_no = int(img_num_list[-1])   # image number from: xxx-000-00.*
                  epoch = 0
      else:
         if not res_dot is None:
            img_num = res_dot.group(1)
            if img_num.isdigit():
               img_no = int(img_num)    # image number from: xxx-000.*
               epoch = 0
         else:
            img_num_list = img_num.rsplit('-', 1)
            if img_num_list[-1].isdigit():
               img_no = int(img_num_list[-1])   # image number from: xxx-000.*
               epoch = 0

      return (img_no, epoch)



# -------------------------------------------------------------------------------------------------
class MpfxMultiEpochDataset(MpfxDataset):

   """!
       Represent a source dataset that takes care of files with multiple exposures (epochs).
       This class is based on the "mono epoch" MpfxDataset class.
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
         Query files matching a set of directory and file patterns, taking care of images with
         multiple epochs (exposures). Files with multiple epochs are expected to have format:
         @code <image_name>.<image_no>.<epoch_no>.<extension> @endcode where @c epoch_no represents an epoch
         index starting with zero.

         @param master master object instance
         @param file_search_pattern Unix-style file search pattern list (like "*.fits")
         @param dir_list [optional] a list of specific directories to search under base directory
         @param image_list a list of specific image numbers to search, like [0, 100]
         @param image_range a range of specific numbers images to search like [50, 100]
         @param sort [optional] tell whether to sort the output file names (default @c True)
         @param recurse [optional] tell whether to walk down directories (default @c True)

         @return a list of file paths matching the search criteria
      """

      # --- Look for files matching the specified criteria
      self._file_matcher = MpfxMultiEpochFileMatcher(master, self, file_search_pattern,
                                                             dir_list, image_list,
                                                             image_range, sort, recurse)

      _ = self._file_matcher.match_files(master)

      # --- Here we are interested in having a set of file paths per image, each path corresponding
      #     to a type of image (e.g. image-xxx-00.fits or galaxy_catalog-xxx.fits)
      #     So we use the file path dictionary from the MpfxFileMatcher class.

      return self._file_matcher.filepath_dico



# --------------------------------------------------------------------------------------------------
class MpfxMultiEpochFileMatcher(MpfxFileMatcher):

   """!
        The responsibility of this class is to locate all files in the source dataset
        that match specified criteria (directory list, file pattern, etc.)
        @note: this should be an inner class but this creates problems when pickling the object
   """

   def __init__(self, master, dataset, file_patterns, dir_list,
                      image_list, image_range, sort=True, recurse=True):
      """!
         Construct a FileMatcher class
         @param master master object instance
         @param dataset the Dataset object that instantiated this class
         @param file_patterns a list of Unix-style file matching patterns, like
                 ["*.fits", "*.txt"]
         @param dir_list [optional] a list of specific directories to search under base directory
         @param image_list a list of specific image numbers to search, like [0, 100]
         @param image_range a range of specific numbers images to search like [50, 100]
         @param sort [optional] tell whether to sort the output file names (default @c True)
         @param recurse [optional] tell whether to walk down directories (default @c True)
      """

      MpfxFileMatcher.__init__(self, master, dataset, file_patterns, dir_list,
                                     image_list, image_range, sort, recurse)


   # ~~~~~~~~~~~~~~
   # Public Methods
   # ~~~~~~~~~~~~~~

   # --------------------------------------------------------------------------------------------
   def build_filepath_dico(self, filepath_dico, path, filename, img_no, epoch_index=None):
      """!
         Build a dictionary with (key=@c tuple(path_tree, img_no) and
         value=@c filepath for each tuple (file_type, extension)
      """

      base_dir = self.dataset.base_dir
      if base_dir.endswith("/"):
         base_dir = base_dir.rstrip("/")

      if len(self.dataset.dir_list) > 0:
         tree_list = path.split(base_dir+"/")[1].split('/')
      else:
         tree_list = ["."]

      tree_list.extend([img_no])

      (file_prefix, ext) = os.path.splitext(filename)
      if '-' in file_prefix:
         # Extracts the type of image e.g. "image in "image-000-0-fits"
         file_components = file_prefix.split('-')
         file_type = file_components[0]

         # Create a file_path_dico dict. with key: (path_tree, img_no, epoch)
         key = tuple(tree_list)
         if not key in filepath_dico:
            filepath_dico[key] = {}

         # Record filepath, keeping track of multiple file extensions for catalogs
         # (e.g. -fits and .txt for the same file type...)
         if not file_type+ext in filepath_dico[key]:
            filepath_dico[key][file_type+ext] = []
         filepath_dico[key][file_type+ext].append(os.path.abspath(os.path.join(path, filename)))


# -- EOF mpfx_data.py
