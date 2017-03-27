"""!
   @package pse.pse_helper helper class
   @author Marc Gentile
   @file pse_help.py
   Helper class
"""

# -- Python imports
import os, sys
import numpy
import math
import pyfits
import glob

# -- External import
from scatalog import *                   # catalog management
from mpfg.mp_helper import *             # base Helper
from pse_version import __version__  # version information


# -- Module-specific imports



# -------------------------------------------------------------------------------------------------
class PseHelper(Helper):

   """! Convenient utility functions """

   # -----------------------------------------------------------------------------------------------
   def __init(self):
      Helper.__init__(self)

   # -----------------------------------------------------------------------------------------------
   def get_version(self):
      """!
          Get the version number of this pse code as text
          @return version number of this pse code as text
      """
      return __version__


   # -----------------------------------------------------------------------------------------------
   def show_config_summary(self, master):

      if master.logging_enabled():

         # --- gfit version
         master.logger.log_info_p("\n*** pse {0} ***\n".format(self.get_version()))

         # --- Python modules
         master.logger.log_info_p("Standard Python modules:")
         master.logger.log_info_p("- numpy {0}\t\t{1}".format(numpy.__version__, numpy.__file__))
         master.logger.log_info_p("- pyfits {0}\t\t{1}".format(pyfits.__version__, pyfits.__file__))
         try:
            mpl = __import__('matplotlib', globals(), locals(), [], -1)
            master.logger.log_info_p("- matplotlib {0}\t{1}".format(mpl.__version__, mpl.__file__))
            asc = __import__('asciidata', globals(), locals(), [], -1)
            master.logger.log_info_p("- asciidata {0}\t{1}".format(asc.__version__, asc.__file__))
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be imported: {0}\n".format(detail))

         master.logger.log_info_p("\nMPF Python modules:")
         try:
            mpfg = __import__('mpfg', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfg {0}\t\t{1}".format(mpfg.__version__, mpfg.__file__))
            mpfx = __import__('mpfx', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfx {0}\t\t{1}".format(mpfx.__version__, mpfx.__file__))
            slog = __import__('slogger', globals(), locals(), [], -1)
            master.logger.log_info_p("- slogger {0}\t\t{1}".format(slog.__version__, slog.__file__))
            sconf = __import__('sconfig', globals(), locals(), [], -1)
            master.logger.log_info_p("- sconfig {0}\t\t{1}".format(sconf.__version__, sconf.__file__))
            scat = __import__('scatalog', globals(), locals(), [], -1)
            master.logger.log_info_p("- scatalog {0}\t{1}".format(scat.__version__, scat.__file__))
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be imported: {0}\n".format(detail))

         master.logger.log_info_p("\n")

   # -----------------------------------------------------------------------------------------------
   def extract_stamp_around_centroid(self, xc, yc, half, image):
      """!
         Cut out a square stamp around the specified centroid (@c xc, @c yc).
         @param xc x coordinate of object centroid in postage stamp
         @param yc y coordinate of object centroid in postage stamp
         @param half half of the postage stamp size (of a side)
         @param image the field image from where the postage stamp has to be cut out

         @note stamp size is supposed to be a even number of pixels.
      """

      return image[xc-half+1:xc+half+1, yc-half+1:yc+half+1]

   # -----------------------------------------------------------------------------------------------
   def write_as_fits(self, data, output_filepath, header=None):
      """!
         Write a two-dimensional data numpy array as a .fits file to some given path
         @param data data to write
         @param header header to add to the .FITS file
         @param output_filepath full path of the -FITS file
      """

      if os.path.exists(output_filepath):
         os.remove(output_filepath)

      pyfits.writeto(output_filepath, data, header=header)

   # -----------------------------------------------------------------------------------------------
   def _open_catalog(self, catalog_filepath, job, worker, hdu_no=1):

      catalog = None
      se_output_cat_type = worker.config.get_as_string("SE_OUTPUT_CATALOG_TYPE", "SEXTRACTOR")
      if se_output_cat_type == "FITS_1.0":
         # --- Assume .FITS format
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no,
                                                 open_mode=FITSCatalog.OpenMode.ReadWrite)
      else:
         # --- Assume ASCII SE format
         import inspect
         print 'SF: sexy type', inspect.getfile(SExCatalog)
         catalog = SExCatalog(catalog_filepath)

      if not catalog is None:
         catalog.open()

#     print "OPEN:", "catalog:", catalog, "nb_rows:", catalog.get_nb_rows(), "format:", catalog.format

      return catalog

   # -----------------------------------------------------------------------------------------------
   def _get_catalog_nb_entries(self, catalog_filepath, file_type, job, worker):
      """! Find the actual number of entries in the catalog """

      catalog = self._open_catalog(catalog_filepath, job, worker)
      if not catalog is None:
         entry_count = catalog.get_nb_rows()
         catalog.close()
         return entry_count
      else:
         if worker.logging_enabled():
               worker.logger.log_error_p(
                   "{0} - /{1}/img-{2:03}-{3:1d} - {4} - "\
                  "Error opening catalog {5}".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch,
                   file_type, catalog_filepath) )
         return 0


   # -----------------------------------------------------------------------------------------------
   def _get_nb_objects(self, image_filepath, stamp_size):
      """! Find the actual number of objects (i.e. postage stamps) in the image """

      header = pyfits.getheader(image_filepath)
      nb_pixels = header.get("NAXIS1", 0) * header.get("NAXIS2", 0)
      return int(nb_pixels / stamp_size**2)

   # -----------------------------------------------------------------------------------------------
   def mark_fits_stamps_file(self, fits_path, coordinates, stamp_size,
                                   marking_value, marking_shape="cross"):
      """! Mark fragmented or flagged objects in the check .fits files """

      half_size = self.get_stamp_center(stamp_size)
      check_image = pyfits.getdata(fits_path)   # SE-generated chech .fits image as numpy array

      for (row, col) in coordinates:

         if marking_shape == "cross":
            # --- Draw a cross
            nb_dots = half_size-1    # size of cross

            check_image[row+half_size, col+half_size] = marking_value
            for i in xrange(1, nb_dots+1):
               check_image[row+half_size-i, col+half_size-i] = marking_value
               check_image[row+half_size+i, col+half_size+i] = marking_value
               check_image[row+half_size-i, col+half_size+i] = marking_value
               check_image[row+half_size+i, col+half_size-i] = marking_value

         elif marking_shape == "square":
            # --- Draw a square
            sq_size = half_size-1     # size of square

            check_image[row+half_size-sq_size:row+half_size+sq_size,
                        col+half_size-sq_size:col+half_size-(sq_size-1)] = marking_value
            check_image[row+half_size-sq_size:row+half_size+sq_size,
                        col+half_size+sq_size:col+half_size+sq_size+1] = marking_value
            check_image[row+half_size-sq_size:row+half_size-(sq_size-1),
                        col+half_size-sq_size:col+half_size+sq_size+1] = marking_value
            check_image[row+half_size+sq_size:row+half_size+sq_size+1,
                        col+half_size-sq_size:col+half_size+sq_size+1] = marking_value

      self.write_as_fits(check_image, fits_path)

   # -----------------------------------------------------------------------------------------------
   def mark_centroids(self, fits_path, centroids, marking_value):
      """! Mark centroid locations in the check .fits files """

      check_image = pyfits.getdata(fits_path)   # SE-generated chech .fits image as numpy array
      for (row, col) in centroids:
         check_image[math.floor(row + 0.5), math.floor(col + 0.5)] = marking_value

      self.write_as_fits(check_image, fits_path)

   # -----------------------------------------------------------------------------------------------
   def collect_catalog_data(self, object_per_type_dico, job, worker):
      data_dico = {}

      # --- Forbidden quantities
      forbidden_col_names = ["xc", "yc", "x", "y"]
      forbidden_col_names.extend([
           worker.config.get_as_string("NUMBER_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("CLASS_STAR_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("X_IMAGE_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("Y_IMAGE_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("A_IMAGE_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("B_IMAGE_PARAM", "PARAMETER_MAPPING"),\
           worker.config.get_as_string("FLAGS_PARAM", "PARAMETER_MAPPING") ] )

      # --- Collect data for each file type (image)
      for file_type in object_per_type_dico.keys():

         se_run_dico = object_per_type_dico[file_type]["se_run_dico"]

         # --- Open SE Catalog and collect column data
         se_output_cat_filepath = se_run_dico["se_output_cat_filepath"]

         se_catalog = self._open_catalog(se_output_cat_filepath, job, worker)

         if se_catalog.get_nb_rows() > 0:

            # --- Extra, derived column names
            data_dico[file_type] = {}
            flux_colname    = worker.config.get_as_string("FLUX_PARAM", "PARAMETER_MAPPING")
            fluxerr_colname = worker.config.get_as_string("FLUX_ERR_PARAM", "PARAMETER_MAPPING")

            flux_data = se_catalog.get_named_col_data(flux_colname)
            fluxerr_data = se_catalog.get_named_col_data(fluxerr_colname)
            data_dico[file_type]["SNR"] = numpy.nan_to_num(flux_data / fluxerr_data)

            # --- Allowed SExtractor column data
            for col_name in se_catalog.get_col_names():
               if not col_name in forbidden_col_names:
                  data_dico[file_type][col_name] = se_catalog.get_named_col_data(col_name)

         se_catalog.close()

      return data_dico

   # -----------------------------------------------------------------------------------------------
   def compute_stats(self, object_per_type_dico, data_dico, job, master):
      """! Compute statistical operations on the SExtractor catalog data """

      stats_dico = {}

      if len(data_dico) > 0:

         # --- Statistics to calculate
         stats_oper = [('Len','len'),('Sum','numpy.sum'),('Avg','numpy.mean'),('Med','numpy.median'),\
                       ('Stdev','numpy.std'),('Min','numpy.min'),('Max','numpy.max')]
         stats_oper_dict = dict(stats_oper)
         stats_oper_keys = stats_oper_dict.keys()
         stats_oper_keys.sort()

         # --- Collect data for each file type (image)
         for file_type in object_per_type_dico.keys():

            # --- Compute stats for variables recorded in dictionary
            var_names = data_dico[file_type].keys()
            var_names.sort()

            stats_dico[file_type] = {}
            for var_name in var_names:
               stats_data = numpy.asarray(data_dico[file_type][var_name])
               oper_result_dico = {}
               for oper in stats_oper_keys:
                  o_ptr = eval(stats_oper_dict[oper])
                  oper_result_dico[oper] = o_ptr(stats_data)
               stats_dico[file_type][var_name] = oper_result_dico

      return stats_dico

   # -----------------------------------------------------------------------------------------------
   def make_stats(self, job_result, stat_output_dir, master, stat_prefix=""):

      if len(job_result.stats_dico) > 0:

         job = job_result.job                # associated job
         object_per_type_dico = job_result.result

         file_types = object_per_type_dico.keys()   # all results for each image types

         for file_type in file_types:

            object_dico = object_per_type_dico[file_type]

            # --- Where plots will be stored
            branch_tree = job.get_branch_tree()

            if master.logging_enabled():
               master.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Generating statistics...".format(
                 master, branch_tree, job.img_no, job.epoch, file_type))

            # --- Stats catalog
            stats_filename = stat_prefix + \
                             "stats_{0}_img_{1:03d}-{2:1d}_{3}.txt".format(
                "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)
            stats_filepath = os.path.join(stat_output_dir, stats_filename)

            # --- Produce stats catalogs...
            fd = open(stats_filepath, "w+")

            fd.write('\n')
            title = "STATISTICS for: Image {0:03d}-{1:1d} - /{2} - {3}\n".format(
                                               job.img_no, job.epoch, branch_tree, file_type)
            fd.write(title)
            fd.write("=" * len(title) + "\n\n")


            for col_name in sorted(job_result.stats_dico[file_type].keys()):

               col_stats_dico = job_result.stats_dico[file_type][col_name]

               fd.write('--- Statistics on {0} ---\n'.format(col_name))
               oper_keys = col_stats_dico.keys()
               #fd.write('{0}:\t'.format(col_name))
               for oper in oper_keys:
                  fd.write('{0}={1:.6f}\t'.format(oper, col_stats_dico[oper]))
               fd.write('\n\n')

            fd.close()

   # -----------------------------------------------------------------------------------------------
   def get_stamp_center(self, stamp_size):
      """ ! @return the postage stamp geometrical center pixel no, indexed from zero """

      if stamp_size % 2 == 0:
         # Even size
         center_pixel = int(stamp_size / 2.0 - 1.0)   # e.g. 48x48 -> 23x23
      else:
         # Odd size
         center_pixel = int(stamp_size/ 2.0)          # e.g. 47x47 -> 23x23

      return center_pixel

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

# -- EOF pse_helper.py
