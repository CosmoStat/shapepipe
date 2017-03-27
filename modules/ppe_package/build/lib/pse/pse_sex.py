"""! 
   pse_sex.py - SEXtractor catalog generation and processing
"""

# -- Python imports
import sys, os
import time
import numpy
import copy
import shutil

# -- External imports
from scatalog import *              # catalog management

# --- Module-specific imports
from pse_help import *              # helper utility functions
from pse_plot import *              # plotter 
from pse_counter import *           # counter class 
from operator import itemgetter

# -------------------------------------------------------------------------------------------------
class SExtractorRunner(object):
   
   """! 
      Run SExtractor and create a SExtractor catalog from an input image
   """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, job_processor):

      self._job_processor = job_processor    # the underlying JobProcessor object
      self._helper  = job_processor.helper   # utility class for convenience

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def job_processor(self):
      """! @return the JobProcessor object. """
      return self._job_processor

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def run_SExtractor(self, file_type, job, worker):

      se_object_dico = {}     # Result dictionary to return populated

      try:

         # --- Keep track of execution time
         start_time = time.clock() 

         # --- Image path to analyze
         image_filepath = job.get_file_path(file_type)

         # --- Get the relevant configuration information
         config = worker.config

         # --- Input directory
         input_dir  = os.path.abspath(worker.base_input_dir)

         # SExtractor configuration
         se_config_filepath = self._get_sex_config_filepath(image_filepath, file_type, 
                                                            job, worker)

         if se_config_filepath is not None:
            if worker.logging_enabled():
               worker.logger.log_info_p(
                   "{0} - /{1}/img-{2:03}-{3:1d} - {4} - "\
                   "Using SExtractor .sex configuration file: {5} ...".format(
                     worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                     se_config_filepath))

         # --- Make a copy of the used .sex file to the log directory for record
         shutil.copy(se_config_filepath, worker.log_output_dir)         

         # --- Target directory where to store files
         se_output_path = os.path.join(worker.result_output_dir, job.get_branch_tree()) 

         # --- SExtractor catalog output filepath
         se_output_cat_filename = self._get_output_se_catalog_filename(image_filepath, job, config)

         # --- Output catalog type
         se_output_cat_type = config.get_as_string("SE_OUTPUT_CATALOG_TYPE", "SEXTRACTOR")
         if se_output_cat_type == "FITS_1.0":
            # --- If .fits format specified, force the output filename extension to .fits
            _, fileext = os.path.splitext(se_output_cat_filename)
            if fileext in [".cat", ".txt"]:   
               se_output_cat_filename = se_output_cat_filename.replace(fileext, ".fits")
         elif se_output_cat_type in ["ASCII_SKYCAT", "ASCII_VOTABLE", "FITS_LDAC"]:
            if worker.logging_enabled():
               worker.logger.log_warning_p(
                  "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Output catalog formats 'ASCII SKYCAT', "\
                  "'ASCII VOTABLE' and 'FITS_LDAC' are not supported - "\
                  "Using ASCII_HEAD format...".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
               worker.logger.flush()
            se_output_cat_type = "ASCII_HEAD"   
         else:
            se_output_cat_type = "ASCII_HEAD"   

         se_output_cat_filepath = os.path.abspath(os.path.join(
                                                          se_output_path, se_output_cat_filename))

         if worker.logging_enabled():
            worker.logger.log_info_p(
                "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Generating SExtractor catalog {5} ...".format(
                worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                se_output_cat_filepath))
            worker.logger.flush()

         # --- SExtractor check .fits output filepath
         se_output_check_filename = self._get_output_se_check_filename(image_filepath, job, config)
         se_output_check_filepath = os.path.abspath(os.path.join(se_output_path, 
                                                                 se_output_check_filename))

         # --- SExtractor redir .txt output filepath
         se_must_redir_output = config.get_as_boolean("SE_REDIR_OUTPUT", "SEXTRACTOR")
         se_output_redir_filename = self._get_output_se_redir_filename(image_filepath, job, config)
         se_output_redir_filepath = os.path.abspath(os.path.join(se_output_path, 
                                                                 se_output_redir_filename))

         # --- Check image type
         se_check_image_type = config.get_as_string("SE_CHECKIMAGE_TYPE", "SEXTRACTOR")
         if se_check_image_type != "NONE":
            if worker.logging_enabled():
               worker.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Generating 'Check Image' file {5}...".format(
                 worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                 se_output_check_filepath))
   
         # --- Verbose type
         se_verbose_type = config.get_as_string("SE_VERBOSE_TYPE", "SEXTRACTOR")
            
         # --- SExtractor Execution line
         se_exec_path = config.get_as_string("SE_EXEC_PATH", "SEXTRACTOR")
         if se_must_redir_output:
            se_exec_line =\
             "{0} {1} -c {2} -CATALOG_NAME {3} -CHECKIMAGE_NAME {4} -CHECKIMAGE_TYPE {5} "\
             "-VERBOSE_TYPE {6} -CATALOG_TYPE {7} > {8}".format(
               se_exec_path, image_filepath, se_config_filepath, se_output_cat_filepath, 
               se_output_check_filepath, se_check_image_type, se_verbose_type, 
               se_output_cat_type, se_output_redir_filepath)
         else:
            se_exec_line =\
            "{0} {1} -c {2} -CATALOG_NAME {3} -CHECKIMAGE_NAME {4} -CHECKIMAGE_TYPE {5} "\
            "-VERBOSE_TYPE {6} -CATALOG_TYPE {7} ".format(
            se_exec_path, image_filepath, se_config_filepath, se_output_cat_filepath, 
            se_output_check_filepath, se_check_image_type, se_verbose_type, se_output_cat_type) 
      
         # --- Execute SEXtractor
         cur_dir = os.getcwd()
         os.chdir(input_dir)
         os.system(se_exec_line)
         os.chdir(cur_dir)

         # --- Dictionary with the relevant information for later processing
         se_object_dico = {}
         se_object_dico["se_output_cat_filepath"] = se_output_cat_filepath
         se_object_dico["se_output_check_filepath"] = se_output_check_filepath
         se_object_dico["se_output_redir_filepath"] = se_output_redir_filepath
         se_object_dico["se_check_image_type"] = se_check_image_type
         se_object_dico["elapsed_time"] = time.clock() - start_time

         if worker.logging_enabled():
            if os.path.exists(se_output_cat_filepath):
             worker.logger.log_info_p(
               "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Catalog {5} generated successfully".format(
                worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                  se_output_cat_filepath))
            else:
               worker.logger.log_error_p(
                 "{0} - An error occurred while generating SEXtractor catalog: {1}".format(
                                                                              worker.name, job))
            worker.logger.flush()

      except Exception as detail:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                "{0} - An error occurred while generating SEXtractor catalog: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                detail))
      return se_object_dico


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _get_sex_config_filepath(self, image_file_path, file_type, job, worker):
      """! 
         Find and return the .sex SExtrator cofiguration filepath to use. 
         - First look for a .sex file having the same filename as that of the image file in the 
           @c input directory. 
         - If no file is found, look for a file with the name specified by the
           @c SEXTRACTOR.SE_DEFAULT_SEX_FILENAME key in the configuration file. 
         - If still not found, look for that file in the parent directories. 
         - If the search fails, return @c None.
      """

      sex_filepath = None

      # --- Look for a .sex file that complies with the same filename as that of the image but 
      #     with a -sex extension
      file_dir, file_name = os.path.split(image_file_path)
      file_main, file_ext = os.path.splitext(file_name)

      search_path = os.path.abspath(os.path.join(worker.base_input_dir, job.branch))
      sex_file_candidate = os.path.join(search_path, file_name.replace(file_ext, ".sex"))

      if self.helper.file_exists(sex_file_candidate):
         # --- Found a .sex file for this image
         sex_filepath = sex_file_candidate
      else:
         # --- Look for a .sex file having the same file_type
         sex_file_candidate = os.path.join(search_path, file_type.replace(file_ext, ".sex"))
          
         if self.helper.file_exists(sex_file_candidate):
            # --- Found .sex file applicable to all images with the same filetype
            sex_filepath = sex_file_candidate
         else:
            # --- Looking for a default .sex file somehwre in the base input directory tree
            default_sex_filename = worker.config.get_as_string("SE_DEFAULT_SEX_FILENAME",
                                                               "SEXTRACTOR")
            found_files = self.helper.locate_files(worker, 
                                                   [default_sex_filename], worker.base_input_dir)
            if len(found_files) > 0:
               sex_filepath = found_files[0]
            else:
               sex_filepath = os.path.abspath(os.path.join(worker.base_input_dir, "default.sex"))

               if worker.logging_enabled():
                  worker.logger.log_warning_p(
                      "{0} - /{1}/img-{2:03}-{3:1d} - {4} - "\
                      "Could not find Sextractor .sex file {5}- Trying {6}".format(
                        worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                        default_sex_filename, sex_filepath))

      return sex_filepath   


   # -----------------------------------------------------------------------------------------------
   def _get_output_se_catalog_filename(self, image_filepath, job, config):
      """! 
         Build the SE output catalog filename based on the image filepath  
      """
      _, image_filename = os.path.split(image_filepath)
      se_catalog_prefix, _ = os.path.splitext(image_filename)
      se_catalog_filename_pattern = config.get_as_string("SE_OUTPUT_CATALOG_FILENAME", "SEXTRACTOR")
      return se_catalog_filename_pattern.format(se_catalog_prefix)
                     
   # -----------------------------------------------------------------------------------------------
   def _get_output_se_check_filename(self, image_filepath, job, config):
      """! 
         Build the SE check .fits filename based on the image filepath  
      """
      _, image_filename = os.path.split(image_filepath)
      se_check_prefix, _ = os.path.splitext(image_filename)
      se_check_filename_pattern = config.get_as_string("SE_OUTPUT_CHECK_FILENAME", "SEXTRACTOR")
      return se_check_filename_pattern.format(se_check_prefix)

   # -----------------------------------------------------------------------------------------------
   def _get_output_se_redir_filename(self, image_filepath, job, config):
      """! 
         Build the SE redir .txt filename based on the image filepath  
      """
      _, image_filename = os.path.split(image_filepath)
      se_redir_prefix, _ = os.path.splitext(image_filename)
      se_redir_filename_pattern = config.get_as_string("SE_OUTPUT_REDIR_FILENAME", "SEXTRACTOR")
      return se_redir_filename_pattern.format(se_redir_prefix)



# -------------------------------------------------------------------------------------------------
class SExtractorProcessor(object):
   
   """! 
      Process a SExtractor catalog
   """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, job_processor):

      self._job_processor = job_processor    # the underlying JobProcessor object
      self._helper = PseHelper()             # ref. to utility class for convenience
      self._plotter = job_processor.plotter  # utility class for convenience   

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def plotter(self):
      """! @return the Helper class """
      return self._plotter

   @property
   def job_processor(self):
      """! @return the JobProcessor object. """
      return self._job_processor

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def process_catalog(self, se_catalog_filepath, se_output_check_filepath, file_type, job, worker):
      """!
         Apply various transformations to a SEXtractor catalog
         @param se_catalog_filepath full filepath of the catalog
         @param se_output_check_filepath full path og check file to generate (if requested)
         @param file_type type of file to process
         @param job underlying Job object   
         @param worker associated worker process 
      """   

      se_xform_dico = {}   # result dictionary to populate

      try:

         se_xform_dico = {}         # result dictionary to populate
         deleted_coords = []

         start_time = time.clock()  # Keep track of execution time
         
         # --- Image and associated catalog paths
         image_filepath = job.get_file_path(file_type) 

         # --- Actual number of objects in the image
         initial_object_count = self.helper._get_catalog_nb_entries(se_catalog_filepath, 
                                                                    file_type, job, worker)
             
         if initial_object_count > 0:    
             
            # --- Add extra columns to the catalog
            actual_object_count = self._add_extra_colunmns(se_catalog_filepath, file_type,
                                     se_output_check_filepath, job, worker)

            if actual_object_count > 0:
   
               # --- Handle splitted objects, i.e. multiple detection of the same object
               duplicate_object_count, duplicated_coords = self._handle_splitted_objects(
                                                                      se_catalog_filepath,
                                                                      se_output_check_filepath, 
                                                                      file_type, job, worker) 
               deleted_coords.extend(duplicated_coords)

               # --- Filter out flagged entries
               flagged_object_count, deleted_flagged_coords = self._filter_flagged_entries(
                                                                   se_catalog_filepath, file_type,
                                                                   se_output_check_filepath, 
                                                                   job, worker)
               deleted_coords.extend(deleted_flagged_coords)
      
      #         # --- Filter out entries with excessive centroid shifts
      #         bad_centroid_shift_count, deleted_shifts = self._filter_excessive_centroid_shifts(
      #                                                                          se_catalog_filepath,
      #                                                                          se_output_check_filepath,  
      #                                                                          file_type, job, worker)
      #         deleted_coords.extend(deleted_shifts)
               deleted_coords = list(set(deleted_coords))
      
      
               # --- Sort by (x,y) in ascending order
               final_object_count = self._sort_by_coordinates(se_catalog_filepath, 
                                                              file_type, job, worker)
      
               if final_object_count > 0:
                  # --- If requested, create a mosaic containing the remaining (non deleted) objects 
                  #     referenced in the final SExtractor catalog. The relative positions of the objects
                  #     are maintained. 
                  if worker.config.get_as_boolean("CREATE_STAMP_MOSAIC", "DEBUGGING"):
                     self._create_stamp_mosaic(image_filepath, 
                                               se_catalog_filepath, deleted_coords, 
                                             file_type, job, worker)
      
               #if worker.config.get_as_boolean("CREATE_STAMP_PLOT_MOSAIC", "DEBUGGING"):
               #   self._create_stamp_plot_mosaic(image_filepath, se_catalog_filepath, 
               #                                  file_type, job, worker)
      
               # --- Dictionary with the relevant information for later processing
               #se_xform_dico["actual_object_count"]  = actual_object_count
               se_xform_dico["initial_object_count"] = initial_object_count
               se_xform_dico["duplicate_object_count"] = duplicate_object_count
               se_xform_dico["flagged_object_count"] = flagged_object_count
               #se_xform_dico["bad_centroids_object_count"] = bad_centroid_shift_count
               se_xform_dico["final_object_count"] = final_object_count
               se_xform_dico["elapsed_time"] = time.clock() - start_time
      
               if worker.logging_enabled():
                  if se_xform_dico["final_object_count"] > 0:
                     worker.logger.log_info_p(
                         "{0} - /{1}/img-{2:03}-{3:1d} - {4} - "\
                         "Catalog {5} processed successfully - Final object count: {6}".format(
                         worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                         se_catalog_filepath, se_xform_dico["final_object_count"]))
                  else:
                     worker.logger.log_warning_p(
                         "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Catalog {5} is empty: "\
                         "all {6} entries have been deleted".format(
                         worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                         se_catalog_filepath, initial_object_count))
      
                  worker.logger.flush()
         else:
            # --- No objects were detecte by SExtractor
            if worker.logging_enabled():
               worker.logger.log_error_p(
                   "{0} - /{1}/img-{2:03}-{3:1d} - {4} - "\
                  "No object could be detected by SExtractor".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type) )

      except Exception as detail:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                "{0} - An error occurred while processing SExtractor catalog: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                detail))
            worker.logger.flush()

      return se_xform_dico

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def _add_extra_colunmns(self, se_catalog_filepath, file_type, se_output_check_filepath,  
                                 job, worker):
      """! Append columns to the SE catalog """


      catalog = None
      try:

         # --- Open SE catalog
         catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)

         initial_object_count = catalog.get_nb_rows()

         if initial_object_count > 0:
         
            if worker.logging_enabled():
               worker.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Catalog contains {5} entries - "\
                 "Adding extra columns...".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                      catalog.get_nb_rows())) 
               worker.logger.flush()

            # --- Get coordinate data
            x_image_colname = worker.config.get_as_string("X_IMAGE_PARAM", "PARAMETER_MAPPING")
            y_image_colname = worker.config.get_as_string("Y_IMAGE_PARAM", "PARAMETER_MAPPING")
            xc_array = catalog.get_named_col_data(x_image_colname) - 1.0   # SE coordinates start from 1, not 0
            yc_array = catalog.get_named_col_data(y_image_colname) - 1.0   # SE coordinates start from 1, not 0
            stamp_size = job.get_stamp_size()
            x_array = numpy.floor(xc_array / stamp_size) * stamp_size
            y_array = numpy.floor(yc_array / stamp_size) * stamp_size
      
            rel_xc_array = numpy.floor(xc_array - x_array + 0.5)
            rel_yc_array = numpy.floor(yc_array - y_array + 0.5)
      
      #      half_size = stamp_size / 2.0
      #     #todo: Disabled for now. Must be done in gfit
      #      rel_xc_array = numpy.where(half_size - rel_xc_array <= 1.0, 
      #                                 numpy.ceil(half_size), rel_xc_array)
      #      rel_yc_array = numpy.where(half_size - rel_yc_array <= 1.0, 
      #                                 numpy.ceil(half_size), rel_yc_array)
      
            # --- Append columns
            catalog.add_col("x", col_format="% 5.1f", 
                                 col_comment="stamp x lower-left corner", col_data=x_array)
		            
            catalog.add_col("y", col_format="% 5.1f",
                                 col_comment="stamp y lower-left corner", col_data=y_array)
      
            catalog.add_col("xc", col_format="% 6.3f", 
                                  col_comment="zero-index x centroid", col_data=xc_array)
            catalog.add_col("yc", col_format="% 6.3f", 
                                  col_comment="zero-index y centroid", col_data=yc_array)
      
            catalog.add_col("rel_xc", col_format="% 5.1f", 
                                      col_comment="relative centroid x pixel in stamp", 
                                      col_data=rel_xc_array)
            catalog.add_col("rel_yc", col_format="% 5.1f", 
                                      col_comment="relative centroid y pixel in stamp", 
                                      col_data=rel_yc_array)
      
            # --- Save changes
            catalog.save()
            #initial_object_count = catalog.get_nb_rows()
            #catalog.close()  

      finally:   
         if not catalog is None:
            catalog.close()  

      return initial_object_count


   # -----------------------------------------------------------------------------------------------
   def _mark_centroids(self, se_output_check_filepath, centroids, file_type, job, worker):

      if worker.config.get_as_string("SE_CHECKIMAGE_TYPE", "SEXTRACTOR") != "NONE" and \
         worker.config.get_as_boolean("MARK_CHECK_FITS_FILES", "DEBUGGING"):
         
         # Marking values   
         mark_centroid_pix_value = worker.config.get_as_float("MARKING_VALUE_CENTROIDS", 
                                                              "DEBUGGING")
         if mark_centroid_pix_value != -1:
            worker.logger.log_info_p(
             "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Marking centroids [value={5}]...".format(
                worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                file_type, mark_centroid_pix_value)) 
            self.helper.mark_centroids(se_output_check_filepath, centroids, 
                                                                 mark_centroid_pix_value)

#   # -----------------------------------------------------------------------------------------------
#   def _mark_stamps(self, se_output_check_filepath, coords, file_type, job, worker):

#      if worker.config.get_as_string("SE_CHECKIMAGE_TYPE", "SEXTRACTOR") != "NONE" and \
#         worker.config.get_as_string("MARK_CHECK_FITS_FILES", "DEBUGGING"):
#         
#         # Marking values   
#         mark_stamp_pix_value = worker.config.get_as_float("MARKING_VALUE_STAMPS", "DEBUGGING")
#         if mark_stamp_pix_value != -1:
#            self.helper.mark_fits_stamps_file(se_output_check_filepath, coords,
#                                              job.get_stamp_size(), 
#                                              mark_stamp_pix_value, marking_shape="cross")


   # -----------------------------------------------------------------------------------------------
   def _handle_splitted_objects(self, se_catalog_filepath, se_output_check_filepath,
                                      file_type, job, worker):
      """! Handle splitted objects: delete all duplicates or keep that with highest flux """

      if worker.logging_enabled():
         worker.logger.log_info_p(
                       "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Handling fragmented objects...".format(
                       worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type)) 
         worker.logger.flush()

      # How to handle splitted objects: 
      # - DELETE_ALL: delete all duplicates (i.e. all fragments except centroid)
      # - DELETE_EXTRA: delete the object with lowest flux (centroid assumed to have greatest flux)
      # - KEEP_ALL: do nothing (keep all objects)
      splitted_obj_mgt = worker.config.get_as_string("METHOD", "CLEANING")

      # --- Maximum allowed shift of centroid from geometrical center of postage stamp
      stamp_center = job.get_stamp_center()

      max_centroid_shift = worker.config.get_as_int("MAX_CENTROID_SHIFT", "CLEANING")

      # --- Fragments assessment
      inner_radius = worker.config.get_as_float("INNER_CIRCLE_RADIUS", 
                                                "CLEANING.DELETE_ALL_FILTERED")
      outer_radius = worker.config.get_as_float("OUTER_CIRCLE_RADIUS", 
                                                "CLEANING.DELETE_ALL_FILTERED")
      min_allowed_mag = worker.config.get_as_float("MIN_FRAGMENT_MAG_IN_RING", 
                                                   "CLEANING.DELETE_ALL_FILTERED")
      max_allowed_mag = worker.config.get_as_float("MAX_FRAGMENT_MAG_IN_RING", 
                                                   "CLEANING.DELETE_ALL_FILTERED")

      # --- Open SE catalog
      catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)

      # Initial Nb of entries
      initial_obj_count = catalog.get_nb_rows()

      # Object coordinates
      x_array = catalog.get_named_col_data("x").astype(int)         
      y_array = catalog.get_named_col_data("y").astype(int)         
      coords = zip(y_array, x_array)    # (row,col) numpy convention

      xc_array = catalog.get_named_col_data("xc")
      yc_array = catalog.get_named_col_data("yc")
      centroid_coords = zip(yc_array, xc_array) # (row,col) numpy convention

      rel_xc_array = catalog.get_named_col_data("rel_xc")         
      rel_yc_array = catalog.get_named_col_data("rel_yc")  
      rel_shifts = numpy.hypot(rel_xc_array - stamp_center - 1, rel_yc_array - stamp_center - 1)       

      # --- Catalog information
      cat_col_names = catalog.get_col_names()

      cat_col_comments = catalog.get_col_comments()
      cat_col_formats = catalog.get_col_formats()

      cat_data_matrix = catalog.get_data()

      # --- Locate duplicated objects (which means uncertain detection)
      counter = Counter(coords)
      c_keys   = numpy.asarray(counter.keys())
      c_values = numpy.asarray(counter.values())      

      repeated_coords = [tuple(e) for e in c_keys[c_values > 1]]
      repeated_indice = []
      entry_indice = []

      flux_column_index = catalog.get_col_index(worker.config.get_as_string(
                                                         "FLUX_PARAM", "PARAMETER_MAPPING"))
      mag_column_index  = catalog.get_col_index(worker.config.get_as_string(
                                                         "MAG_PARAM", "PARAMETER_MAPPING"))
      rel_xc_column_index = catalog.get_col_index('rel_xc')
      rel_yc_column_index = catalog.get_col_index('rel_yc')

      #nb_centroids = len(list(set(repeated_coords)))
 
      # --- Collect information about fragmented entries...
      all_object_indice_dico = {}
      extra_object_centroid_dico = {}
      extra_object_flux_dico = {}
      extra_object_mag_dico = {}
      extra_object_shift_dico = {}

      # For collecting info about fragment distance/magnitude relationship, etc.
      fragment_info_dico = {"dist":[], "mag":[], "count":[], \
                            "nb_inner_removed":0, "nb_between_removed":0}   

      for i in xrange(0, len(coords)):

         # --- Check duplicated catalog entries (same stamp coordinates) 
         entry_indice.append(i)
         if not coords[i] in all_object_indice_dico:
            all_object_indice_dico[coords[i]] = []
         all_object_indice_dico[coords[i]].append(i)

         # --- Entries containing at least one fragment 
         if coords[i] in repeated_coords:
            repeated_indice.append(i)

         # --- Store entry fluxes
         if not coords[i] in extra_object_flux_dico:
            extra_object_flux_dico[coords[i]] = []
         object_flux  = cat_data_matrix[i, flux_column_index]
         extra_object_flux_dico[coords[i]].append(object_flux)

         # --- Store entry centroid shifts
         if not coords[i] in extra_object_shift_dico:
            extra_object_shift_dico[coords[i]] = []
         object_shift = (math.fabs(cat_data_matrix[i, rel_xc_column_index] - stamp_center - 1), 
                         math.fabs(cat_data_matrix[i, rel_yc_column_index] - stamp_center - 1))
         extra_object_shift_dico[coords[i]].append(object_shift)

         # --- Store entry magnitudes
         if not coords[i] in extra_object_mag_dico:
            extra_object_mag_dico[coords[i]] = []
         object_mag  = cat_data_matrix[i, mag_column_index]
         extra_object_mag_dico[coords[i]].append(object_mag)

         # --- Store fragment centroids
         if not coords[i] in extra_object_centroid_dico:
            extra_object_centroid_dico[coords[i]] = []
         object_centroid = (cat_data_matrix[i, rel_yc_column_index], 
                            cat_data_matrix[i, rel_xc_column_index])
         extra_object_centroid_dico[coords[i]].append(object_centroid)

      catalog.close()

      if worker.logging_enabled():
         worker.logger.log_info_p(
             "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Found {5} objects with fragments ".format(
             worker.name, job.get_branch_tree(), job.img_no,job.epoch, file_type, 
             len(repeated_coords))) 
         worker.logger.flush()

      # --- Remove some or all duplicated objects
      deleted_indice = []


      count = 0

      if splitted_obj_mgt == 'DELETE_ALL':   

         # --- Delete all duplicated entries in catalog (i.e. unwanted fragments) 
         deleted_indice = repeated_indice

      elif splitted_obj_mgt == 'DELETE_EXTRA' or splitted_obj_mgt == 'DELETE_ALL_FILTERED':  

         # --- If requested, plot fragment-related data for debugging/tuning
         plot_fragment_info = worker.config.get_as_boolean("CREATE_SE_PLOTS", "DEBUGGING")

         # --- Remove objects from catalog depending on the number/mag/distance of fragments
         coords = extra_object_flux_dico.keys()

         for coord in coords:

            # --- Deletes extra objects with smaller flux than peak flux object
            new_indice_list = all_object_indice_dico[coord]

            list_flux  = extra_object_flux_dico[coord]
            list_mag   = extra_object_mag_dico[coord]
            list_shift = extra_object_shift_dico[coord]
            list_centroid = extra_object_centroid_dico[coord]
            
            if max_centroid_shift != -1:
               # Only keep object index with maximum fluxes and suitable maximum centroid shifts 
               # (the other ones will be removed)

               list_shift_flux = [ (s, f) for (s, f) in zip(list_shift, list_flux) \
                                          if f > 0.0 and \
                                             s[0] <= max_centroid_shift and \
                                             s[1] <= max_centroid_shift ]
               
            else:
               list_shift_flux = [ (s, f) for (s, f) in zip(list_shift, list_flux) if f > 0.0 ]
               #print "list_shift_flux:", list_shift_flux

            if len(list_shift_flux) > 0: 

               # --- Flux and index of most likely centroid

               # --- Find the fragment with maximum flux but with the smmallest shift -> centroid
               if len(list_shift_flux) > 1:

                  # Sorted on flux, descending order 
                  sorted_list_shift_flux = sorted(list_shift_flux, key=itemgetter(1, 0), \
                                                                   reverse=True) 

                  max_flux = sorted_list_shift_flux[0][1]   # shifts sorted ascending order
                  min_shift = min([ s for (s, f) in list_shift_flux if f == max_flux ])
                  max_flux_index = list_shift.index(min_shift) 

               else:

                  # --- Only one possible centroid, take it
                  max_flux_index = list_shift.index(list_shift_flux[0][0])

               max_obj_index = new_indice_list[max_flux_index] 

               # --- Delete all entries that do not correspond to the assumed centroid (found above)
               if splitted_obj_mgt == 'DELETE_EXTRA':

                  # --- Add objects with smaller fluxes in the object list to delete
                  new_indice_list.remove(max_obj_index)  # => now contains all indice to delete

                  deleted_indice.extend(new_indice_list)  

               # --- Delete entries provided specific conditions are met
               elif splitted_obj_mgt == 'DELETE_ALL_FILTERED':

                  # --- Object assumed centroid
                  object_centroid = (cat_data_matrix[max_obj_index, rel_yc_column_index], 
                                     cat_data_matrix[max_obj_index, rel_xc_column_index])

                  # --- Only remove all fragments (catalog entries) if specific conditions are met
                  must_delete_object = self._must_delete_object(list_centroid, list_mag, 
                                                                coord, object_centroid, 
                                                                inner_radius, outer_radius, 
                                                                min_allowed_mag, 
                                                                max_allowed_mag, 
                                                                fragment_info_dico, # to populate 
                                                                job, worker)

                  if must_delete_object: 

                     # --- Add rejected entries to the delete list
                     deleted_indice.extend(new_indice_list)  

                  else:

                     # --- Must keep the object in the catalog => retain entry with located centroid
                     new_indice_list.remove(max_obj_index)

                     # --- But nevertheless removes extra entries besides centroid
                     #     => add rejected entries to the delete list
                     deleted_indice.extend(new_indice_list)  
                     
            else:

               #print coord, new_indice_list, list_shift,  
               #print list_flux, list_flux,
               #print "=> *** will remove:", new_indice_list

               list_zero_flux = [ (s, f) for (s, f) in zip(list_shift, list_flux) if f == 0.0 ]
               if len(list_zero_flux) > 0:
                  if worker.logging_enabled():
                     worker.logger.log_warning_p(
                    "{0} - /{1}/img-{2:03}-{3:1d} - {4} - ({5}, {6}) - "\
                     "Object found with zero flux => deleting this object".format(
                     worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                     coord[1], coord[0]))
                     worker.logger.flush()
                  
               else:   
                  # Could not find an object with max flux that has a centroid shift below
                  # <max_centroid_shift> => remove all fragments
                  if worker.logging_enabled():
                     worker.logger.log_warning_p(
                    "{0} - /{1}/img-{2:03}-{3:1d} - {4} - ({5}, {6}) - "\
                     "Could not find maximum-flux fragments with "\
                     "centroid shift below {7} => deleting all {8} fragments "\
                     "for this object".format(
                     worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                     coord[1], coord[0], max_centroid_shift, len(new_indice_list)))
                     worker.logger.flush()

               # --- Add objects with smaller fluxes in the object list to delete
               deleted_indice.extend(new_indice_list)  

         if splitted_obj_mgt == 'DELETE_ALL_FILTERED':

            # --- Print fragment-related stats
            if worker.logging_enabled():

               worker.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Inner circle: removed {5} fragments".format(
                 worker.name, job.get_branch_tree(), job.img_no,job.epoch, file_type, 
                 fragment_info_dico["nb_inner_removed"])) 
               worker.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Ring: removed {5} fragments".format(
                 worker.name, job.get_branch_tree(), job.img_no,job.epoch, file_type, 
                 fragment_info_dico["nb_between_removed"])) 

            # --- Plot fragment-related data
            if plot_fragment_info:
               self.plotter.make_fragment_plots(fragment_info_dico, 
                                                inner_radius, outer_radius, 
                                                min_allowed_mag, max_allowed_mag,
                                                file_type, job, worker)

      # --- Keep all catalog entries (without removing anything)
      elif splitted_obj_mgt == 'KEEP_ALL':
         if worker.logging_enabled():
            worker.logger.log_info_p(
           "{0} - /{1}/img-{2:03}-{3:1d} - {4} - No filtering of fragmented objects".format(
            worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
            worker.logger.flush()

      # --- Actually delete duplicated objects
      splitted_obj_count = 0       

      # --- Final coordinates (numpy comvention )
      remaining_indice = list(set(entry_indice) - set(deleted_indice))
      remaining_coords =  zip(y_array[remaining_indice], x_array[remaining_indice])

      deleted_coords = list(set(coords) - set(remaining_coords))

      deleted_centroid_coords = zip(yc_array[deleted_indice], xc_array[deleted_indice])  
      remaining_centroid_coords = list(set(centroid_coords) - set(deleted_centroid_coords))

      cat_data_matrix = numpy.delete(cat_data_matrix, deleted_indice, 0)

      if worker.logging_enabled():
         worker.logger.log_info_p(
              "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Removed {5} catalog entries".format(
                 worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type,
                 len(deleted_indice)))
         worker.logger.flush()

      # Save changes by overwriting the catalog
      se_output_cat_type = worker.config.get_as_string("SE_OUTPUT_CATALOG_TYPE", "SEXTRACTOR")
      if se_output_cat_type in  ["FITS_1.0"]:
         # --- FITS format
         catalog = FITSCatalog(se_catalog_filepath)
         if len(remaining_indice) > 0:
            
            # Save catalog with remaining (undeleted) entries
            catalog.open()
            header  = catalog.get_header()
            extname = catalog.get_info()["extname"]
            extver  = int(catalog.get_info()["extver"])
            catalog.create_from_numpy(cat_data_matrix, cat_col_names, 
                                      header=header, ext_name=extname, ext_ver=extver)
            
         else:
            # The catalog has become empty
            catalog.create()
      else:
         # --- SExtractor format
         catalog = SExCatalog(se_catalog_filepath)
         if len(remaining_indice) > 0:
            # Save catalog with remaining (undeleted) entries
            
            catalog.create_from_numpy(cat_data_matrix, 
                                      cat_col_names, cat_col_comments, cat_col_formats)
         else:
            # The catalog has become empty
            catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)
            catalog.delete_rows(0, len(deleted_indice))
            
      catalog.save()
      splitted_obj_count = catalog.get_nb_rows()
      catalog.close()          

      if worker.logging_enabled():
         worker.logger.log_info_p(
               "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Remaining number of objects: {5}".format(
                  worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                 len(remaining_centroid_coords)))
         worker.logger.flush()

      # Mark corresponding stamps in check .fits file if requested
      if worker.config.get_as_string("SE_CHECKIMAGE_TYPE", "SEXTRACTOR") != "NONE" and \
         worker.config.get_as_string("MARK_CHECK_FITS_FILES", "DEBUGGING"):

         try:
            # Marking values   
            mark_duplicated_objs_pix_value = worker.config.get_as_float(
                                                   "MARKING_VALUE_DUPLICATED_OBJS", "DEBUGGING")
            if mark_duplicated_objs_pix_value != -1:

               # Do the marking 
               self.helper.mark_fits_stamps_file(se_output_check_filepath, deleted_coords, 
                                                 job.get_stamp_size(), 
                                                 mark_duplicated_objs_pix_value, 
                                                 marking_shape="cross")   

            # Mark centroids in check .FITS file if requested
            deleted_centroid_coords = zip(yc_array[deleted_indice], xc_array[deleted_indice])  
            remaining_centroid_coords = list(set(centroid_coords) - set(deleted_centroid_coords))
            self._mark_centroids(se_output_check_filepath, remaining_centroid_coords, 
                                 file_type, job, worker)

         except Exception as detail:
            if worker.logging_enabled():
               worker.logger.log_error_p(
                "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Error while marking object: {5}".format(
                 worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, detail))
               worker.logger.flush()

      else:
         if worker.logging_enabled():
            worker.logger.log_info_p(
             "{0} - /{1}/img-{2:03}-{3:1d} - {4} - No duplicated objects found to remove...".format(
              worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
            worker.logger.flush()

      return splitted_obj_count, deleted_coords



   # -----------------------------------------------------------------------------------------------
   def _must_delete_object(self, list_centroid, list_mag, object_coord, 
                                 object_centroid,  # (yc, xc), i.e. numpy convention
                                 inner_radius, outer_radius, 
                                 min_allowed_mag, max_allowed_mag, fragment_info_dico,
                                 job, worker):
      """! 
          Decide wether to delete the object (star, galaxy...) based on the characteristics
          of its fragments
      """
      
      must_delete = False
      (yc, xc) = object_centroid
      
      # --- Find distance and magnitudes of fragments (not including assumed centroids)
      dist_mag_pairs = [ (math.hypot(x-xc, y-yc), mag) \
                                      for ((y, x), mag) in zip(list_centroid, list_mag) \
                                      if not (y, x) == (yc, xc)]

      # --- Collect fragment info
      for (dist, mag) in dist_mag_pairs:
 
         fragment_info_dico["dist"].append(dist)
         fragment_info_dico["mag"].append(mag)
         fragment_info_dico["count"].append(len(dist_mag_pairs))

      # --- Selectively remove fragments
      for (dist, mag) in dist_mag_pairs:
 
         # --- Check fragment locations and magnitudes
         if dist <= inner_radius:
            # --- Fragment found inside inner circle, remove this object (i.e. all its 
            #     fragmented entries, including its centroid) from the catalog and exit
            must_delete = True
            #print "coord:", object_coord, "centroid:", object_centroid, "list_centroids:", list_centroid , "list_mags:", list_mag, "dists:", dist_mag_pairs
            #print("*** dist {0} <= inner_radius  {1} => delete".format(dist, inner_radius))
            fragment_info_dico["nb_inner_removed"] += len(dist_mag_pairs) + 1

            #print "*** INNER: Removing fragment:", (dist, mag)
            break       

         elif dist > inner_radius and dist <= outer_radius:
            # --- Fragment located in ring between inner and outer circle => check its magnitude 
            #     is not brighter than <min_allowed_mag>

            if mag >= max_allowed_mag and mag <= min_allowed_mag:
               # --- Magnitude too high, remove object from catalog and exit 
               must_delete = True
               #print "coord:", object_coord, "centroid:", object_centroid, "list_centroids:", list_centroid , "list_mags:", list_mag, "dists:", dist_mag_pairs
               #print("*** mag {0} < min_allowed_mag {1} => delete".format(mag, min_allowed_mag))
               fragment_info_dico["nb_between_removed"] += len(dist_mag_pairs) + 1

               #print "*** RING: Removing fragment:", (dist, mag)
               break       

            else:
               pass

         elif dist > outer_radius:
            pass

      return must_delete

   # -----------------------------------------------------------------------------------------------
   def _filter_flagged_entries(self, se_catalog_filepath, file_type, se_output_check_filepath,  
                                     job, worker):
      """! Remove entries with flags belonging to invalid-flag list """

      if worker.logging_enabled():
         worker.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Filtering out flagged objects...".format(
                 worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type)) 
         worker.logger.flush()

      # --- Open SE catalog
      catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)

      # --- Check for na empty catalog
      if catalog.get_nb_rows() == 0:
         return 0, []

      # --- Catalog information
      cat_col_names = catalog.get_col_names()
      cat_col_comments = catalog.get_col_comments()
      cat_col_formats = catalog.get_col_formats()
      cat_data_matrix = catalog.get_data()

      # --- List of allowed flags
      obj_array = catalog.get_named_col_data(worker.config.get_as_string(
                                                         "NUMBER_PARAM", "PARAMETER_MAPPING"))
      flag_array = catalog.get_named_col_data(worker.config.get_as_string(
                                                         "FLAGS_PARAM", "PARAMETER_MAPPING"))
      obj_flux  = catalog.get_named_col_data(worker.config.get_as_string(
                                                         "FLUX_PARAM", "PARAMETER_MAPPING"))     

      rows = catalog.get_named_col_data("y").astype(int)         
      cols = catalog.get_named_col_data("x").astype(int) 
      coords = zip(rows, cols)            # (row,col) numpy convention

      catalog.close()

      # --- Which flags to take into account
      error_flags = [int(e) for e in worker.config.get_as_list("SE_ERROR_FLAGS", "CLEANING")]
      error_flag_indice = [numpy.where(flag_array == e)[0] for e in error_flags]

      all_flag_indice = numpy.where(flag_array > 0)[0]
      all_flags_found = list(set(flag_array[all_flag_indice])) 
      flagged_coords = zip(rows[all_flag_indice], cols[all_flag_indice])

      if len(flagged_coords) > 0:
         if worker.logging_enabled():
            worker.logger.log_info_p(
             "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Found a total of {5} objects "\
             "with non-zero flags {6}".format(worker.name, job.get_branch_tree(), job.img_no, 
                                   job.epoch, file_type, len(flagged_coords), all_flags_found) ) 
            worker.logger.flush()

      # --- Locate flagged objects
      to_delete = []
      invalid_flags = []
      flagged_coords = []
      for e in error_flag_indice:
         # List of flagged objects to delete
         to_delete = numpy.append(to_delete, e)

         # List of error flags
         invalid_flags = numpy.append(invalid_flags, flag_array[e])

         # List of coordinates (row, col) with error flags
         [flagged_coords.append(coords[i]) for i in e]

      # --- Delete flagged objects
      flagged_object_count = 0

      if len(to_delete) > 0:

         # Mark corresponding stamps in check .fits file if requested
         if worker.config.get_as_string("SE_CHECKIMAGE_TYPE", "SEXTRACTOR") != "NONE" and \
            worker.config.get_as_boolean("MARK_CHECK_FITS_FILES", "DEBUGGING"):
            
            # Marking values   
            mark_flagged_objs_pix_value = worker.config.get_as_float(
                                                         "MARKING_VALUE_FLAGGED_OBJS", "DEBUGGING")
            if mark_flagged_objs_pix_value != -1:
               # Do the marking 
               self.helper.mark_fits_stamps_file(se_output_check_filepath, flagged_coords, 
                                                 job.get_stamp_size(), 
                                                 mark_flagged_objs_pix_value, 
                                                 marking_shape="square")

         # Delete relevant flagged objects
         cat_data_matrix = numpy.delete(cat_data_matrix, to_delete.astype(int), 0)
         if worker.logging_enabled():
            worker.logger.log_info_p(
            "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Removed {5} objects having "\
            "invalid flags: {6}".format(
               worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
               len(to_delete), list(set(invalid_flags)) ))
            worker.logger.flush()

         # Save changes by overwriting the catalog
         se_output_cat_type = worker.config.get_as_string("SE_OUTPUT_CATALOG_TYPE", "SEXTRACTOR")
         if se_output_cat_type in  ["FITS_1.0"]:
            # .FITS SE format
            catalog = FITSCatalog(se_catalog_filepath)
            catalog.open()
            header = catalog.get_header()
            extname = catalog.get_info()["extname"]
            extver  = int(catalog.get_info()["extver"])
            catalog.create_from_numpy(cat_data_matrix, cat_col_names, 
                                      header=header, ext_name=extname, ext_ver=extver)
         else:   
            # --- SExtractor format
            catalog = SExCatalog(se_catalog_filepath)
            catalog.create_from_numpy(cat_data_matrix, 
                                      cat_col_names, cat_col_comments, cat_col_formats)         
            
         catalog.save()
         flagged_object_count = catalog.get_nb_rows()
         catalog.close()          

      return flagged_object_count, flagged_coords


   # -----------------------------------------------------------------------------------------------
   def _sort_by_coordinates(self, se_catalog_filepath, file_type, job, worker):
      """! Sort the catalog by the x and then y coordinates """

      if worker.logging_enabled():
         worker.logger.log_info_p("{0} - /{1}/img-{2:03}-{3:1d} - {4} - Sorting catalog...".format(
                             worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type)) 
         worker.logger.flush()


      final_object_count = 0

      # --- Open SE catalog
      catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)

      # --- Check for na empty catalog
      if catalog.get_nb_rows() > 0:

         # --- Catalog information
         cat_col_names = catalog.get_col_names()
         cat_col_comments = catalog.get_col_comments()
         cat_col_formats = catalog.get_col_formats()
         cat_matrix = catalog.get_data()

         # Sort by x and then y
         x_col_index = catalog.get_col_index("x")
         y_col_index = catalog.get_col_index("y")

         catalog.close()

         sorted_indice = np.lexsort((cat_matrix[:,y_col_index], cat_matrix[:,x_col_index])) 
         sorted_cat_matrix = cat_matrix[sorted_indice].flatten().reshape(cat_matrix.shape)

         se_output_cat_type = worker.config.get_as_string("SE_OUTPUT_CATALOG_TYPE", "SEXTRACTOR")
         if se_output_cat_type in  ["FITS_1.0"]:
            # --- FITS format
            catalog = FITSCatalog(se_catalog_filepath)
            catalog.open()
            header = catalog.get_header()
            extname = catalog.get_info()["extname"]
            extver  = int(catalog.get_info()["extver"])
            catalog.create_from_numpy(sorted_cat_matrix, cat_col_names, 
                                      header=header, ext_name=extname, ext_ver=extver)
         else:
            # --- SExtractor catalog   
            # Save changes by overwriting the catalog
            catalog = SExCatalog(se_catalog_filepath)
            catalog.create_from_numpy(sorted_cat_matrix, 
                                      cat_col_names, cat_col_comments, cat_col_formats)

         catalog.save(write_header=True)
         final_object_count = catalog.get_nb_rows() 

      catalog.close()  

      return final_object_count

   # -----------------------------------------------------------------------------------------------
   def _create_stamp_mosaic(self, image_filepath, 
                                  se_catalog_filepath, deleted_coords, file_type, job, worker):

      if worker.logging_enabled():
         worker.logger.log_info_p(
                           "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Creating stamp mosaic...".format(
                           worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type)) 
         worker.logger.flush()

      # --- Open SE catalog
      catalog = self.helper._open_catalog(se_catalog_filepath, job, worker)

      # --- Catalog information
      x_coords = catalog.get_named_col_data("x").astype(int) 
      y_coords = catalog.get_named_col_data("y").astype(int) 
      obj_coords = zip(x_coords, y_coords)

      catalog.close()

      # --- Output file paths
      result_output_dir = os.path.join(worker.result_output_dir, job.get_branch_tree())
      mosaic_filename_pattern = worker.config.get_as_string("SE_OUPUT_FINAL_CHECK_FILENAME",  
                                                            "SEXTRACTOR")
      file_main, file_ext = os.path.splitext(file_type)
      mosaic_filename = mosaic_filename_pattern.format(file_main)
      mosaic_filename = "{0}-{1:03}-{2:1d}.fits".format(mosaic_filename, job.img_no, job.epoch)
      output_mosaic_filepath = os.path.join(result_output_dir, mosaic_filename)

      # --- Setup image mosaic 
      pixel_size = worker.config.get_as_string("PIXEL_FLOAT_SIZE", 
                                               "PRIMARY_DATASET.IMAGE_PROPERTIES")    
      stamp_size = job.get_stamp_size()  

      image_data = pyfits.getdata(image_filepath).astype(pixel_size)
      (mosaic_height, mosaic_width) = image_data.shape

      mosaic_data = numpy.zeros((mosaic_height, mosaic_width)).astype(pixel_size)

      splitted_obj_mgt = worker.config.get_as_string("METHOD", "CLEANING")

      # --- Populate mosaic
      for (x, y) in obj_coords:

         # --- if DELETE_EXTRA, only populate non-deleted objects
         if splitted_obj_mgt == "DELETE_EXTRA":
            mosaic_data[y:y+stamp_size, x:x+stamp_size] = image_data[y:y+stamp_size, \
                                                                     x:x+stamp_size]
         else:
            if not (y, x) in deleted_coords:
               mosaic_data[y:y+stamp_size, x:x+stamp_size] = image_data[y:y+stamp_size, \
                                                                        x:x+stamp_size]

      # --- Save as FITS file (taking same header as source image)
      self.helper.write_as_fits(mosaic_data, output_mosaic_filepath) 


# -- EOF pse_sex.py
