"""! 
   sp_helper.py - various helper functions
""" 

# -- Python imports
import os, sys
import math
import numpy
import scipy
import pyfits

# -- External import
from scatalog import *           # file catalog management             
from mpf.mp_helper import *      # base Helper

# -- Module-specific imports  
from sp_helper import *         # this module helper utility functions


# -------------------------------------------------------------------------------------------------
class SpredictHelper(Helper):

   """! Convenient utility functions """

   # -----------------------------------------------------------------------------------------------
   def get_spredict_version(self):
      """! 
          Get the version number of this spredict code as text 
          @return version number of this spredict code as text
      """
      return "v0.5.0" 


   # -----------------------------------------------------------------------------------------------
   def show_config_summary(self, master):
      
      if master.logging_enabled():

         # --- spredict version
         master.logger.log_info_p("\n*** spredict {0} ***\n".format(self.get_spredict_version()))

         # --- configuration
         spredict_config_path = os.path.join(master.default_config_dir, master.default_config_filename)
#         model_config_path   = master.config.get_as_string("MODEL_CONFIG_PATH", "MULTIFIT")
#         method_config_path  = master.config.get_as_string("METHOD_CONFIG_PATH", "MULTIFIT")
#         fitting_config_path = master.config.get_as_string("FITTING_CONFIG_PATH", "MULTIFIT")

         master.logger.log_info_p("Base configuration: {0}".format(spredict_config_path))
#         master.logger.log_info_p("- Fitting methods : {0}".format(method_config_path))
#         master.logger.log_info_p("- Galaxy models   : {0}".format(model_config_path))
#         master.logger.log_info_p("- Custom fitting  : {0}\n".format(fitting_config_path))

#         # --- Selected minimizer & galaxy model
#         galaxy_model = master.config.get_as_string("GALAXY_MODEL_NAME", "GALAXY_FITTING")
#         fitting_method = master.config.get_as_string("GALAXY_FITTING_METHOD", "GALAXY_FITTING")
#         master.logger.log_info_p("Selected galaxy model: {0}".format(galaxy_model))
#         master.logger.log_info_p("Selected minimizer   : {0}\n".format(fitting_method))

         # --- Directories
         #master.logger.log_info_p("Base input directory   : {0}".format(master.base_input_dir))
         master.logger.log_info_p("Base output directory : {0}".format(master.base_output_dir))
         master.logger.log_info_p("Run output directory  : {0}".format(master.run_output_dir))
         master.logger.log_info_p("log output directory  : {0}".format(master.log_output_dir))
         master.logger.log_info_p("plot output directory : {0}".format(master.plot_output_dir))
         master.logger.log_info_p("stats output directory: {0}".format(master.stat_output_dir))
         master.logger.log_info_p("error output directory: {0}\n".format(master.error_output_dir))

         # --- Datasets   
         src_dataset_dir = master.config.get_as_string("BASE_DIR", "DATASET")
#         gal_se_dataset_dir = master.config.get_as_string("BASE_DIR", "SEXTRACTOR_DATASET.GALAXY")
#         psf_se_dataset_dir = master.config.get_as_string("BASE_DIR", "SEXTRACTOR_DATASET.PSF")
#         psf_dataset_dir  = master.config.get_as_string("BASE_DIR", "PSF_DATASET")

         master.logger.log_info_p("Main source dataset: {0}".format(src_dataset_dir))
#         master.logger.log_info_p("SExtractor galaxy dataset: {0}".format(gal_se_dataset_dir))
#         master.logger.log_info_p("SExtractor psf dataset   : {0}".format(psf_se_dataset_dir))
#         master.logger.log_info_p("PSF dataset              : {0}\n".format(psf_dataset_dir))

         

   # -----------------------------------------------------------------------------------------------
   def get_nb_objects(self, image_filepath, stamp_size):
      """! Find the actual number of objects (i.e. postage stamps) in the image """
      return int(len(pyfits.getdata(image_filepath)) / stamp_size)**2
   

   # -----------------------------------------------------------------------------------------------
   def get_catalog_filepath(self, filepath_dico, cat_extension, job):
      """! 
         Get the full path of the input catalog containing at least the coordinates of the objects
         for interpolation
         @param job SpredictJob object
         @param worker Worker process object
         @param cat_extension catalog extension (".fits" or ".txt")
         @return full path of the catalog containing the galaxy coordinates
      """

      match = None

      # --- looking for catalog filepaths with key (img_no, None) in path dico
      #     (only images have a non-None epoch number, so we get only catalogs)
      if len(filepath_dico) > 0 and \
         (job.img_no, None) in filepath_dico:
         filepaths = filepath_dico.get((job.img_no, None), None)
         if filepaths is not None:
            matches = [e for e in filepaths if cat_extension in e]             
            if len(matches) > 0:
               match = matches[0] 
      return match

   # -----------------------------------------------------------------------------------------------
   def get_nb_tiles(self, catalog_filepath, catalog_file_extension, catalog_is_SE,
                           hdu_no, job, worker):
      """!
         Count the number of tiles in a catalog file, assuming the number of tiles on the x and y
         axes are the same
      """

      # --- Open catalog
      if catalog_file_extension == ".fits":
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
      elif catalog_file_extension == ".txt":
         if catalog_is_SE:
            catalog = SExCatalog(catalog_filepath)
         else:
            catalog = TextCatalog(catalog_filepath)
      else:
         return None

      # --- Count the number of tiles
      catalog.open()
      nb_tiles = len(numpy.unique(catalog.get_named_col_data("x_tile_index")))      
      catalog.close()

      return nb_tiles 

   # -----------------------------------------------------------------------------------------------
   def get_IDs_per_tile_dico(self, catalog_filepath, catalog_file_extension, 
                                   catalog_is_SE, hdu_no, nb_tiles, job, worker):

      # --- Open catalog
      if catalog_file_extension == ".fits":
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
      elif catalog_file_extension == ".txt":
         if catalog_is_SE:
            catalog = SExCatalog(catalog_filepath)
         else:
            catalog (catalog_filepath)
      else:
         return None

      # --- Create dictionary of tile positions indexed by tile no
      catalog.open()

      x_tile_index = catalog.get_named_col_data("x_tile_index") 
      y_tile_index = catalog.get_named_col_data("y_tile_index") 

      # --- Stamp psotions within an image
      IDs = catalog.get_named_col_data("ID")

      catalog.close()

      ids_per_tile_dico   = {}    # positions within tile per tile

      for x_tile_no in xrange(nb_tiles):
         for y_tile_no in xrange(nb_tiles):

            tile_filter = numpy.logical_and(x_tile_index == x_tile_no, y_tile_index == y_tile_no)
            pos_indice = numpy.ravel(numpy.where(tile_filter))
            #y_pos_indice = numpy.ravel(numpy.where(tile_filter))

            #print (x_tile_no, y_tile_no), "len(pos_indice):", len(pos_indice)

            ids_per_tile_dico[(x_tile_no, y_tile_no)] = IDs[pos_indice]

      return ids_per_tile_dico


   # -----------------------------------------------------------------------------------------------
   def get_positions_per_tile_dico(self, catalog_filepath, catalog_file_extension, 
                                         catalog_is_SE, hdu_no, XY_labels, nb_tiles, 
                                         job, worker):

      # --- Open catalog
      if catalog_file_extension == ".fits":
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
      elif catalog_file_extension == ".txt":
         if catalog_is_SE:
            catalog = SExCatalog(catalog_filepath)
         else:
            catalog = TextCatalog(catalog_filepath)
      else:
         return None

      # --- Create dictionary of tile positions indexed by tile no
      catalog.open()

      x_tile_index = catalog.get_named_col_data("x_tile_index") 
      y_tile_index = catalog.get_named_col_data("y_tile_index") 

      if not catalog.col_exists(XY_labels[0]) or not catalog.col_exists(XY_labels[1]):
         worker.helper.print_error("Column {0} or {1} not found in catalog: {2}".format(
                                                    XY_labels[0], XY_labels[1], catalog_filepath))
         return None

      # --- Positions (degs) within a tile
      tile_x_pos_deg = catalog.get_named_col_data(XY_labels[0]) 
      tile_y_pos_deg = catalog.get_named_col_data(XY_labels[1]) 

      catalog.close()

      pos_per_tile_dico   = {}    # positions within tile per tile

      for x_tile_no in xrange(nb_tiles):
         for y_tile_no in xrange(nb_tiles):

            tile_filter = numpy.logical_and(x_tile_index == x_tile_no, y_tile_index == y_tile_no)
            pos_indice = numpy.ravel(numpy.where(tile_filter))
            #y_pos_indice = numpy.ravel(numpy.where(tile_filter))

            #print (x_tile_no, y_tile_no), "len(pos_indice):", len(pos_indice)

            pos_per_tile_dico[(x_tile_no, y_tile_no)] = zip(tile_x_pos_deg[pos_indice],
                                                            tile_y_pos_deg[pos_indice])

      return pos_per_tile_dico


   # -----------------------------------------------------------------------------------------------
   def get_values_per_tile_dico(self, param_names, catalog_filepath, 
                                      catalog_file_extension, catalog_is_SE,
                                      hdu_no, nb_tiles, job, worker):

      # --- Open catalog
      if catalog_file_extension == ".fits":
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
      elif catalog_file_extension == ".txt":
         if catalog_is_SE:
            catalog = SExCatalog(catalog_filepath)
         else:
            catalog = TextCatalog(catalog_filepath)
      else:
         return None

      # --- Create dictrionary of tile positions indexed by tile no
      catalog.open()

      x_tile_index = catalog.get_named_col_data("x_tile_index") 
      y_tile_index = catalog.get_named_col_data("y_tile_index") 

      values_per_tile_dico = {}

      for param_name in param_names:

         values_per_tile_dico[param_name] = {}
         param_values = catalog.get_named_col_data(param_name) 

         for x_tile_no in xrange(nb_tiles):
            for y_tile_no in xrange(nb_tiles):

               tile_filter = numpy.logical_and(x_tile_index == x_tile_no, y_tile_index == y_tile_no)
               pos_indice = numpy.ravel(numpy.where(tile_filter))

               values_per_tile_dico[param_name][(x_tile_no, y_tile_no)] = param_values[pos_indice]

               #print param_name, (x_tile_no, y_tile_no), values_per_tile_dico[param_name]

      catalog.close()

      return values_per_tile_dico

   # -----------------------------------------------------------------------------------------------
   def get_field_positions_per_tile_dico(self, catalog_filepath, catalog_file_extension, 
                                         catalog_is_SE, hdu_no, XY_labels, nb_tiles, 
                                         job, worker):

      # --- Open catalog
      if catalog_file_extension == ".fits":
         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
      elif catalog_file_extension == ".txt":
         if catalog_is_SE:
            catalog = SExCatalog(catalog_filepath)
         else:
            catalog = TextCatalog(catalog_filepath)
      else:
         return None

      # --- Create dictionary of tile positions indexed by tile no
      catalog.open()

      # --- Tile indice
      x_tile_index = catalog.get_named_col_data("x_tile_index") 
      y_tile_index = catalog.get_named_col_data("y_tile_index") 

      if not catalog.col_exists(XY_labels[0]) or not catalog.col_exists(XY_labels[1]):
         worker.helper.print_error("Column {0} or {1} not found in catalog: {2}".format(
                                                    XY_labels[0], XY_labels[1], catalog_filepath))
         return None

      # --- True field positions
      field_x_pos_deg = catalog.get_named_col_data(XY_labels[0]) 
      field_y_pos_deg = catalog.get_named_col_data(XY_labels[1])  

      catalog.close()

      field_pos_dico   = {}    # field positions across tiles per tile

      for x_tile_no in xrange(nb_tiles):
         for y_tile_no in xrange(nb_tiles):

            tile_filter = numpy.logical_and(x_tile_index == x_tile_no, y_tile_index == y_tile_no)
            pos_indice = numpy.ravel(numpy.where(tile_filter))
            #y_pos_indice = numpy.ravel(numpy.where(tile_filter))

            #print (x_tile_no, y_tile_no), "len(pos_indice):", len(pos_indice)

            #print (x_tile_no, y_tile_no), field_x_pos_deg[pos_indice]

#            print "*** field_pos_dico:", (x_tile_no, y_tile_no), field_pos_dico[(x_tile_no, y_tile_no)]

            field_pos_dico[(x_tile_no, y_tile_no)] = zip(field_x_pos_deg[pos_indice],
                                                         field_y_pos_deg[pos_indice])




      return field_pos_dico

#   # -----------------------------------------------------------------------------------------------
#   def get_field_positions_per_tile_dico(self, catalog_filepath, catalog_file_extension, 
#                                         catalog_is_SE, hdu_no, XY_labels, nb_tiles, 
#                                         job, worker):

#      # --- Open catalog
#      if catalog_file_extension == ".fits":
#         catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no, mem_map=True)
#      elif catalog_file_extension == ".txt":
#         if catalog_is_SE:
#            catalog = SExCatalog(catalog_filepath)
#         else:
#            catalog = TextCatalog(catalog_filepath)
#      else:
#         return None

#      # --- Create dictionary of tile positions indexed by tile no
#      catalog.open()

#      # --- Tile indice
#      x_tile_index = catalog.get_col_data("x_tile_index") 
#      y_tile_index = catalog.get_col_data("y_tile_index") 

#      if not catalog.col_exists(XY_labels[0]) or not catalog.col_exists(XY_labels[1]):
#         worker.helper.print_error("Column {0} or {1} not found in catalog: {2}".format(
#                                                    XY_labels[0], XY_labels[1], catalog_filepath))
#         return None

#      # --- Positions (degs) within a tile
#      tile_x_pos_deg = catalog.get_col_data(XY_labels[0]) 
#      tile_y_pos_deg = catalog.get_col_data(XY_labels[1]) 

#      # --- True field positions
#      field_x_pos_deg = catalog.get_col_data("x_field_true_deg") 
#      field_y_pos_deg = catalog.get_col_data("y_field_true_deg") 

#      catalog.close()

#      field_pos_dico = dict(zip(zip(x_tile_index, y_tile_index, tile_x_pos_deg, tile_y_pos_deg), 
#                                zip(field_x_pos_deg, field_y_pos_deg)))

#      print "*** field_pos_dico:", field_pos_dico.keys(), len(field_pos_dico)

#      return field_pos_dico


   # -----------------------------------------------------------------------------------------------
   def save_from_list_dico(self, data_dico, output_directory, output_filename, col_list = [], 
                                 key_index_map=[], key_fmt_map=[], default_fmt="%.18e"):
      """! 
         Save a dictionary to disk as a catalog file. It is assumed a list of values is attached to
         each first-level key in the dictionary (a "list" dicrionary)
         @param data_dico dictionary with the data
         @param output_directory directory where to create the file
         @param output_filename name of the file to create
         @param col_list list of column names. If empty, take all the keys of the dictionary 
         @param key_index_map if not empty, contains a map with the preferred order for some keys
         @param key_fmt_map if not empty, contains the preferred output format of some key values
      """

      output_catalog = None    # output catalog

      try:

         #print data_dico

         # --- Create output file
         output_filepath = os.path.join(output_directory, output_filename)

         # --- Build the list of columns in the required order
         cat_col_list = []
         cat_col_fmt = []
         cat_col_comments = ""

         if len(col_list) == 0:
            col_names = data_dico.keys()
         else:
            col_names = col_list

         # --- Check the column names are indeed in the dictionary
         col_names_to_remove = []
         for col_name in col_names:
            if not col_name in data_dico:
               self.print_warning( "column: {0} not found in the dictionary".format(col_name) )
               col_names_to_remove.append(col_name)
         for col_name in col_names_to_remove:
            self.print_warning( "Removing column: {0}".format(col_name) )
            col_names.remove(col_name)

         for col_name in col_names:
            col_no = col_names.index(col_name)
            if col_name in key_index_map:
               cat_col_list.insert(key_index_map[col_name], col_name)
               if col_name in key_fmt_map:
                  cat_col_fmt.insert(key_index_map[col_name], key_fmt_map[col_name])
               else:
                  cat_col_fmt.insert(key_index_map[col_name], default_fmt)
            else:
               cat_col_list.append(col_name)
               cat_col_fmt.append(default_fmt)

         # --- Insert the columns in the catalog
         col_data_list = []   
         for col_name in cat_col_list:
           cat_col_comments += col_name + " " 
           col_data_list.append(numpy.asarray([]))

         for col_name in cat_col_list:
            col_no = cat_col_list.index(col_name)
            col_data_list[col_no] = numpy.concatenate(  
                                                    (col_data_list[col_no], data_dico[col_name] ) ) 

         data_matrix = numpy.asmatrix(col_data_list).transpose().squeeze()

         numpy.savetxt(output_filepath, data_matrix, fmt=cat_col_fmt, header=cat_col_comments)

      except Exception as detail:
         self.print_error("could not create catalog from dictionary ({0})".format(detail))


# -- EOF sp_helper.py

