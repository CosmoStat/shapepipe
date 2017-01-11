"""! 
   spredict - Moffat model - moffat.py
"""

# --- Python import
import numpy

# --- External imports

# --- Module-specific imports
from sp_interp import *       # definition of the BaseModel class
from sp_plot import *         # plotter


# --------------------------------------------------------------------------------------------------
class Model(BaseModel):

   """! 
      Moffat model
      @note extends parent BaseModel class
   """

   def __init__(self, config_dir):
      """! 
         Construct a @c gbdsersic FittingModel object 
         @param config_dir configuration directory of the method
      """  

      BaseModel.__init__(self, name="moffat", config_dir=config_dir)

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def adjust_results(self, data_dico, master):
         
      """! 
         Optionally post-process interpolation results
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
      """
      pass

#      if self.config.get_as_boolean("AUTO_DERIVE_BETA", "MOFFAT.INTERPOLATE"):

#         # --- Apply the 'AUTO_DERIVE_BETA' computation option is selected

#         print data_dico.keys()

#         file_types = data_dico.keys()
#         for file_type in file_types:

#            param_data_dico = data_dico[file_type]

#            # --- Estimate beta values from FWHM values using polynomial regression
#            if "FWHM" in param_data_dico and "beta" in param_data_dico:

#               print param_data_dico["FWHM"][(0,0)].keys()

##               param_data_dico["beta"] = self._get_betas_from_FWHMs(param_data_dico, master)

#      else:
#         # --- Just return the results unaltered
#         return data_dico

#   # -----------------------------------------------------------------------------------------------
#   def self._get_betas_from_FWHMs(self, data_dico, master):

#      """! 
#        Fit a polynomial to beta and FWHM values, in order to derive beta values from FWHM values 
#      """

#      poly_deg = self.config.get_as_int("AUTO_DERIVE_BETA_POLY_DEG", "MOFFAT.INTERPOLATE")    
#      poly_coeffs = numpy.polyfit(data_dico["FWHM"], data_dico["beta"], poly_deg)
#      fitted_betas = numpy.polyval(poly_coeffs, data_dico["FWHM"])

##      if self.config.get_as_boolean("AUTO_DERIVE_BETA_PLOTS", "MOFFAT.INTERPOLATE"):
##         self._plot_scatter_betas_FWHMs(data_dico, fitted_betas, job, master)
#     
#      return fitted_betas




   # -----------------------------------------------------------------------------------------------
   def create_plots(self, method, result_dico, file_type, job, master):

      """! 
         Create plots that are specific to this model 
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
      """

      file_main, file_ext = os.path.splitext(file_type)

      if master.logging_enabled():
         master.logger.log_info_p(
          "{0} - /{1}/{2}-{3:03d}.{4} - {5} ({6}) - Generating input plots... ".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              self.name, job.method_name))

      # --- Create histograms for model parameters
      self._create_model_parameter_histograms(method, result_dico, file_type, 
                                              job, master, data_type="input", pcolor="#2e8b57")

      if master.logging_enabled():
         master.logger.log_info_p(
          "{0} - /{1}/{2}-{3:03d}.{4} - {5} ({6}) - Generating target plots... ".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              self.name , job.method_name))

      # --- Create histograms for model parameters
      self._create_model_parameter_histograms(method, result_dico, file_type, 
                                              job, master, data_type="target", pcolor="#6495ed")


   # -----------------------------------------------------------------------------------------------
   def create_maps(self, method, result_dico, file_type, job, master):
      """! 
         Create maps that are specific to this model 
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
      """

      file_main, file_ext = os.path.splitext(file_type)

      if master.logging_enabled():
         master.logger.log_info_p(
          "{0} - /{1}/{2}-{3:03d}.{4} - {5} ({6}) - Generating (e1, e2) input maps... ".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              self.name , job.method_name))

      # --- Create ellipticity vector maps   
      self._create_e1_e2_vector_maps(method, result_dico, file_type, 
                                     job, master, data_type="input", scale=5, width=0.001) 

      if master.logging_enabled():
         master.logger.log_info_p(
          "{0} - /{1}/{2}-{3:03d}.{4} - {5} ({6}) - Generating (e1, e2) target maps... ".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              self.name , job.method_name))

      # --- Create ellipticity vector maps   
      self._create_e1_e2_vector_maps(method, result_dico, file_type,  
                                     job, master, data_type="target", scale=2.5, width=0.001) 

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _flatten_tiled_param_dico(self, tiled_dico, param_name):

      data_list = []

      tile_tuples = sorted(tiled_dico.keys()) # sorted by tile tuple (0,0), (0,1)...(1,0), (1.1)...
      for (x_tile_no, y_tile_no) in tile_tuples:
         data_list.extend(
                       tiled_dico[(x_tile_no, y_tile_no)][param_name]["results"]["target_values"])

      return data_list

   # -----------------------------------------------------------------------------------------------
   def _create_e1_e2_vector_maps(self, 
                                 method, result_dico, file_type, job, master, 
                                 data_type, scale, width):
      """! 
         Create vector "quiver" maps that are specific to this model
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
      """

      plotter = SpredictPlotter()

      # --- Which maps to create
      map_section_name = "{0}.MAPS.{1}".format(self.name.upper(), data_type.upper())
        
      create_e1_e2_field_vector_maps = self.config.get_as_boolean(
                                                       "CREATE_FIELD_ELLIPTICITY_VECTOR_MAPS", 
                                                       map_section_name)
      create_e1_e2_tile_vector_maps  = self.config.get_as_boolean(
                                                       "CREATE_TILE_ELLIPTICITY_VECTOR_MAPS", 
                                                       map_section_name)

      # --- Where maps will be stored
      branch_tree = job.get_branch_tree()
      map_output_dir = os.path.join(master.map_output_dir, branch_tree)

      # --- Coordinates and valuesfor the global map
      glob_x_coords, glob_y_coords = [], []
      glob_x_values, glob_y_values = [], []

      # --- Create a map for each tile...
      tile_tuples = sorted(result_dico.keys())
      
      for (x_tile_no, y_tile_no) in tile_tuples:

         if "e1" in result_dico[(x_tile_no, y_tile_no)] and \
            "e2" in result_dico[(x_tile_no, y_tile_no)]:

            # --- Access the data to map
            x_data = result_dico[(x_tile_no, y_tile_no)]["e1"]["results"]
            y_data = result_dico[(x_tile_no, y_tile_no)]["e2"]["results"]

            x_coords = x_data["x_"+data_type+"_coords"]
            y_coords = y_data["y_"+data_type+"_coords"]
            x_values = x_data[data_type+"_values"]
            y_values = y_data[data_type+"_values"]

            glob_x_coords.extend(x_coords)
            glob_y_coords.extend(y_coords)
            glob_x_values.extend(x_values)
            glob_y_values.extend(y_values)

            # --- Create tile maps if requested 
            if create_e1_e2_tile_vector_maps:

               # --- Ellipticity components (e1, e2) for the current tile

               map_title = "{0}-{1} ({2}) - {3} tile map - Image {4:03d}-{5} - Tile ({6},{7})".format(
                  "e1", "e2", method.name, data_type, job.img_no, job.epoch, x_tile_no, y_tile_no)

               map_filename = "{0}_{1}-{2}_tile_map_{3}_{4:03d}_{5:1d}_tile_{6}_{7}.png".format(
                  data_type, "e1", "e2", method.name, job.img_no, job.epoch, x_tile_no, y_tile_no)

               plotter.create_image_quiver_plots(x_coords, y_coords, x_values, y_values,
                                                 "e1", "e2", file_type, 
                                                 map_title, map_output_dir, map_filename,
                                                 job, master,
                                                 scale=scale, width=width, pfilter=None)

      if len(glob_x_values) > 0 and len(glob_y_values) > 0: 

         # --- Create global field maps if requested 
         if create_e1_e2_field_vector_maps:

            # --- Create a global map from true field positions...

            # --- Aggregation of the field coordinates from all tiles
            glob_x_coords = numpy.asarray(glob_x_coords)
            glob_y_coords = numpy.asarray(glob_y_coords)
            glob_x_values = numpy.asarray(glob_x_values)  
            glob_y_values = numpy.asarray(glob_y_values)  

            # --- Sorted coordinates and associated z-values
            sorting_indice = numpy.lexsort((glob_y_coords, glob_x_coords)) # sorted
            sorted_glob_x_coords = glob_x_coords[sorting_indice]
            sorted_glob_y_coords = glob_y_coords[sorting_indice]
            sorted_glob_x_values = glob_x_values[sorting_indice]
            sorted_glob_y_values = glob_y_values[sorting_indice]

            # --- Setup map titles and filenames
            map_title = "{0}-{1} ({2}) - {3} field vector map - Image {4:03d}-{5}".format(
               "e1", "e2",  method.name, data_type, job.img_no, job.epoch)

            map_filename = "{0}_{1}-{2}_field_vector_map_{3}_{4:03d}_{5}.png".format(
               data_type, "e1", "e2", method.name, job.img_no, job.epoch)

            # --- Create the global field map 
            plotter.create_image_quiver_plots(
                                           sorted_glob_x_coords, sorted_glob_y_coords, 
                                           sorted_glob_x_values, sorted_glob_y_values,
                                           "e1", "e2", file_type, 
                                           map_title, map_output_dir, map_filename,
                                           job, master, scale=scale, width=width, pfilter=None)

   # -----------------------------------------------------------------------------------------------
   def _create_model_parameter_histograms(self, method, result_dico, file_type, 
                                               job, master, data_type, pcolor='b'):
      """! Create model parameters histograms """

      # --- Create histogram plots?
      plot_section_name = "{0}.PLOTS.{1}".format(self.name.upper(), data_type.upper())
      if self.config.get_as_boolean("CREATE_TILE_PARAMETER_HISTOGRAMS", plot_section_name):

         plotter = SpredictPlotter()

         # --- Create a plot for each tile...
         branch_tree = job.get_branch_tree()
         plot_output_dir = os.path.join(master.plot_output_dir, branch_tree)

         tile_tuples = sorted(result_dico.keys())
         for (x_tile_no, y_tile_no) in tile_tuples:
            
            for param_name in result_dico[(x_tile_no, y_tile_no)].keys():

               # --- Access the data to plot
               param_data = result_dico[(x_tile_no, y_tile_no)][param_name]["results"]
               param_values = param_data[data_type+"_values"]                

               # --- Mean_value of main variable
               var_mean = numpy.mean(param_values)                
               var_std = numpy.std(param_values) 

               plot_title = "Tile {0}-{1} ({2}) - {3} Parameter {4} - Image {5:03d}-{6} - "\
                            "[{7:.3f}, {8:.3f}]".format(
                  x_tile_no, y_tile_no, method.name, data_type, param_name, job.img_no, job.epoch, 
                  var_mean, var_std)

               plot_filename = "{0}_tile_{1}-{2}_hist_{3}_{4}_{5:03d}_{6}_{7:.3f}_{8:.3f}.png".format(
                  data_type, x_tile_no, y_tile_no, method.name, param_name, job.img_no, job.epoch, 
                  var_mean, var_std)

               plotter.plot_histogram(param_data, [data_type+"_values"],
                                   plot_filename, plot_output_dir, plot_title, pfilter=None,
                                   pnb_bins=25.0, pmin_xlim=None, pmax_xlim=None,
                                   pcolor=pcolor, is_normed=False)


