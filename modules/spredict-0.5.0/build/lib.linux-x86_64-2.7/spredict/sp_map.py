"""! 
   sp_map.py - Mapping functions.
"""

# -- Python imports
import os, sys
import math
import numpy
import matplotlib

matplotlib.use("AGG")
import matplotlib.pyplot as plt
import pylab
import mpl_toolkits.mplot3d.axes3d as pylab3
from matplotlib.font_manager import FontProperties

# -- External imports


# --- Module-specific imports
from sp_plot import *         # plotter 

# -------------------------------------------------------------------------------------------------
class SpredictMapper(object):
   
   """! 
      Class with a number of convenient mapping methods.   
   """

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def create_maps(self, model, method, result_dico, file_type, job, master):
      """!
         Create two-dimensional maps
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
      """

      file_main, file_ext = os.path.splitext(file_type)

      # --- Model-independent maps (input - output)

      if master.logging_enabled():
         master.logger.log_info_p(
          "{0} - /{1}/{2}-{3:03d}.{4} - Generating Input maps for parameter(s) "\
          "{5}\t({6}) ...".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              job.param_names, job.method_name))

      self.create_scatter_maps(model, method, result_dico, job, master, data_type="input")        

      # --- Scatter Target maps
      if master.logging_enabled():
         master.logger.log_info_p(
           "{0} - /{1}/{2}-{3:03d}.{4} - Generating Target maps for parameter(s) "\
           "{5}\t({6}) ...".format(
              master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
              job.param_names, job.method_name))

      self.create_scatter_maps(model, method, result_dico, job, master, data_type="target")        
 
      # --- Model-specific maps (input - output)
      model.create_maps(method, result_dico, file_type, job, master) 

   


   # -----------------------------------------------------------------------------------------------
   def create_scatter_maps(self, model, method, result_dico, job, master, data_type):
      """!
         Create scatter maps
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
         @param data_type type of data to map: "input" or "target"
      """

      # --- Which maps to create
      map_section_name = "MAPS.{0}".format(data_type.upper())
        
      create_field_maps = master.config.get_as_boolean("CREATE_FIELD_SCATTER_MAPS", 
                                                       map_section_name)
      create_tile_maps  = master.config.get_as_boolean("CREATE_TILE_SCATTER_MAPS", 
                                                       map_section_name)

      # --- Where maps will be stored
      branch_tree = job.get_branch_tree()
      map_output_dir = os.path.join(master.map_output_dir, branch_tree)

      # --- Coordinates and valuesfor the global map
      glob_x_coords, glob_y_coords = [], []
      glob_z_values = []

      # --- Create a map for each tile...
      tile_tuples = sorted(result_dico.keys())
      all_param_names = []      

      for (x_tile_no, y_tile_no) in tile_tuples:

         param_names = sorted(result_dico[(x_tile_no, y_tile_no)].keys()) 
         for param_name in param_names:

            # --- Access the data to map
            all_param_names.append(param_name)
            param_data = result_dico[(x_tile_no, y_tile_no)][param_name]["results"]
            x_coords = param_data["x_"+data_type+"_coords"]
            y_coords = param_data["y_"+data_type+"_coords"]
            z_values = param_data[data_type+"_values"]

            glob_x_coords.extend(x_coords)
            glob_y_coords.extend(y_coords)
            glob_z_values.extend(z_values)

            # --- Create tile maps if requested 
            if create_tile_maps:
               # --- Setup map titles and filenames
               map_title = "{0} ({1}) - {2} tile map - Image {3:03d}-{4} - Tile ({5},{6})".format(
                  param_name, method.name, data_type, job.img_no, job.epoch, x_tile_no, y_tile_no)

               map_filename = "{0}_{1}_tile_map_{2}_{3:03d}_{4}_tile_{5}_{6}.png".format(
                  data_type, param_name, method.name, job.img_no, job.epoch, x_tile_no, y_tile_no)

               # --- Create a field map for the current tile
               self.create_scatter_map(x_coords, y_coords, z_values, map_title=map_title, 
                                   has_grid=False, output_dir=map_output_dir, output_file=map_filename, 
                                   pfigure=None, show_levels=False, pcmap=None, pcbar=None, 
                                   pmeshcolor=True, show=False)

      # --- Create global field maps if requested 
      if create_field_maps:

         for param_name in numpy.unique(all_param_names):

            # --- Create a global map from true field positions...

            # --- Aggregation of the field coordinates from all tiles
            glob_x_coords = numpy.asarray(glob_x_coords)
            glob_y_coords = numpy.asarray(glob_y_coords)
            glob_z_values = numpy.asarray(glob_z_values)  

            # --- Sorted coordinates and associated z-values
            sorting_indice = numpy.lexsort((glob_y_coords, glob_x_coords)) # sorted
            sorted_glob_x_coords = glob_x_coords[sorting_indice]
            sorted_glob_y_coords = glob_y_coords[sorting_indice]
            sorted_glob_z_values = glob_z_values[sorting_indice]

            # --- Setup map titles and filenames
            map_title = "{0} ({1}) - {2} field map - Image {3:03d}-{4}".format(
               param_name, method.name, data_type, job.img_no, job.epoch)

            map_filename = "{0}_field_map_{1}_{2}_{3:03d}_{4}.png".format(
               data_type, param_name, method.name, job.img_no, job.epoch)

            # --- Create the global field map 
            self.create_scatter_map(
                                sorted_glob_x_coords, sorted_glob_y_coords,
                                sorted_glob_z_values, map_title=map_title, 
                                has_grid=False, output_dir=map_output_dir, output_file=map_filename, 
                                pfigure=None, show_levels=False, pcmap=None, pcbar=None, 
                                pmeshcolor=True, show=False)


   # -----------------------------------------------------------------------------------------------
   def create_scatter_map(self, x_coords, y_coords, z_values, map_title='Scatter Map', 
                                has_grid=False, output_dir='.', output_file='scatter.map.png', 
                                pfigure=None, show_levels=False, pcmap=None, pcbar=None, 
                                pmeshcolor=False, show=False):
      
      """! Create a scatter map (X,Y,Z) """
  
      # --- Coordinates boundaries
      pmin_xlim, pmin_ylim = min(x_coords), min(y_coords)
      pmax_xlim, pmax_ylim = max(x_coords), max(y_coords)

      data_dico = {}
      data_dico['X'], data_dico['Y'] = x_coords, y_coords

      if pmeshcolor:
         xi, yi, zi = self.map_interpolated_data_on_grid(x_coords, y_coords, z_values, 
                                                         x_step=1000, y_step= 1000, 
                                                         interp_method='nn')
      if pfigure is None:
         pfigure = pylab.figure()

      # --- Default color map
      if pcmap is None:
         pcmap = pylab.get_cmap('jet')

      if pmeshcolor:
         # --- Color mesh
         pylab.pcolormesh(xi, yi, zi, cmap=pcmap, edgecolors='face', shading='gouraud')

      facecolors = z_values

      # --- Create the scatter
      plotter = SpredictPlotter()
      plot = plotter.plot_scatter(data_dico, ['X', 'Y'],
                           output_file, output_dir, plot_title=map_title, pfilter=None,
                           pxmin_lim=pmin_xlim, pxmax_lim=pmax_xlim , 
                           pymin_lim=pmin_ylim, pymax_lim=pmax_ylim, 
                           x_label='X (degs)', y_label='Y (degs)', 
                           pcolor=facecolors, pmarker='o', pmarkersize=12, 
                           pmarkerfacecolor=facecolors, pmarkeredgecolor='face', pcmap=pcmap,
                           is_normed=False, has_equalaxes=False, has_grid=False, pfigure=pfigure)

      # --- Default color bar
      if pcbar is None:
         pylab.colorbar() 

      if not pfigure is None:
         # --- Save map and free memory
         pylab.savefig(os.path.join(output_dir, output_file))
         pylab.clf()
         pylab.close()


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Interpolate a subset of data on a grid  
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def map_interpolated_data_on_grid(self, x_coords, y_coords, z_values, x_step, y_step, 
                                           interp_method='linear'):
      """!
         Interpolate a subset of data on a grid  
      """

      # --- Interpolate over the target grid 
      coord_tuples = zip(x_coords, y_coords)
      x_grid_size, y_grid_size = max(x_coords), max(y_coords)

      x_lin_space = numpy.linspace(min(x_coords), x_grid_size, x_step)
      y_lin_space = numpy.linspace(min(y_coords), y_grid_size, y_step)

      xi, yi = numpy.meshgrid(x_lin_space, y_lin_space)
      zi = pylab.griddata(x_coords, y_coords, z_values, xi, yi, interp=interp_method)

      return xi, yi, zi


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~


# -- EOF sp_map.py
