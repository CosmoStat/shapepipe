"""! 
   spredict - Interpolation using Radial Basis Functionc - RBF.py
"""

# --- Python imports
import numpy
import scipy.interpolate

# -- External imports

# -- Module-specific imports
from sp_interp import *    # definition of the BaseInterpolator class
from sp_spatial import *   # distance calculations, nearest neighbors, etc.

# --------------------------------------------------------------------------------------------------
class Method(BaseMethod):

   """! 
      Inerpolation using RBF (Radial basis Functions)
      @note extends parent BaseInterpolator class
   """

   def __init__(self, method_config_dir, interp_config_dir):
      """! 
         Construct an Interpolator object 
         @param method_config_dir configuration directory of the method
         @param interp_config_dir directory of extra configuration files for the method
      """      
      BaseMethod.__init__(self, name="RBF", 
                                method_config_dir=method_config_dir,
                                interp_config_dir=interp_config_dir)
      
      self.is_local = self.method_config.get_as_boolean("IS_LOCAL", 
                                                        "{0}.METHOD".format(self.name.upper()))


   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

   # ----------------------------------------------------------------------------------------------- 
   def predict(self, model, param_name, input_coords, target_coords, input_values, job, worker):
      """!
         Spatially interpolate input values @c input_values of input coordinates 
         @c input_coords to a target grid @ target_coords, using Radial Basis Functions

         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param input_coords input coordinate tuples (x_input_coords, y_input_coords)
         @param target_coords target coordinate tuples (x_target_coords, y_target_coords)
         @param model parameter values at positions @c input_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
      """

      result_dico = {}

      # --- Method-specific default configuration
      config_dico = self._setup_method_config(model, job, worker)

      # --- Read extra configuration information (from interpolate.cfg) and override previous 
      #     definitions if found     
      self.adjust_method_config(config_dico, model, param_name, job, worker)

      result_dico = self._interpolate(model, param_name, config_dico, 
                                      input_coords, target_coords, input_values, job, worker) 

      return result_dico


   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   def _interpolate(self, model, param_name, config_dico, 
                    input_coords, target_coords, input_values, job, worker):


      result_dico = {}

      # --- Input and Target coordinates
      x_input_coords, y_input_coords = zip(*input_coords)
      x_input_coords, y_input_coords = numpy.asarray(x_input_coords), numpy.asarray(y_input_coords)

      x_target_coords, y_target_coords = zip(*target_coords)
      x_target_coords, y_target_coords = numpy.asarray(x_target_coords), \
                                         numpy.asarray(y_target_coords)

      # --- Input values
      input_values = numpy.asarray(input_values)

      # --- Since RBF is a "local" methods, create a SpatialCalculator to efficiently compute
      #     distances find neighbors, etc. on the input data grid
      spatial_calc = SpatialCalculator(input_coords)

      # --- Distance & neighbor search configuration
      min_neighbors, max_neighbors = config_dico["MIN_OBJECTS"], config_dico["MAX_OBJECTS"]
      max_dist = config_dico["MAX_DISTANCE"]
      dist_tol, dist_incr = config_dico["DISTANCE_TOL"], config_dico["DIST_INCREMENT"] 

      # --- Option to force re-evaluation of target points that are located at the same coordinates
      #     as the input points
      force_reevaluation = self.method_config.get_as_boolean("REEVALUATE_IF_EXACT", 
                                                             "{0}.METHOD".format(self.name.upper()))
      # --- May collect extra data for diagnostic purpose
      collect_extra_data = self.method_config.get_as_boolean("COLLECT_EXTRA_DATA", 
                                                             "{0}.METHOD".format(self.name.upper()))

      # --- Compute the RBFs based on the neighorhood settings specified in the configuration
      if worker.logging_enabled():
         worker.logger.log_info(
           "{0} - {1} - Computing RBFs for parameter: {2} ({3}) "\
           "[kernel={4} smooth={5}]...".format(worker.name, self.name, param_name, model.name, 
                                            config_dico["KERNEL"], config_dico["SMOOTH"]))

      target_values = []   # values to predict

      ###for (x,y) in input_coords:      # test
      for (x,y) in target_coords:

         # --- Find a sufficient number of neighboring points
         [obj_list, obj_indice] = spatial_calc.locate_neighbors_within_range(
                                                       (x,y), min_neighbors, max_neighbors, 
                                                       max_dist, dist_tol, dist_incr, dist_type=2)
         
         if len(obj_indice) == 0:

            # --- No neighbor found           
            if worker.logging_enabled():
               err_msg = "{0} - {1} - {2} - No neighbour found at coord: ({3:.2f},{4:.2f}), "\
                         "aborting interpolation process...".format(
                                      worker.name, self.name, param_name, x, y) 
               worker.logger.log_error_p(err_msg)

               raise SpatialInterpolator.InterpolationError(err_msg)

         # --- Too few neighbors found   
         if len(obj_indice) < min_neighbors:
            
            if worker.logging_enabled():
               worker.logger.log_warning_p(
                 "{0} - {1} - {2} - The number of neighbors ({3} < {4}) is too low "\
                 "for accurate RBF evaluation".format(
                      worker.name, self.name, param_name, len(obj_indice), min_neighbors))

         # --- More neighbors found than requested
         if len(obj_indice) > max_neighbors:

            if worker.logging_enabled():
               worker.logger.log_info(
                 "{0} - {1} - {2} - Found ({3} > {4}) neighbors within distance d={5}. "\
                 "Only retaining the first {6} nearest neighbors...".format(
                                        worker.name, self.name, param_name, len(obj_indice), 
                                        max_neighbors, max_dist, max_neighbors))

            obj_list   = obj_list[:max_neighbors]         
            obj_indice = obj_indice[:max_neighbors]          

         # --- Data values of neighboring points
         value_vector = numpy.take(input_values, obj_indice)

         # --- Coordinates of neighboring points
         neighbor_x_coords = numpy.take(x_input_coords, obj_indice)
         neighbor_y_coords = numpy.take(y_input_coords, obj_indice)

         # --- Force reevaluation at location even if the target point was found in the range 
         #     (then we do not return the same point)
         if (x,y) in obj_list and force_reevaluation:
            value_vector = numpy.delete(value_vector, obj_list.index((x, y))) 
            neighbor_x_coords = numpy.delete(neighbor_x_coords, obj_list.index((x,y))) 
            neighbor_y_coords = numpy.delete(neighbor_y_coords, obj_list.index((x,y))) 
            obj_list.remove((x,y))

         # --- Compute the RBFs around (x,y) given the neighborhood found above
         if config_dico["EPSILON"] is None:
            rbfi = scipy.interpolate.Rbf(neighbor_x_coords, neighbor_y_coords, value_vector, 
                                         function=config_dico["KERNEL"], 
                                         smooth=config_dico["SMOOTH"])
         else:
            rbfi = scipy.interpolate.Rbf(neighbor_x_coords, neighbor_y_coords, value_vector, 
                                         function=config_dico["KERNEL"], 
                                         epsilon=config_dico["EPSILON"],
                                         smooth=config_dico["SMOOTH"])

         # --- Now, can evaluate the RBFs at the target locations
         target_values.append(rbfi(x, y))
   

      # --- End for

      # --- Store the results
      result_dico["results"] = {} 
      result_dico["results"]["x_input_coords"]  = x_input_coords
      result_dico["results"]["y_input_coords"]  = y_input_coords
      result_dico["results"]["x_target_coords"] = x_target_coords
      result_dico["results"]["y_target_coords"] = y_target_coords
      result_dico["results"]["input_values"]  = input_values
      result_dico["results"]["target_values"] = target_values

      # --- Extra data
      # TODO: for now, don't record anything
      if collect_extra_data:
         result_dico["extra_data"] = {}

      # --- Summary stats
      # TODO: for now, don't record anything
      result_dico["stats"] = {}
      #result_dico["stats"]["RSS"] = bispl.get_residual()

      return result_dico


   # -----------------------------------------------------------------------------------------------
   def _setup_method_config(self, model, job, worker):
      """! Read interpolation options from method configuration """   

      section_name = "{0}.{1}".format(self.name.upper(), "OPTIONS")

      config_dico = {}

      # --- Distance calculation, neighbor search
      for var in ["MIN_OBJECTS", "MAX_OBJECTS"]:
         config_dico[var] = self.method_config.get_as_int(var, section_name)

      for var in ["DIST_INCREMENT", "MAX_DISTANCE", "DISTANCE_TOL"]:
         config_dico[var] = self.method_config.get_as_float(var, section_name)
   
      config_dico["KERNEL"] = self.method_config.get_as_string("KERNEL", section_name)

      for var in ["EPSILON", "SMOOTH"]:
         config_dico[var] = self.method_config.get_as_string(var, section_name)
         if eval(config_dico[var]) is not None:
             config_dico[var] = self.method_config.get_as_float(var, section_name)  

      return config_dico


# --- EOF RBF.py 
