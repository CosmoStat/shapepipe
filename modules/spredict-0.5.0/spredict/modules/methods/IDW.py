"""! 
   spredict - Interpolation using Radial Basis Functionc - IDW.py
"""

# --- Python imports
import numpy

# -- External imports

# -- Module-specific imports
from sp_interp import *    # definition of the BaseInterpolator class
from sp_spatial import *   # distance calculations, nearest neighbors, etc.

# --------------------------------------------------------------------------------------------------
class Method(BaseMethod):

   """! 
      Interpolation using IDW (Inverse Distance Weighting)
      @note extends parent BaseInterpolator class
   """

   def __init__(self, method_config_dir, interp_config_dir):
      """! 
         Construct an Interpolator object 
         @param method_config_dir configuration directory of the method
         @param interp_config_dir directory of extra configuration files for the method
      """      
      BaseMethod.__init__(self, name="IDW", 
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
         @c input_coords to a target grid @ target_coords, using Inverse Distance Weighting

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

      # --- Since IDW is a "local" methods, create a SpatialCalculator to efficiently compute
      #     distances find neighbors, etc. on the input data grid
      spatial_calc = SpatialCalculator(input_coords)

      # --- Option to force re-evaluation of target points that are located at the same coordinates
      #     as the input points
      force_reevaluation = self.method_config.get_as_boolean("REEVALUATE_IF_EXACT", 
                                                             "{0}.METHOD".format(self.name.upper()))

      collect_extra_data = self.method_config.get_as_boolean("COLLECT_EXTRA_DATA", 
                                                             "{0}.METHOD".format(self.name.upper()))

      # --- Distance & neighbor search configuration
      min_neighbors, max_neighbors = config_dico["MIN_OBJECTS"], config_dico["MAX_OBJECTS"]
      max_dist = config_dico["MAX_DISTANCE"]
      dist_tol, dist_incr = config_dico["DISTANCE_TOL"], config_dico["DIST_INCREMENT"] 

      # --- Distance exponent
      dist_exponent = config_dico["DISTANCE_EXPONENT"] 

      # --- Compute the IDWs based on the neighorhood settings specified in the configuration

      if worker.logging_enabled():
         worker.logger.log_info(
           "{0} - {1} - Conputing IDWs for parameter: {2} ({3}) "\
           "[exponent={4}]...".format(worker.name, self.name, param_name, model.name, 
                                      config_dico["DISTANCE_EXPONENT"]))

      target_values = []   # values predicted at a given location
      weight_list = []     # list of weigjts used for interpolation at a given location

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
                 "for accurate IDW evaluation".format(
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

         # --- Build distance vector from point <location>
         dist_vector = spatial_calc.build_distance_vector_from_location((x,y), obj_list)

         # --- Data values of neighboring points
         value_vector = numpy.take(input_values, obj_indice)

         # --- Coordinates of neighboring points
         neighbor_x_coords = numpy.take(x_input_coords, obj_indice)
         neighbor_y_coords = numpy.take(y_input_coords, obj_indice)

         # --- Compute the IDWs weights around (x,y) given the neighborhood found above
         weights = dist_vector**(-dist_exponent)
         total_weights = numpy.sum(weights)
         weights /= total_weights

         #print 'weights:', weights
         #print 'total_weights:', total_weights, weights.sum()

         # --- Forcing or not reevaluation at location even if the same target point was found in 
         #     the input data
         if (x,y) in obj_list:

            print "IDENTICAL LOCATION", (x,y)    

            if force_reevaluation:

               # --- Compute the IDWs weights around (x,y) given the neighborhood found above
               value_vector = numpy.delete(value_vector, obj_list.index((x, y))) 
               neighbor_x_coords = numpy.delete(neighbor_x_coords, obj_list.index((x,y))) 
               neighbor_y_coords = numpy.delete(neighbor_y_coords, obj_list.index((x,y))) 
               obj_list.remove((x,y))

               # --- Estimated value based on neoghbors values and distances
               target_value = numpy.sum(value_vector * weights) 

            else:

               # --- Estimated value: no re-evaluation, just return the input value at that point
               target_value = value_vector[obj_list.index((x,y))]

         else: 

            # --- Estimated value based on neoghbors values and distances
            target_value = numpy.sum(value_vector * weights) 

         # --- Store results
         target_values.append(target_value)

         # --- Keep tracks of some extra data
         if collect_extra_data:
            weight_list.append(list(weights))    # list of weights

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
         result_dico["extra_data"]

#      # --- Summary stats
#      # TODO: for now, don't record anything
#      result_dico["stats"] = {}


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

      # --- Distance exponent
      config_dico["DISTANCE_EXPONENT"] = self.method_config.get_as_int("DISTANCE_EXPONENT", 
                                                                       section_name)

      return config_dico


# --- EOF IDW.py 
