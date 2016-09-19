"""! 
   Object filtering
"""

# -- Python imports
import os, sys
import numpy

# -- External imports

# --- Module-specific imports
from wgfit_plot import *          # plotter 
from wgfit_helper import *        # helper utility functions



# -------------------------------------------------------------------------------------------------
class ObjectFilter(object):

   """! 
      Galaxy and PSF filtering    
   """

   def __init__(self, request):
      """! 
         GalaxyFilter constructor            
         @param request A ShapeMeasurementRequest object
      """

      self._helper  = WGfitHelper()              # helper utility functions
      self._plotter = WGfitPlotter()             # plotter

      self._request = request                  # shear measurement request
      self._pre_filtered_coords = None         # object coordinates matching pre-filter criteria

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the WGfitHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the WGfitPlotter instance. """
      return self._plotter

   @property
   def request(self):
      """! @return the ShapeMeasurementRequest instance. """
      return self._request

   @property
   def pre_filtered_coords(self):
      """! @return the set of pre_filtered_coordinates. """
      return self._pre_filtered_coords


# -------------------------------------------------------------------------------------------------
class GalaxyFilter(ObjectFilter):
   
   """! 
      Galaxy filtering    
   """

   def __init__(self, request):
      """! GalaxyFilter constructor 
           @param request A ShapeMeasurementRequest object
      """

      ObjectFilter.__init__(self, request)
         
      # --- Object coordinates mattching pre-filtering criteria
      pre_filtered_coords = self._init_pre_filter(request.se_galaxy_info_dico)
      pre_filtered_coords.extend(self._init_pre_filter(request.galaxy_info_dico))
      self._pre_filtered_coords = numpy.unique(pre_filtered_coords)   

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def match(self, coords, custom_data, stamp, job, worker):
      """!
         Predicate indicating whether input data match the filter
         @param coords galaxy coordinate tuple
         @param custom_data specific data corrected after the galaxy fit (such as the residuals)
         @param stamp the postage stamp of the observed galaxy 
         @param job an object of class WGfitJob representing tne current job
         @param worker instance of the worker process
      """

      matched = False
      
      # --- Check pre-filter data, known before fitting
      matched = coords in set(self.pre_filtered_coords) 

      if not matched:

         # --- Check data that depend on fitting results
         if not custom_data is None and len(custom_data) > 0:        
            [convolved_galaxy, residuals] = custom_data[-1]
            if not residuals is None and \
               not numpy.any(numpy.isnan(convolved_galaxy)) and \
               not numpy.any(numpy.isnan(residuals)):

               # --- Check residuals
               norm_residuals = numpy.sum(residuals) / (residuals.shape[0]*residuals.shape[0])    

               [lower, upper] = self.request.galaxy_filter_dico["NORM_RESIDUALS"]
               if lower is None and upper is not None:
                  matched = (norm_residuals > upper)
               elif lower is not None and upper is None: 
                  matched = (norm_residuals < lower)
               elif lower is not None and upper is not None:         
                  matched = (norm_residuals < lower) or (norm_residuals > upper)

               if not matched:
                  max_residuals  = numpy.max(residuals) / (residuals.shape[0]*residuals.shape[0])   
                  [lower, upper] = self.request.galaxy_filter_dico["MAX_RESIDUALS"]
                  if lower is None and upper is not None:
                     matched = (max_residuals > upper)
                  elif lower is not None and upper is None: 
                     matched = (max_residuals < lower)
                  elif lower is not None and upper is not None:         
                     matched = (max_residuals < lower) or (max_residuals > upper)

            else:
               if worker.logging_enabled():
                  worker.logger.log_warning(
                                   "{0} - /{1}/img {2:03}-{3:1d} - NAN value in filter".format(
                                        worker.name, job.get_branch_tree(), job.img_no, job.epoch))
               matched = True    # eliminate nans

      return matched

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _init_pre_filter(self, data_dico):
      """!
         Initialize filtering of data that are available before fitting. 
         The data essentially comes from SEXtractor catalogs. 
         @return galaxy coordinates that match the filter
      """         

      filtered_coords = []

      filtered = [False]

      # --- Go through filtering criteria and collect available data that may be used for filtering
      for var in self.request.galaxy_filter_dico.keys():

         # --- Access the data
         if var in data_dico:

            x_coords, y_coords = zip(*data_dico[var].keys()) 
            x_array = numpy.asarray(x_coords)  
            y_array = numpy.asarray(y_coords)  
            values = numpy.asarray(data_dico[var].values())

            nan_values = numpy.any(numpy.isnan(values))
            if not numpy.any(nan_values):

               # --- Apply filter
               [lower, upper] = self.request.galaxy_filter_dico[var]

               if lower is None and upper is not None:
                  filtered = numpy.logical_or(filtered, values > upper)
               elif lower is not None and upper is None:         
                  filtered = numpy.logical_or(filtered, values < lower)
               elif lower is not None and upper is not None:        
                  filtered = numpy.logical_or(filtered, numpy.logical_or(
                                                                values < lower, values > upper))

               filtered_coords.extend(zip(x_array[filtered], y_array[filtered]))

            else:
               if worker.logging_enabled():
                  worker.logger.log_warning(
                      "{0} - /{1}/img {2:03}-{3:1d} - NAN value detected in filter, "\
                      " variable: {4}".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, var))

               filtered_coords.extend(zip(x_array[nan_values], y_array[nan_values]))

      return list(set(filtered_coords))


# -------------------------------------------------------------------------------------------------
class PSFFilter(ObjectFilter):

   """! 
      PSF (star) filtering    
   """

   def __init__(self, request):
      """! 
         PSFFilter constructor          
         @param request A ShapeMeasurementRequest object
      """

      ObjectFilter.__init__(self, request)

      # --- Object coordinates mattching pre-filtering criteria
      pre_filtered_coords = self._init_pre_filter(request.se_psf_info_dico)
      pre_filtered_coords.extend(self._init_pre_filter(request.psf_info_dico))
      self._pre_filtered_coords = numpy.unique(pre_filtered_coords)   

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def match(self, coords, custom_data, stamp, job, worker):
      """!
         Predicate indicating whether input data match the filter
         @param coords galaxy coordinate tuple
         @param custom_data specific data corrected after the galaxy fit (such as the residuals)
         @param stamp the postage stamp of the observed galaxy 
         @param job an object of class WGfitJob representing tne current job
         @param worker instance of the worker process
      """

      # --- Check pre-filter data, known before fitting
      return coords in set(self.pre_filtered_coords) 

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _init_pre_filter(self, data_dico):
      """!
         Initialize filtering of data that are available before fitting. 
         This data essentially comes from SEXtractor catalogs. 
         @return galaxy coordinates that match the filte
      """      

      filtered_coords = []

      filtered = [False]

      # --- Go through filtering criteria and collect available data that may be use for filtering
      for var in self.request.psf_filter_dico.keys():

         # --- Access the data
         if var in data_dico:

            x_coords, y_coords = zip(*data_dico[var].keys()) 
            x_array = numpy.asarray(x_coords)  
            y_array = numpy.asarray(y_coords)  
            values = numpy.asarray(data_dico[var].values())

            nan_values = numpy.any(numpy.isnan(values))
            if not numpy.any(nan_values):

               # --- Apply filter
               [lower, upper] = self.request.psf_filter_dico[var]

               if lower is None and upper is not None:
                  filtered = numpy.logical_or(filtered, values > upper)
               elif lower is not None and upper is None:         
                  filtered = numpy.logical_or(filtered, values < lower)
               elif lower is not None and upper is not None:        
                  filtered = numpy.logical_or(filtered, numpy.logical_or(values < lower, values > upper))

               filtered_coords.extend(zip(x_array[filtered], y_array[filtered]))

            else:
               if worker.logging_enabled():
                  worker.logger.log_warning(
                      "{0} - /{1}/img {2:03}-{3:1d} - NAN value detected in filter, "\
                      " variable: {4}".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, var))

               filtered_coords.extend(zip(x_array[nan_values], y_array[nan_values]))

      return list(set(filtered_coords))


# -- EOF wgfit_filter.py
