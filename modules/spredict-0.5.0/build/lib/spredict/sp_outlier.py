"""! 
   sp_outlier.py - Outlier detection and removal
"""

# -- Python imports
import os, sys
import string
import time
import math
import numpy

# --- Module-specific imports
from sp_helper import *       # helper utility functions
from sp_plot import *         # plotter 


# -------------------------------------------------------------------------------------------------
class DataCleaner(object):
   """! Detect and Remove outliers from a 2D square data grid """

   def __init__(self, master):
      """! DataCleaner constructor """

      self._helper  = SpredictHelper()              # helper utility functions
      self._plotter = SpredictPlotter()             # plotter

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the GfitHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the GfitPlotter instance. """
      return self._plotter

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # ----------------------------------------------------------------------------------------------- 
   def remove_outliers(self, model, param_name, input_coords, input_values, job, worker):
      """!
         Remove outliers from input values @c input_values of input coordinates 
         @c input_coords

         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param input_coords input coordinate tuples (x_input_coords, y_input_coords)
         @param model parameter values at positions @c input_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
      """

      new_input_coords = input_coords
      new_input_values = input_values

      # --- Input and Target coordinates
      x_coords, y_coords = zip(*input_coords)
      x_coords, y_coords = numpy.asarray(x_coords), numpy.asarray(y_coords)

      # --- Input values
      values = numpy.asarray(input_values)

      #print "### len coords:", len(x_coords), len(y_coords)
      #print "### len values:", len(values)
      
      outlier_indice = self._locate_outliers(model, param_name, 
                                             x_coords, y_coords, values, job, worker)
      if len(outlier_indice) > 0:
         new_input_coords, new_input_values = self._remove_outliers(model, param_name, 
                                                                    x_coords, y_coords, values, 
                                                                    outlier_indice, job, worker)  

      #print "### len new values:", len(new_input_coords), len(new_input_values)

      return new_input_coords, new_input_values, len(input_coords) - len(new_input_coords)     
      
   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _locate_outliers(self, model, param_name, x_coords, y_coords, values, job, worker):
      """!
         Locate outliers from @c values of coordinates  @c coords

         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param x_coords coordinate along the x axis 
         @param y_coords coordinate along the y axis 
         @param model parameter values at positions @c x_coords and @c y_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
         @return the indice of the coordinates and values to remove from the data
      """

      outlier_indice = []     # list of outlier indice to return

      # --- Retrieve relevant configuration info
      model_outlier_section_name = "{0}.{1}.{2}".format(model.name.upper(), "OUTLIERS", param_name)
      nb_sectors = math.sqrt(model.config.get_as_float("NUMBER_OF_SECTORS", 
                             model_outlier_section_name))

      # --- Divide the data grid into sectors
      max_x_coords = int(math.ceil(numpy.max(x_coords)))
      max_y_coords = int(math.ceil(numpy.max(y_coords)))
      x_sector_size = int(max_x_coords / nb_sectors)   
      y_sector_size = int(max_y_coords / nb_sectors)   

      for x_size in xrange(0, max_x_coords, x_sector_size):
         for y_size in xrange(0, max_y_coords, y_sector_size):
            
            #print "x_size:", x_size, "y_size:", y_size, "x max:", x_size + x_sector_size, "ymax:", y_size + y_sector_size

            # --- Apply cooordinate filtering on X coords
            coords_indice = numpy.where(numpy.logical_and(
              numpy.logical_and(x_coords >= x_size, x_coords < x_size + x_sector_size),
              numpy.logical_and(y_coords >= y_size, y_coords < y_size + y_sector_size) ) )

            # --- Find indice of outliers if any
            sector_outlier_indice_list = self._find_outliers(model, param_name, 
                                                             values[coords_indice], 
                                                             job, worker)

            # --- Keep a record of outlier indice for subsequent removal
            #     Note: there may be duplicated indice because the sector size may be overestimated
            if len(sector_outlier_indice_list) > 0:
               outlier_indice_list = numpy.take(coords_indice, sector_outlier_indice_list)
               outlier_indice.extend(outlier_indice_list)

      #print "### Nb outlier indice:", len(numpy.unique(outlier_indice))

      return numpy.unique(outlier_indice) 

   # -----------------------------------------------------------------------------------------------
   def _find_outliers(self, model, param_name, values, job, worker):
      """! Locate outliers based on allowed range and/or variance """

      # --- Retrieve relevant configuration info
      model_outlier_section_name = "{0}.{1}.{2}".format(model.name.upper(), "OUTLIERS", param_name)
      median_sigma  = model.config.get_as_float("MEDIAN_ABS_SIGMA", model_outlier_section_name)
      value_range = model.config.get_as_list("ALLOWED_RANGES", model_outlier_section_name)

      #print "*** Initial nb values:", len(values) 


      # --- Check allowed range
      outlier_indice = list(self._check_range(model, param_name, values, value_range, job, worker))

      #print "*** outlier_values after allowed range:", len(outlier_indice)

      # --- Check data variance 
      if median_sigma > 0:
         outlier_indice.extend(list(self._check_variance(model, param_name, values,
                                                                median_sigma, job, worker))) 

      #print "*** outlier_values after check variance:", len(outlier_indice)

      return outlier_indice

   # -----------------------------------------------------------------------------------------------
   def _check_range(self, model, param_name, values, value_range, job, worker):

      outlier_indice = []

      if value_range[0] is None and value_range[1] is not None:
         outlier_indice = list(numpy.ravel(numpy.where(values > value_range[1])))
      elif value_range[0] is not None and value_range[1] is None:
         outlier_indice = list(numpy.ravel(numpy.where(values < value_range[0])))
      elif value_range[0] is not None and value_range[1] is not None:
         outlier_indice = list( numpy.ravel( numpy.where(
                             numpy.logical_or(values < value_range[0], values > value_range[1])) ))

      return outlier_indice


   # -----------------------------------------------------------------------------------------------
   def _check_variance(self, model, param_name, values, nsigma, job, worker):

      predictable_median = numpy.median(values)   
      predictable_std_median = numpy.sqrt(numpy.mean((values - predictable_median)**2))

      predictable_diff_values = numpy.absolute(values - predictable_median)
      #predictable_diff_values = numpy.absolute(values - predictable_mean)

      outlier_indice = numpy.ravel(
                           numpy.where(predictable_diff_values > predictable_std_median * nsigma))
      #outlier_indice = numpy.where(predictable_diff_values > predictable_std_mean * nsigma)   

      #print 'Outliers:', predictable, values[outlier_indice]
      #print predictable, values, values[outlier_indice]
      #print 'median:', predictable_median, 'std:', numpy.std(values), 'std_median:', predictable_std_median

      return outlier_indice

   # -----------------------------------------------------------------------------------------------
   def _remove_outliers(self, model, param_name, 
                              x_coords, y_coords, values, outlier_indice, job, worker):
      """!
         Remove outliers from @c values of coordinates  @c coords

         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param x_coords coordinate along the x axis 
         @param y_coords coordinate along the y axis 
         @param model parameter values at positions @c x_coords and @c y_coords
         @param outlier_indice indice of outliers to remove in coordinates and values 
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
      """

      # Remove outliers
      cleaned_values = numpy.delete(values.copy(), numpy.ravel(outlier_indice))   
         
      new_x_coords = numpy.delete(x_coords.copy(), numpy.ravel(outlier_indice))   
      new_y_coords = numpy.delete(y_coords.copy(), numpy.ravel(outlier_indice))   

      return zip(new_x_coords, new_y_coords), cleaned_values


      


# -- EOF sp_outlier.py
