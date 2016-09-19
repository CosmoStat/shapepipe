"""! 
   Galaxy shape estimator
"""

# -- Python imports
import os, sys
import time
import math
import numpy
from numpy import NaN
import scipy, scipy.signal
import pyfits
from operator import itemgetter
import shutil

# -- External imports
from scatalog import *           # catalog management
from multifit import *           # multifit

# --- Module-specific imports
from gfit_plot import *          # plotter 
from gfit_helper import *        # helper utility functions
from gfit_filter import *        # object filtering


# -------------------------------------------------------------------------------------------------
class GfitShapeEstimator(object):
   
   """! 
      gfit galaxy shape estimator    
   """

   def __init__(self, master):
      """! GfitShapeEstimator constructor """

      self._helper  = GfitHelper()              # helper utility functions
      self._plotter = GfitPlotter()             # plotter

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



   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def measure_galaxy_shapes(self, request, job, worker):
      """!
         Measure the shapes of galaxies contained in the postage stamps of an entire galaxy image 
         file
         @param request shape measurement request
         @param job the GfitJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary with the estimates for relevant features (model parameters and more)
      """

      # --- Check validity of request
      self._check_shape_measurement_request(request, job, worker)

      # --- Process the shape measurement request. The result dicitonary contains 
      #     the measurements and additional information
      result_dico = self._process_shape_measurement_request(request, job, worker)

      return result_dico

   # -----------------------------------------------------------------------------------------------
   def compute_average_ellipticity(self, result_dico, bad_entry_marker):
      """!     
         Compute the mean ellipticity (e1, e2), taking bad entries into account
         @param result_dico result dictionary
         @param bad_entry_marker value used for marking bad ellipticities
      """

      e1_array = numpy.asarray(result_dico["e1"])
      e2_array = numpy.asarray(result_dico["e2"])
      good_e1 = e1_array[e1_array < bad_entry_marker]
      good_e2 = e2_array[e2_array < bad_entry_marker]
      good_e = numpy.hypot(good_e1, good_e2)

      return [numpy.mean(good_e1), numpy.mean(good_e2), numpy.mean(good_e)]

#   # -----------------------------------------------------------------------------------------------
#   def compute_average_result_values(self, result_dico, param_names, bad_entry_marker):

#      param_values = []
#      e1_array = numpy.asarray(result_dico["e1"])
#      good_ei = (e1_array < bad_entry_marker)     # assuming e2 is alwyas marked as bad as well

#      param_values = [numpy.mean(numpy.asarray(result_dico["result"][p])[good_ei]))\
#                                                                              for p in param_names]

#      return param_values
      
   # -----------------------------------------------------------------------------------------------
   class ShapeMeasurementError(Exception):
      """! 
         Exception thrown when a shape measurement error has occurred
      """
      def __init__(self, msg):
         """!
            @param msg error message 
         """
         self._msg = msg

      def __str__(self):
         """! String representation of the ShapeMeasurementError class """
         return "*** ERROR ***: {0}".format(self._msg)


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _check_shape_measurement_request(self, request, job, worker):
      """! Check validity of shape measurement request, paths to image & catalogs in particular """

      # --- Check path to SE catalog file   
      if request.galaxy_se_catalog_filepath is None:
         msg = "{0} - Could not find SExtractor catalog for: /{1}/image-{2:03d}-{3:1d} - {4}"\
             " => Job will not be processed".format(
             worker.name, job.get_branch_tree(), job.img_no, job.epoch, request.image_file_type)  
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      if request.psf_se_catalog_filepath is None:
         msg = "{0} - Could not find SExtractor catalog for: /{1}/image-{2:03d}-{3:1d} - {4}"\
             " => Job will not be processed".format(
             worker.name, job.get_branch_tree(), job.img_no, job.epoch, request.image_file_type)  
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      # --- Check path to PSF image file   
      if request.psf_image_filepath is None:
         msg = "{0} - Could not find PSF image for: /{1}/image-{2:03d}-{3:1d} - {4}"\
                " => Job will not be processed".format(
                worker.name, job.get_branch_tree(), job.img_no, job.epoch, request.image_file_type)
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         


   # -----------------------------------------------------------------------------------------------
   def _process_shape_measurement_request(self, request, job, worker):
      """! Process a shape measurement request """

      # TODO: may be useful to externalize tasks that acts on the postage stamp (e.g. denoising)
      #       instead of hard-wiring such tasks here. 
      # TODO: modularize this code further

      # --- Keep track of execution time
      start_time = time.clock() 


      # --------------------------------------- GALAXY DATA ----------------------------------------

      # --- Output directory where to store the results
      result_output_dir = os.path.join(worker.result_output_dir, job.get_branch_tree())
      error_output_dir = os.path.join(worker.error_output_dir, job.get_branch_tree())

      # --- Get the Galaxy data
      galaxy_image_data, galaxy_image_header = self.helper.read_fits_image(
                                                        request.galaxy_image_filepath, 
                                                        hdu=request.galaxy_image_hdu_no, 
                                                        float_size=request.galaxy_pixel_float_size)  

      # --- Load the galaxy catalog and get all galaxy IDs
      galaxy_id_dico = request.galaxy_id_dico
      if galaxy_id_dico is None:
         if worker.logging_enabled():
            worker.logger.log_error_p(           
               "{0} - Galaxy IDs could not be read - /{1}/image-{2:03d}-{3:1d} - {4}]".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                   request.image_file_type))
         return {}

      # --- Load the SExtractor catalog and get galaxy coordinates and centroids
      galaxy_centroid_dico = request.se_galaxy_info_dico["centroids"]
      galaxy_rel_centroid_dico = request.se_galaxy_info_dico["rel_centroids"]

      # --- Galaxy coordinates sorted by with x varying faster than y
      galaxy_coords = sorted(galaxy_centroid_dico.keys(), key=itemgetter(1, 0))

      # --- 1/2 Stamp_size
      galaxy_half_size = self.helper.get_stamp_center(request.galaxy_stamp_size) + 1

      # ------------------------------------------ PSF DATA ----------------------------------------

      # --- Get the PSF data
      psf_image_data, psf_image_header = self.helper.read_fits_image(
                                                            request.psf_image_filepath, 
                                                            hdu=request.psf_image_hdu_no, 
                                                            float_size=request.psf_pixel_float_size)  

      # --- Load the SExtractor catalog and get PSF coordinates and centroids
      psf_centroid_dico = request.se_psf_info_dico["psf_centroids"]
      psf_rel_centroid_dico = request.se_psf_info_dico["psf_rel_centroids"]

      # --- PSF coordinates sorted by y then x
      psf_coords = sorted(psf_centroid_dico.keys(), key=itemgetter(1, 0))

      # --- 1/2 Stamp_size
      psf_half_size = self.helper.get_stamp_center(request.psf_stamp_size) + 1

      # --- Check for a constant PSF
      psf_stamp = None
      is_constant_psf = (len(psf_coords) == 1)
      if is_constant_psf:
         
         # --- Assume the PSF is constant over the galaxy field =>  extract the PSF once for all
         [psf_stamp, psf_actual_stamp_size, full_psf_stamp] = self._extract_psf_stamp(request, 
                        psf_image_data, psf_coords[0], psf_rel_centroid_dico[psf_coords[0]],
                        psf_half_size, job, worker)

         # --- PSF sky
         psf_sky_noise = self.helper.compute_stamp_sigma_back_noise(full_psf_stamp)
         if request.substract_psf_sky:
            psf_stamp -= psf_sky_noise
            #print "PSF Sky:", psf_sky, self.helper.estimate_sigma_noise(psf_st

         # --- PSF denoising
         psf_sigma_noise = self.helper.compute_stamp_noise(psf_stamp)

         # --- If requested extract postage stamps for debugging
         if __debug__:

            if worker.config.get_as_boolean("EXTRACT_PSF_STAMPS", "DEBUGGING"):
               self._write_stamp_for_debugging(psf_stamp, "PSF", 0, 
                                               psf_rel_centroid_dico[psf_coords[0]],
                                               result_output_dir, job, worker)
            if worker.config.get_as_boolean("PLOT_PSF_STAMPS", "DEBUGGING"):
               plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
               self._plot_stamp_for_debugging(psf_stamp, "PSF", 0, 
                                              request.se_psf_info_dico["psf_centroids"][psf_coords[0]],
                                              plot_output_dir, False, job, worker)

      # ------------------------------------- GALAXY MODEL FITTER ----------------------------------
      
      # --- Instantiate, initialize and configure the GalaxyFitter   
      galaxy_fitter = self._setup_galaxy_fitter(request, job, worker)


      # -------------------------------------- OBJECT FILTERING ------------------------------------

      galaxy_filter = GalaxyFilter(request)

      # ------------------------------------ MAIN PROCESSING LOOP ----------------------------------

      # --- Number of galaxies to process
      total_nb_galaxies = len(galaxy_coords)
      max_nb_galaxies = worker.config.get_as_int("MAX_NB_STAMPS", "DEBUGGING")
      if max_nb_galaxies == -1:
         max_nb_galaxies = total_nb_galaxies 

      # --- Range of galaxies to process
      if worker.config.has_key("OBJECT_RANGE", "DEBUGGING"):
         object_range = worker.config.get_as_list("OBJECT_RANGE", "DEBUGGING")
         if len(object_range) == 0:
            object_range = [-1, -1]

         min_object_no = object_range[0] 
         max_object_no = object_range[1] 

         if min_object_no == -1 or min_object_no is None:    
            min_object_no = 1
         if max_object_no == -1 or max_object_no is None:
            max_object_no = total_nb_galaxies
      else:
         min_object_no = 1
         max_object_no = total_nb_galaxies

      # --- Effective numbers of objects processed
      actual_nb_galaxies = 0    # actual nb. of galaxies sucessfully processed
      failed_nb_galaxies = 0    # nb. of galaxies unsucessfully processed
      warned_nb_galaxies = 0    # nb of galaxy whose extraction raised a warning 
      filtered_nb_galaxies = 0  # nb of filtered galaxies 
      
      result_dico = OrderedDict() # will contain the shape measurement results and other information

      # --- Names & bounds of parameters to fit according to the underlying galaxy model
      param_names  = galaxy_fitter.model.get_param_names()
      param_bounds = galaxy_fitter.model.get_bounds()
   
      # --- Dictionaries where to record results
      result_dico["result"]   = OrderedDict() # values of fitted parameters + gal ID
      result_dico["error"]    = OrderedDict() # values of the fitted parameters + gal ID
      result_dico["filtered"] = OrderedDict() # values of the filtered parameters + gal ID
      for param_name in param_names:
         result_dico["result"][param_name] = []    # values of the fitted parameters + gal ID  
         result_dico["error"][param_name] = []     # always collect error data 
         result_dico["filtered"][param_name] = []  # filtered objects  
      result_dico["result"]["GAL_id"] = []
      result_dico["error"]["GAL_id"] = []         
      result_dico["filtered"]["GAL_id"] = []         
      result_dico["result"]["flag"] = []
      result_dico["error"]["flag"] = []         
      result_dico["filtered"]["flag"] = []         
      result_dico["result"]["weight"] = []

      # --- Flag values to mark entries as either successul (0) failed (1), or filtered (2)
      success_flag  = 0  
      failed_flag   = 1  
      filtered_flag = 2  

      # --- Weighting factors
      success_weight  = 1.0  
      failed_weight   = 0.0  
      filtered_weight = 0.0  

      # --- Create a mosaic of convolved stamps if requested
      
      if worker.config.has_key("CREATE_CONVOLVED_GALAXY_MOSAIC", "DEBUGGING"):
         create_convolved_galaxy_mosaic = worker.config.get_as_boolean(
                                                             "CREATE_CONVOLVED_GALAXY_MOSAIC",
                                                             "DEBUGGING")
      else:
         create_convolved_galaxy_mosaic = False         
                                                         
      if create_convolved_galaxy_mosaic:

         max_common_stamp_size = self._get_max_possible_common_stamp_size(
                                                                      request.se_galaxy_info_dico, 
                                                                      request.galaxy_stamp_size, 
                                                                      job, worker)
         conv_galaxy_mosaic_data, mosaic_stamp_size = self._setup_stamp_mosaic(
                                                            galaxy_image_data,
                                                            request.galaxy_stamp_size,
                                                            request.galaxy_effective_stamp_size,  
                                                            max_common_stamp_size,
                                                            request.galaxy_pixel_float_size,
                                                            job, worker)             
         (mosaic_height, mosaic_width) = conv_galaxy_mosaic_data.shape


      # --- Create a mosaic of residual stamps if requested
      if worker.config.has_key("CREATE_RESIDUALS_GALAXY_MOSAIC", "DEBUGGING"):
         create_residuals_galaxy_mosaic = worker.config.get_as_boolean(
                                                               "CREATE_RESIDUALS_GALAXY_MOSAIC",
                                                               "DEBUGGING")
      else:
         create_residuals_galaxy_mosaic = False

      if create_residuals_galaxy_mosaic:
         max_common_stamp_size = self._get_max_possible_common_stamp_size(
                                                                      request.se_galaxy_info_dico, 
                                                                      request.galaxy_stamp_size, 
                                                                      job, worker)
         res_galaxy_mosaic_data, mosaic_stamp_size = self._setup_stamp_mosaic(
                                                         galaxy_image_data, 
                                                         request.galaxy_stamp_size,
                                                         request.galaxy_effective_stamp_size,
                                                         max_common_stamp_size,
                                                         request.galaxy_pixel_float_size,
                                                         job, worker)    
         (mosaic_height, mosaic_width) = res_galaxy_mosaic_data.shape

      # --- Tell if we keep a record of failed fits in the errors directory
      record_failed_fits = worker.config.get_as_boolean("RECORD_FAILED_FITS", "DEBUGGING")

      # --- Progress count
      progress_count = worker.config.get_as_int("PROGRESS_COUNT", "DEBUGGING")


      logging_enabled = worker.logging_enabled()

      if logging_enabled:
         # --- Stamp size
         worker.logger.log_info_p(
            "{0} - /{1}/image-{2:03d}-{3:1d} - "\
            "Original stamp size: {4} - Effective stamp size: {5}".format(
                worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                request.galaxy_stamp_size, request.galaxy_effective_stamp_size))
   
         # --- How many objects to process?
         if total_nb_galaxies > max_nb_galaxies: 
            worker.logger.log_info_p(
               "{0} - /{1}/image-{2:03d}-{3:1d} - Processing {4}/{5} galaxies".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                   max_nb_galaxies, total_nb_galaxies))
         else:
            worker.logger.log_info_p(
               "{0} - /{1}/image-{2:03d}-{3:1d} - Processing {4} galaxies".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, total_nb_galaxies))

      if logging_enabled:
         worker.logger.log_info("=" * 165)

      # --- Iterate over all galaxy stamps and measure the corresponding galaxy shapes 

      igal = 1
      for (x,y) in galaxy_coords:

         # --- Skip objects according to requested object range (if any)   
         if igal < min_object_no:
            if igal > max_nb_galaxies or igal > max_object_no:
               break
            else:
               continue

         # --- Extract a galaxy postage stamp of the requested size around the centroid
         [galaxy_stamp, galaxy_actual_stamp_size, full_galaxy_stamp] =\
                 self._extract_galaxy_stamp(request, 
                                     galaxy_image_data, (x,y), galaxy_rel_centroid_dico[(x,y)], 
                                     galaxy_half_size, job, worker)
         
         if galaxy_stamp is None:
            igal += 1
            if igal > max_nb_galaxies or igal > max_object_no:
               break
            else:
               continue

         # --- Sky
         gal_sky_noise = self.helper.compute_stamp_sigma_back_noise(full_galaxy_stamp) 
         if request.substract_galaxy_sky:
            galaxy_stamp -= gal_sky_noise

         # --- Galaxy denoising

         # --- Object Noise and Denoising

         gal_sigma_noise = self.helper.compute_stamp_noise(galaxy_stamp)
         full_gal_sigma_noise = self.helper.compute_stamp_noise(full_galaxy_stamp)

         #print "Galaxy noise:", gal_sigma_noise, full_gal_sigma_noise 
         #print "gal_sigma_noise:",  gal_sigma_noise, "full gal_sigma_noise",  full_gal_sigma_noise, 
         #      "sky:", self.helper.compute_stamp_sigma_back_noise(galaxy_stamp)

         # --- Extract a PSF postage stamp of the requested size around the centroid
         if not is_constant_psf:

            [psf_stamp, psf_actual_stamp_size, full_psf_stamp] = self._extract_psf_stamp(request, 
                                               psf_image_data, (x, y), psf_rel_centroid_dico[(x,y)], 
                                               psf_half_size, job, worker)

            if psf_stamp is None:
               igal += 1
               if igal > max_nb_galaxies or igal > max_object_no:
                  break
               else:
                  continue

            # --- PSF sky
            psf_sky_noise = self.helper.compute_stamp_sigma_back_noise(full_psf_stamp)
            if request.substract_psf_sky:
               psf_stamp -= psf_sky_noise
               #print "PSF Sky:", psf_sky, self.helper.estimate_sigma_noise(psf_stamp)

            # --- PSF denoising
            psf_sigma_noise = self.helper.compute_stamp_noise(psf_stamp)


            # --- If requested extract postage stamps for debugging
            if __debug__:
               if worker.config.get_as_boolean("EXTRACT_PSF_STAMPS", "DEBUGGING"):
                  #print "EXTRACTING STAMP", igal, x,y, psf_centroid_dico[(x,y)] 
                  self._write_stamp_for_debugging(psf_stamp, "PSF", igal, 
                                                                   psf_centroid_dico[(x,y)],
                                                                   result_output_dir, job, worker)
               if worker.config.get_as_boolean("PLOT_PSF_STAMPS", "DEBUGGING"):
                  plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
                  self._plot_stamp_for_debugging(psf_stamp, "PSF", igal, 
                                         request.se_psf_info_dico["psf_centroids"][psf_coords[0]],
                                         plot_output_dir, False, job, worker)
                  enable_3D = worker.config.get_as_boolean("PLOT_ENABLE_3D", "DEBUGGING")
                  if enable_3D:
                     self._plot_stamp_for_debugging(psf_stamp, "PSF", igal, 
                                         request.se_psf_info_dico["psf_centroids"][psf_coords[0]],
                                         plot_output_dir, True, job, worker)

         # --- If requested extract postage stamps for debugging
         if __debug__:         
            if worker.config.get_as_boolean("EXTRACT_GALAXY_STAMPS", "DEBUGGING"):
               self._write_stamp_for_debugging(galaxy_stamp, "galaxy", igal,
                                               galaxy_centroid_dico[(x,y)],
                                               result_output_dir, job, worker)

            if worker.config.get_as_boolean("PLOT_GALAXY_STAMPS", "DEBUGGING"):
               plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
               enable_3D = worker.config.get_as_boolean("PLOT_ENABLE_3D", "DEBUGGING")
                     
               self._plot_stamp_for_debugging(galaxy_stamp, "galaxy", igal, 
                                              request.se_galaxy_info_dico["centroids"][(x,y)],
                                              plot_output_dir, False, job, worker)
               if enable_3D:
                  self._plot_stamp_for_debugging(galaxy_stamp, "galaxy", igal, 
                                              request.se_galaxy_info_dico["centroids"][(x,y)],
                                              plot_output_dir, True, job, worker)
         
         # --- Check the stamp sizes
         if galaxy_actual_stamp_size < request.galaxy_effective_stamp_size or \
            psf_actual_stamp_size < request.psf_effective_stamp_size:

            # Update request data and increase warning count
            #request.galaxy_actual_stamp_size = galaxy_actual_stamp_size
            #request.psf_actual_stamp_size = psf_actual_stamp_size
            
            warned_nb_galaxies += 1

         # --- Normalize the galaxy image ?
         if request.normalize_galaxy:
            galaxy_stamp /= numpy.sum(galaxy_stamp)
            full_galaxy_stamp /= numpy.sum(full_galaxy_stamp)

         # --- Normalize the PSF image ?
         if request.normalize_psf:
            psf_stamp /= numpy.sum(psf_stamp)

         # --- Measure the galaxy shape
         [guess_params, fitted_params, status] = self._measure_galaxy_shape(
                                                                   galaxy_fitter, request, igal,
                                                                   (x,y), 
                                                                   galaxy_stamp, full_galaxy_stamp,
                                                                   psf_stamp, full_psf_stamp,
                                                                   gal_sigma_noise, psf_sigma_noise,
                                                                   job, worker)

         # --- Check Status and parameter validity
         if status.errno >= 0 and \
           galaxy_fitter.model.check_fitted_params((x, y), galaxy_stamp,
                                                           request.galaxy_info_dico, 
                                                           request.se_galaxy_info_dico, 
                                                           guess_params, fitted_params):

            # --- Apply filtering if requested
            filtered = False
            if request.enable_galaxy_filtering:

               filtered = galaxy_filter.match((x,y), status.custom_data, galaxy_stamp, job, worker)
               if filtered:

                  # --- Record filtered data
                  for (param_name, param_value) in zip(param_names, fitted_params):               
                     result_dico["filtered"][param_name].append(param_value)    
                  result_dico["filtered"]["GAL_id"].append(galaxy_id_dico[(x,y)])   
                  result_dico["filtered"]["flag"].append(filtered_flag) 
                  
                  # --- Apply filtering policy, i.e. decide what to do with the filtered objects
                  #     - remove the objects alltogether, i.e. mack them as failed (ei=10)
                  #     - apply a weighting scheme to minimize the influence of such objects
                  ###self._apply_galaxy_filtering_policy(result_dico, request, job, worker)               

                  # --- For now, we mark as "failed" (i.e. "remove") filtered objects
                  if request.galaxy_filtering_policy in ["remove", "weight"]:

                     for (param_name, param_value) in zip(param_names, fitted_params):
                        if param_name in ["e1", "e2"]:
                           result_dico["result"][param_name].append( 
                                                               request.failed_ellipticity_value)
                        else:            
                           result_dico["result"][param_name].append(param_value)
                     result_dico["result"]["GAL_id"].append(galaxy_id_dico[(x,y)]) 
                     result_dico["result"]["flag"].append(filtered_flag) 
                     result_dico["result"]["weight"].append(filtered_weight) 

                  # --- Record additional quantities besides model parameters
                  self._collect_extra_data(result_dico["filtered"], request, (x,y), igal,
                                           full_galaxy_stamp, galaxy_stamp, 
                                           status.custom_data, job, worker)

                  filtered_nb_galaxies += 1  # one more filtered object

            if not filtered:

               # --- Fitting successful and not filtered...
               actual_nb_galaxies += 1       # update nb. of successully processed galaxies

               # --- Update the "result" section of the result dictionary
               for (param_name, param_value) in zip(param_names, fitted_params):
                  result_dico["result"][param_name].append(param_value)    

               result_dico["result"]["GAL_id"].append(galaxy_id_dico[(x,y)]) 
               result_dico["result"]["flag"].append(success_flag) 
               result_dico["result"]["weight"].append(success_weight) 

               try:   

                  # --- Record fitted stanmps and/or residuals if requested
                  if create_convolved_galaxy_mosaic or create_residuals_galaxy_mosaic:

                     # --- Equivalent object coordinates
                     xm = x / request.galaxy_stamp_size * mosaic_stamp_size 
                     ym = y / request.galaxy_stamp_size * mosaic_stamp_size 

                     # --- Just record (model - estimate); 
                     #     <residuals> contains (model - estimate)^2/<sigma-noise>^2

                     # --- Record fitted convolved stamps in a mosaic if requested
                     if create_convolved_galaxy_mosaic:

                        if not filtered and not status.custom_data is None:
                           [convolved_galaxy_stamp, residuals] = status.custom_data[-1]

                           rel_xc = rel_yc = mosaic_stamp_size / 2.0  
                           new_fitted_galaxy_stamp = self.helper.cut_stamp_around_centroid(
                                               convolved_galaxy_stamp,
                                               mosaic_stamp_size,
                                               (rel_yc, rel_xc))

                           conv_galaxy_mosaic_data[ym:ym+mosaic_stamp_size,\
                                                 xm:xm+mosaic_stamp_size] = new_fitted_galaxy_stamp


                     # --- Record residual stamps in a mosaic if requested   
                     if create_residuals_galaxy_mosaic:
                        if not filtered and not status.custom_data is None:
                           [convolved_galaxy_stamp, residuals] = status.custom_data[-1]
                           galaxy_residuals = numpy.sqrt(residuals * gal_sigma_noise**2)

                           rel_xc = rel_yc = mosaic_stamp_size / 2.0  
                           residuals_stamp = self.helper.cut_stamp_around_centroid(
                                                                                 galaxy_residuals, 
                                                                                 mosaic_stamp_size,  
                                                                                 (rel_yc, rel_xc)) 
                           res_galaxy_mosaic_data[ym:ym+mosaic_stamp_size,\
                                                  xm:xm+mosaic_stamp_size] = residuals_stamp


               except Exception as detail:
                  if worker.logging_enabled():
                     worker.logger.log_warning_p(
                               "{0} - Exception while populating mosaic: {1} ({2})".format(
                                                                                       worker.name, 
                                                                                       job,
                                                                                       detail))
            if logging_enabled:
               guess_param_str  = ["{0}={1:.9f}".format(n, v) for (n,v) \
                                                              in zip(param_names, guess_params)]
               fitted_param_str = ["{0}={1:.9f}".format(n, v) for (n,v) \
                                                              in zip(param_names, fitted_params)]
               if not filtered:
                  worker.logger.log_info(
                     "{0:>5d} - Guess Params : {1} \n {2} Fitted Params: {3}".format(igal, 
                         str(["{0}={1:+.9f}".format(n, v) for (n,v) \
                            in zip(param_names, guess_params)]).replace(",","").replace("'",""),
                         " "*15,
                         str(["{0}={1:+.9f}".format(n, v) for (n,v) \
                            in zip(param_names, fitted_params)]).replace(",","").replace("'","") ) )
               else:
                  worker.logger.log_info(
                     "{0:>5d} - *** Filtered ***\n                 Guess Params : {1} \n {2} "\
                     "Fitted Params: {3}".format(
                       igal, str(["{0}={1:+.9f}".format(n, v) for (n,v)\
                       in zip(param_names, guess_params)]).replace(",","").replace("'",""), " "*15,
                       str(["{0}={1:+.9f}".format(n, v) for (n,v) in zip(param_names, 
                                               fitted_params)]).replace(",","").replace("'","") ) )
               worker.logger.log_info("-" * 165)


         else:

            # --- Fitting failed

            failed_nb_galaxies += 1     # update nb. of failuress

            # --- Update the "result" section of the result dictionary
            for (param_name, param_value) in zip(param_names, fitted_params):
               if param_name in ["e1", "e2"]:
                  if request.failed_ellipticity_value != -1:
                     result_dico["result"][param_name].append(request.failed_ellipticity_value)
                  else:
                     result_dico["result"][param_name].append(param_value)                   
               else:            
                  result_dico["result"][param_name].append(param_value)    
            result_dico["result"]["GAL_id"].append(galaxy_id_dico[(x,y)]) 
            result_dico["result"]["flag"].append(failed_flag)
            result_dico["result"]["weight"].append(failed_weight) 
 

            # --- Keep track of failed measurements for later analysis
            for (param_name, param_value) in zip(param_names, fitted_params):
               result_dico["error"][param_name].append(param_value)    
            result_dico["error"]["GAL_id"].append(galaxy_id_dico[(x,y)]) 
            result_dico["error"]["flag"].append(failed_flag) 

            # --- Record additional quantities besides model parameters
            self._collect_extra_data(result_dico["error"], request, (x,y), igal,
                                     full_galaxy_stamp, galaxy_stamp, 
                                     status.custom_data, job, worker)

            # --- Keep a record of the galaxy in the errors directory
            if record_failed_fits:
               self._write_stamp_for_debugging(galaxy_stamp, "galaxy", igal,
                                               galaxy_centroid_dico[(x,y)],
                                               error_output_dir, job, worker)
               self._plot_stamp_for_debugging(galaxy_stamp, "galaxy", igal, 
                                              request.se_galaxy_info_dico["centroids"][(x,y)],
                                              error_output_dir, False, job, worker)
            
            if logging_enabled:
               guess_param_str  = ["{0}={1:.9f}".format(n, v) for (n,v) \
                                                              in zip(param_names, guess_params)]
               fitted_param_str = ["{0}={1:.9f}".format(n, v) for (n,v) \
                                                              in zip(param_names, fitted_params)]
               worker.logger.log_info(
                  "{0:>5d} - *** Fitting FAILED ***\n                 Guess Params : {1} \n {2} "\
                  "Fitted Params: {3}".format(igal, 
                      str(["{0}={1:+.9f}".format(n, v) for (n,v) in zip(param_names, 
                                         guess_params)]).replace(",","").replace("'",""), " "*15,
                      str(["{0}={1:+.9f}".format(n, v) for (n,v) in zip(param_names, 
                                         fitted_params)]).replace(",","").replace("'","") ) )
               worker.logger.log_info("-" * 165)
           
         # --- Endif status

         # --- Keep track of fitted parameters and other quantities
         self._collect_extra_data(result_dico["result"], request, (x,y), igal,
                                  full_galaxy_stamp, galaxy_stamp, 
                                  status.custom_data, job, worker)

         # --- Display progress
         if igal % progress_count == 0:
            print(
              "{0} - /{1}/image-{2:03d}-{3:1d} - "\
              "Galaxies processed: {4}/{5} ({6}-{7}-{8}-{9})".format(
              worker.name, job.get_branch_tree(), job.img_no, job.epoch, igal, max_nb_galaxies,
              actual_nb_galaxies, warned_nb_galaxies, failed_nb_galaxies, filtered_nb_galaxies) )      
            worker.logger.flush()      

         # --- Check if we are done
         igal += 1         
         if igal > max_nb_galaxies or igal > max_object_no:
            break


      # --- end for (x,y)

#      branch_key = (job.branch, job.obs_type, job.data_type)
#      image_key = (job.img_no, job.epoch)
#      print "COLLECTED:", data_collector.query_data_for(variable="flux", branch_key=branch_key, image_key=image_key)
#      print "COLLECTED:", data_collector.query_data_for(variable="flux", branch_key=branch_key)

      # --- Update the "information" section of the result dictionary
      result_dico["info"] = {}
      result_dico["info"]["success_count"]  = actual_nb_galaxies
      result_dico["info"]["failure_count"]  = failed_nb_galaxies
      result_dico["info"]["warning_count"]  = warned_nb_galaxies
      result_dico["info"]["filtered_count"] = filtered_nb_galaxies
      result_dico["info"]["total_count"]    = total_nb_galaxies
      result_dico["info"]["elapsed_time"]   = time.clock() - start_time

      # --- Information about the fitted model, methods and fitted parameters
      result_dico["model"] = {}
      result_dico["model"]["model_name"]   = galaxy_fitter.method.name
      result_dico["model"]["method_name"]  = galaxy_fitter.model.name
      result_dico["model"]["param_names"]  = param_names
      result_dico["model"]["param_bounds"] = param_bounds


      # --- Collect unprocessed galaxy IDs (not in processed SExtractor but in galaxy catalog)
      result_dico["info"]["unprocessed_ids"] = self._find_unprocessed_values(
                                                                   result_dico["result"]["GAL_id"],
                                                                   request.galaxy_id_dico.values()) 
      # --- Collect unprocessed field positions
      if "GAL_Xc" in result_dico["result"] and "GAL_Yc" in result_dico["result"]:

         Xc_Yc_coords = zip(result_dico["result"]["GAL_Xc"], 
                                  result_dico["result"]["GAL_Yc"])   

         unprocessed_Xc_Yc_coords = self._find_unprocessed_values(
                                                    Xc_Yc_coords, 
                                                    request.galaxy_info_dico["GAL_Xc_Yc"].values())
      else:
         unprocessed_Xc_Yc_coords = []

      if len(unprocessed_Xc_Yc_coords) > 0:	
         result_dico["info"]["unprocessed_GAL_Xc"],\
         result_dico["info"]["unprocessed_GAL_Yc"] = zip(*unprocessed_Xc_Yc_coords) 
      else:
         result_dico["info"]["unprocessed_GAL_Xc"] = []
         result_dico["info"]["unprocessed_GAL_Yc"] = []

      if "GAL_RA" in result_dico["result"] and \
         "GAL_DEC" in result_dico["result"]:
         RA_DEC_coords = zip(result_dico["result"]["GAL_RA"], 
                                     result_dico["result"]["GAL_DEC"])   
         unprocessed_RA_DEC_coords = self._find_unprocessed_values(
                                                   RA_DEC_coords, 
                                                   request.galaxy_info_dico["GAL_RA_DEC"].values())
      else:
         unprocessed_RA_DEC_coords = []

      if len(unprocessed_RA_DEC_coords) > 0:
         result_dico["info"]["unprocessed_GAL_RA"],\
         result_dico["info"]["unprocessed_GAL_DEC"] = zip(*unprocessed_RA_DEC_coords) 
      else:
         result_dico["info"]["unprocessed_GAL_RA"] = []
         result_dico["info"]["unprocessed_GAL_DEC"] = []

      # --- Output summary statistics 
      if result_dico["info"]["success_count"] > 0:

         if logging_enabled:
            good_ei = (numpy.asarray(result_dico["result"]["e1"]) < request.failed_ellipticity_value)
            fitted_param_str = str(
               ["<{0}>={1:.15f}".format(
                    p, numpy.mean(numpy.asarray(result_dico["result"][p])[good_ei])) \
                                                                           for p in param_names])
            worker.logger.log_info_p("{0}".format(
                                                fitted_param_str.replace(",","").replace("'","")) )

            worker.logger.log_info("-" * 165)
            msg = "{0} - /{1}/img {2:03}-{3:1d} - Galaxies processed: {4}/{5} ({6}-{7}-{8}-{9})"
            worker.logger.log_info_p(msg.format(
               worker.name, job.get_branch_tree(), job.img_no, job.epoch,
               result_dico["info"]["success_count"], result_dico["info"]["total_count"], 
               result_dico["info"]["success_count"], result_dico["info"]["warning_count"],
               result_dico["info"]["failure_count"], result_dico["info"]["filtered_count"]))

            worker.logger.log_info("=" * 165)

      # --- Mosaics
      if create_convolved_galaxy_mosaic:

         # --- Save the mosaic with convolved stamps
         conv_galaxy_mosaic_filename = "image-{0:03}-{1:1d}.fits".format(
                                                                          job.img_no, job.epoch)
         conv_galaxy_mosaic_filepath = os.path.join(result_output_dir, 
                                                    conv_galaxy_mosaic_filename)
         source_image_header = pyfits.getheader(request.galaxy_image_filepath)
         self.helper.write_as_fits(conv_galaxy_mosaic_data, conv_galaxy_mosaic_filepath, 
                                   header=source_image_header) 

      if create_residuals_galaxy_mosaic:

         # --- Save the mosaic with convolved stamps
         res_galaxy_mosaic_filename = "residuals_galaxy_mosaic-{0:03}-{1:1d}.fits".format(
                                                                          job.img_no, job.epoch)
         res_galaxy_mosaic_filepath = os.path.join(result_output_dir, 
                                                    res_galaxy_mosaic_filename)
         source_image_header = pyfits.getheader(request.galaxy_image_filepath)
         self.helper.write_as_fits(res_galaxy_mosaic_data, res_galaxy_mosaic_filepath, 
                                   header=source_image_header) 

      return result_dico


   # -----------------------------------------------------------------------------------------------
   def _measure_galaxy_shape(self, fitter, request, igal, (x,y), galaxy_stamp, full_galaxy_stamp, 
                                                           psf_stamp, full_psf_stamp,
                                                           gal_sigma_noise, psf_sigma_noise,
                                                           job, worker):

#      # --- Substract background noise
#      if request.substract_galaxy_sky:
#         galaxy_sky = self.helper.compute_stamp_sigma_back_noise(galaxy_stamp)
#         #galaxy_sky = self.helper.compute_stamp_sky_background(full_galaxy_stamp)
#         galaxy_stamp -= galaxy_sky
#         #print "Galaxy Sky:", galaxy_sky, self.helper.estimate_sigma_noise(galaxy_stamp)

#      if request.substract_psf_sky:
#         psf_sky = self.helper.compute_stamp_sigma_back_noise(psf_stamp)
#         #psf_sky = self.helper.compute_stamp_sky_background(full_psf_stamp)
#         psf_stamp -= psf_sky
#         #print "PSF Sky:", psf_sky, self.helper.estimate_sigma_noise(psf_stamp)

      #print "gal_sigma_back_noise :", gal_sigma_back_noise , "psf_sigma_back_noise", psf_sigma_back_noise

      #full_gal_SNR = self.helper.compute_stamp_SNR(full_galaxy_stamp)
      #gal_SNR = self.helper.compute_stamp_SNR(galaxy_stamp)

      #full_gal_noise = self.helper.compute_stamp_noise(full_galaxy_stamp)
      #gal_sigma_noise = self.helper.compute_stamp_noise(galaxy_stamp)
      #psf_sigma_noise = self.helper.compute_stamp_noise(psf_stamp)

      sigma_noise_sq = gal_sigma_noise**2 + psf_sigma_noise**2
      #sigma_noise_sq = psf_sigma_back_noise**2 + gal_sigma_back_noise**2

      ####### TEMP #######
      ### sigma_noise_sq = 1.0
      ####### TEMP #######

      # --- Objective function: we use the Chi2 template from multifit
      fitter.method.objective_func = fitter.method.get_chi2_objective_func()

      # --- Least-Squares regularization term    
      #rlambda = request.regularization_lambda

      # --- Arguments to pass to the objective function
      #     Note: the first 4 parameters are mandatory
      args = [self._create_galaxy_model, fitter.model, galaxy_stamp,  sigma_noise_sq,
              [(x,y), galaxy_stamp, psf_stamp, request, job, worker]]

      # --- Ask the model to estimate the parameter guess values
      fitter.model.reset_params_from_config()
      guess_params = fitter.model.get_best_guess_params( (x,y), galaxy_stamp, 
                                                         request.galaxy_info_dico, 
                                                         request.se_galaxy_info_dico)

      # --- Setup diagnostic options for debugging
      if __debug__:
         self._setup_diag_options(fitter, igal, job, worker)

      # --- Fit the galaxy from the postage stamp image
      status = fitter.fit(args)

      # --- Fitted parameters
      fitted_params = fitter.fitted_params

      #--- Residuals
      residuals = status.residuals

      # --- Check for invalid fits
      if not status.custom_data is None:

         [convolved_galaxy_stamp, residuals] = status.custom_data[-1]
         if not convolved_galaxy_stamp is None:

            # --- Check validity of the model data
            if type(convolved_galaxy_stamp) is float or numpy.any(numpy.isnan(
                                                                          convolved_galaxy_stamp)):
               status.errno = -1 # => make the fit fail        

            if status.errno != -1:
               convolved_galaxy_stamp = convolved_galaxy_stamp.reshape(galaxy_stamp.shape)
            else:
               convolved_galaxy_stamp = None

         if not residuals is None:
            residuals = residuals.reshape(galaxy_stamp.shape)

      if __debug__:

         # --- Write residuals to disk if requested
         if not residuals is None:  
            if worker.config.get_as_boolean("EXTRACT_GALAXY_RESIDUALS", "DEBUGGING"):
                  result_output_dir = os.path.join(worker.result_output_dir, job.get_branch_tree())
                  self._write_stamp_for_debugging(residuals, "residuals", igal, 
                                                  request.se_galaxy_info_dico["centroids"][(x,y)],
                                                  result_output_dir, job, worker)

            # --- Plot residuals to disk if requested
            if not residuals is None:  
               if worker.config.get_as_boolean("PLOT_GALAXY_RESIDUALS", "DEBUGGING"):
                  plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
                  enable_3D = worker.config.get_as_boolean("PLOT_ENABLE_3D", "DEBUGGING")
                  self._plot_stamp_for_debugging(residuals, "residuals", igal, 
                                                 request.se_galaxy_info_dico["centroids"][(x,y)],
                                                 plot_output_dir, False, job, worker)
                  if enable_3D:
                     self._plot_stamp_for_debugging(residuals, "residuals", igal, 
                                                    request.se_galaxy_info_dico["centroids"][(x,y)],
                                                    plot_output_dir, True, job, worker)            

         # --- Write model data to disk if requested
         if not status.custom_data is None:

            if not convolved_galaxy_stamp is None:  
               if worker.config.get_as_boolean("EXTRACT_CONVOLVED_GALAXY_STAMPS", 
                                                          "DEBUGGING"):
                  result_output_dir = os.path.join(worker.result_output_dir, job.get_branch_tree())
                  self._write_stamp_for_debugging(convolved_galaxy_stamp, "convolved_galaxy", igal, 
                                                  request.se_galaxy_info_dico["centroids"][(x,y)],
                                                  result_output_dir, job, worker)
         

            # --- Plot model data to disk if requested
            if not convolved_galaxy_stamp is None:  
               if worker.config.get_as_boolean("PLOT_CONVOLVED_GALAXY_STAMPS", "DEBUGGING"):
                  plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
                  enable_3D = worker.config.get_as_boolean("PLOT_ENABLE_3D", "DEBUGGING")
                  self._plot_stamp_for_debugging(convolved_galaxy_stamp, "convolved_galaxy", igal, 
                                                 request.se_galaxy_info_dico["centroids"][(x,y)],
                                                 plot_output_dir, False, job, worker)
                  if enable_3D:
                     self._plot_stamp_for_debugging(convolved_galaxy_stamp, "convolved_galaxy", igal, 
                                                    request.se_galaxy_info_dico["centroids"][(x,y)],
                                                    plot_output_dir, True, job, worker)

      return [[p.guess_value for p in guess_params], fitted_params, status]


   # -----------------------------------------------------------------------------------------------
   def _collect_extra_data(self, data_dico, request, (x,y), igal,
                                 full_galaxy, galaxy, custom_data, job, worker):   

      # --- Galaxy catalog data
      gal_data_dico = request.galaxy_info_dico

      # Coordinates in galaxy catalog
      if not "GAL_x" in data_dico:
         data_dico["GAL_x"] = []   
      if not "GAL_y" in data_dico:
         data_dico["GAL_y"] = []   
      data_dico["GAL_x"].append(x)
      data_dico["GAL_y"].append(y)

      # Physical coordinates
      (gal_Xc, gal_Yc) = gal_data_dico["GAL_Xc_Yc"][(x,y)]
      if not "GAL_Xc" in data_dico:
         data_dico["GAL_Xc"] = []   
      if not "GAL_Yc" in data_dico:
         data_dico["GAL_Yc"] = []   
      data_dico["GAL_Xc"].append(gal_Xc)
      data_dico["GAL_Yc"].append(gal_Yc)
      
      # RA and DEC if available
      if "GAL_RA_DEC" in gal_data_dico:   
         (gal_RA, gal_DEC) = gal_data_dico["GAL_RA_DEC"][(x,y)]
         if not "GAL_RA" in data_dico:
            data_dico["GAL_RA"] = []   
         if not "GAL_DEC" in data_dico:
            data_dico["GAL_DEC"] = []   
         data_dico["GAL_RA"].append(gal_RA)
         data_dico["GAL_DEC"].append(gal_DEC)

      # Extra columns
      excluded_columns = ["GAL_id", "GAL_xc_yc","GAL_Xc_Yc", "GAL_RA_DEC"]       
      for colname in gal_data_dico.keys():
         if not colname in excluded_columns:
            if not colname in data_dico:
               data_dico[colname] = []
            data_dico[colname].append(gal_data_dico[colname][(x,y)]) 
           
      # --- Galaxy Sextractor catalog data
      se_gal_data_dico = request.se_galaxy_info_dico

      # Galaxy centroids estimated by SExtractor
      (se_gal_xc, se_gal_yc) = se_gal_data_dico["centroids"][(x,y)]
      if not "SE_GAL_xc" in data_dico:
         data_dico["SE_GAL_xc"] = []   
      if not "SE_GAL_yc" in data_dico:
         data_dico["SE_GAL_yc"] = []   
      data_dico["SE_GAL_xc"].append(se_gal_xc)
      data_dico["SE_GAL_yc"].append(se_gal_yc)

      # Galaxy relative centroids estimated by SExtractor
      (se_gal_rel_xc, se_gal_rel_yc) = se_gal_data_dico["rel_centroids"][(x,y)]
      if not "SE_GAL_rel_xc" in data_dico:
         data_dico["SE_GAL_rel_xc"] = []   
      if not "SE_GAL_rel_yc" in data_dico:
         data_dico["SE_GAL_rel_yc"] = []   
      data_dico["SE_GAL_rel_xc"].append(se_gal_rel_xc)
      data_dico["SE_GAL_rel_yc"].append(se_gal_rel_yc)

      # Extra SExtractor columns
      se_excluded_columns = ["SE_GAL_X_IMAGE", "SE_GAL_Y_IMAGE", "centroids", "rel_centroids"]
#      se_excluded_columns = ["SE_GAL_X_IMAGE", "SE_GAL_Y_IMAGE", "centroids", "rel_centroids", \
#                             "SE_GAL_rel_xc", "SE_GAL_rel_yc"]
      for se_key in se_gal_data_dico.keys():
         if not se_key in se_excluded_columns:
            if not se_key in data_dico:
               data_dico[se_key] = []
            data_dico[se_key].append(se_gal_data_dico[se_key][(x,y)])

      # --- PSF SExtractor catalog data (assume PSF are at galaxy positions after interpolation)
      se_psf_data_dico = request.se_psf_info_dico
      
      if len(se_psf_data_dico["psf_centroids"]) > 1:

         # PSF centroids estimated by SExtractor
         (se_psf_xc, se_psf_yc) = se_psf_data_dico["psf_centroids"][(x,y)]
         if not "SE_PSF_xc" in data_dico:
            data_dico["SE_PSF_xc"] = []   
         if not "SE_PSF_yc" in data_dico:
            data_dico["SE_PSF_yc"] = []   
         data_dico["SE_PSF_xc"].append(se_psf_xc)
         data_dico["SE_PSF_yc"].append(se_psf_yc)

         # Extra SExtractor columns
         se_excluded_columns = ["SE_PSF_X_IMAGE", "SE_PSF_Y_IMAGE",\
                                "psf_centroids", "psf_rel_centroids"]       
         for se_key in se_psf_data_dico.keys():
            if not se_key in se_excluded_columns:
               if not se_key in data_dico:
                  data_dico[se_key] = []
               data_dico[se_key].append(se_psf_data_dico[se_key][(x,y)])

      else:
         # --- This may be a unique constant PSF, we don't track it here
         pass

      # --- Galaxy number (useful for plotting)
      if not "iobj" in data_dico:
         data_dico["iobj"] = []
      data_dico["iobj"].append(igal)

      # --- Actual stamp size
      if not "gal_stamp_size" in data_dico:
         data_dico["gal_stamp_size"] = []
      data_dico["gal_stamp_size"].append(galaxy.shape[0])

    
      # ---Record data that depends on the residuals, if they have been computed
      #    *** Note: if the fit has failed, <convolved_galaxy> may be NaN or and the <residuals> 
      #              may be nans. We nevertheless record the data.
      qvars = []
      if not custom_data is None and len(custom_data) > 0:  

         # --- Compute normalized residuals (per pixel)
         [convolved_galaxy, residuals] = custom_data[-1]
         norm_residuals = residuals / (residuals.shape[0]*residuals.shape[1]) 
#         norm_residuals_12_12 = self.helper.cut_stamp_around_centroid( 
#                                        norm_residuals, 12, (residuals.shape[0]/2.0, 
#                                                             residuals.shape[1]/2.0))
#         norm_residuals_16_16 = self.helper.cut_stamp_around_centroid( 
#                                        norm_residuals, 16, (residuals.shape[0]/2.0, 
#                                                             residuals.shape[1]/2.0))
#         norm_residuals_20_20 = self.helper.cut_stamp_around_centroid( 
#                                        norm_residuals, 20, (residuals.shape[0]/2.0, 
#                                                             residuals.shape[1]/2.0))
         ### norm_residuals_2 = residuals / residuals.sum() 

         qvars = ["norm_residuals"]
         #qvars = ["norm_residuals", \
         #         "norm_residuals_12_12", "norm_residuals_16_16", "norm_residuals_20_20"]
         ### qvars = ["norm_residuals", "norm_residuals_2"]
         qvals = [eval(q) for q in qvars]
         self.helper.record_data_items(data_dico, qvars, qvals)

 
      # *** Note: Depending on the minimizer, the residuals may not be available for analysis.
      #           In such a case, <custom_data> is None.   


   # -----------------------------------------------------------------------------------------------
   def _apply_galaxy_filtering_policy(self, result_dico, request, job, worker):

      if request.galaxy_filtering_policy in ["remove", "weight"]:

         for (param_name, param_value) in zip(param_names, fitted_params):
            if param_name in ["e1", "e2"]:
               if request.failed_ellipticity_value != -1:
                  result_dico["result"][param_name].append(request.failed_ellipticity_value)
               else:
                  result_dico["result"][param_name].append(param_value)                   
            else:            
               result_dico["result"][param_name].append(param_value)



         for param_name in ["e1", "e2"]:
            result_dico["result"][param_name].append(request.failed_ellipticity_value)

      elif request.galaxy_filtering_policy == "weight":
         # TODO: for now we do the same as for "remove"
         for param_name in ["e1", "e2"]:
            result_dico["result"][param_name].append(request.failed_ellipticity_value)
         

   # -----------------------------------------------------------------------------------------------
   def _get_max_possible_common_stamp_size(self, se_info_dico, stamp_size, job, worker):
      """! Return the biggest possible common stamp size """

      half_stamp_size = self.helper.get_stamp_center(stamp_size) + 1
      (rel_xc_values, rel_yc_values) = zip(*se_info_dico["rel_centroids"].values()) 
      rel_xc_values = numpy.asarray(rel_xc_values)
      rel_yc_values = numpy.asarray(rel_yc_values)

      max_shift = max( numpy.max(numpy.absolute(rel_xc_values - half_stamp_size)), 
                       numpy.max(numpy.absolute(rel_yc_values - half_stamp_size)) )

      return 2 * (half_stamp_size - max_shift)


   # -----------------------------------------------------------------------------------------------
   def _extract_galaxy_stamp(self, request, 
                                   image_data, (x, y), (rel_xc, rel_yc), half_galaxy_stamp_size, 
                                   job, worker):
      """!
         Extract a galaxy stamp of size galaxy_effective_stamp_size around the relative 
         centroid (rel_xc,rel_yc) inside the full-size galaxy stamp. The actual size of the
         extracted stamp is possibly reduced in order to remain square.   
      """

      # --- Compute the largest stamp that can be extracted, given the centroid shifts in x and y
      max_shift = max(math.fabs(rel_xc - half_galaxy_stamp_size), 
                      math.fabs(rel_yc - half_galaxy_stamp_size))
      max_allowed_stamp_size = 2 * (half_galaxy_stamp_size - max_shift)

      # --- Extract the postage stamp of the requested size around the centroid
      try:
         galaxy_actual_stamp_size = int(min(
                                      request.galaxy_effective_stamp_size, max_allowed_stamp_size))
         full_galaxy_stamp = image_data[y:y+request.galaxy_stamp_size, 
                                        x:x+request.galaxy_stamp_size].astype(
                                                                   request.galaxy_pixel_float_size)

         galaxy_stamp = self.helper.cut_stamp_around_centroid(full_galaxy_stamp, 
                                                              galaxy_actual_stamp_size,  
                                                              (rel_yc, rel_xc)) 

         # --- Check the validity of the cut-out stamp 
         if galaxy_actual_stamp_size < request.galaxy_effective_stamp_size:
            if worker.logging_enabled():
               worker.logger.log_warning_p(           
                  "{0} - Reduced the galaxy postage size from {1} to {2} in: "\
                  "/{3}/image-{4:03d}-{5:1d} - {6} - [Centroid: {7}]".format(
                      worker.name, (request.galaxy_effective_stamp_size, 
                                    request.galaxy_effective_stamp_size), 
                                    galaxy_stamp.shape, job.get_branch_tree(), job.img_no, job.epoch, 
                      request.image_file_type, (rel_xc, rel_yc)))

         if galaxy_stamp.shape != (galaxy_actual_stamp_size, galaxy_actual_stamp_size):

            # ... Probable centroid error, try to extract from the geometrical stamp center
            centroid = (galaxy_half_size, galaxy_half_size) 
            [galaxy_stamp, galaxy_actual_stamp_size, full_galaxy_stamp] =\
                    self._extract_galaxy_stamp(request, 
                                        galaxy_image_data, (x,y), centroid, 
                                        galaxy_half_size, job, worker)

            (Xc, Yc) = request.galaxy_info_dico["GAL_Xc_Yc"][(x,y)]
            raise Exception("Invalid galaxy shape at "\
                            "(x,y)=({0:.2f}, {1:.2f}) <=> (Xc, Yc)=({2:.2f}, {3:.2f})".format(
                                                                                   x, y, Xc, Yc))

         return [galaxy_stamp, galaxy_actual_stamp_size, full_galaxy_stamp]

      except Exception as detail:

         if worker.logging_enabled():
            worker.logger.log_warning_p(           
               "{0} - Error while extracting galaxy postage in: "\
               "/{1}/image-{2:03d}-{3:1d} - {4} - [Centroids: {5}] - {6}".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                   request.image_file_type, (rel_xc, rel_yc), detail))

         return [None, 0, None]

   # -----------------------------------------------------------------------------------------------
   def _extract_psf_stamp(self, request, 
                                image_data, (x, y), (rel_xc, rel_yc), half_psf_stamp_size, 
                                job, worker):
      """!
         Extract a PSF stamp of size psf_effective_stamp_size around the relative 
         centroid (rel_xc,rel_yc) inside the full-size psf stamp. The actual size of the
         extracted stamp is possibly reduced in order to remain square.   
      """

      # --- Compute the largest stamp that can be extracted, given the centroid shifts in x and y
      max_shift = max(math.fabs(rel_xc - half_psf_stamp_size), 
                      math.fabs(rel_yc - half_psf_stamp_size))
      max_allowed_stamp_size = 2 * (half_psf_stamp_size - max_shift)

      try:
         # --- Extract the postage stamp of the requested size around the centroid
         psf_actual_stamp_size = int(min(
                                        request.psf_effective_stamp_size, max_allowed_stamp_size))
         full_psf_stamp = image_data[y:y+request.psf_stamp_size, 
                                     x:x+request.psf_stamp_size].astype(
                                                                   request.psf_pixel_float_size)

         psf_stamp = self.helper.cut_stamp_around_centroid(full_psf_stamp, 
                                                           psf_actual_stamp_size,  
                                                           (rel_yc, rel_xc)) 

         # --- Check the validity of the cut-out stamp 
         if psf_actual_stamp_size < request.psf_effective_stamp_size:
            if worker.logging_enabled():
               worker.logger.log_warning_p(           
                  "{0} - Reduced the PSF postage size from {1} to {2} in: "\
                  "/{3}/image-{4:03d}-{5:1d} - {6} - [Centroid: {7}]".format(
                      worker.name, (request.psf_effective_stamp_size, 
                                    request.psf_effective_stamp_size), 
                                    psf_stamp.shape, job.get_branch_tree(), job.img_no, job.epoch, 
                                   request.image_file_type, (rel_xc, rel_yc)))

         if psf_stamp.shape != (psf_actual_stamp_size, psf_actual_stamp_size):

            # ... Probable centroid error, try to extract from the geometrical stamp center
            centroid = (psf_half_size, psf_half_size) 
            [psf_stamp, psf_actual_stamp_size, full_psf_stamp] =\
                    self._extract_psf_stamp(request, 
                                        psf_image_data, (x,y), centroid, 
                                        psf_half_size, job, worker)

            (Xc, Yc) = request.galaxy_info_dico["GAL_Xc_Yc"][(x,y)]
            raise Exception("Invalid PSF shape at "\
                            "(x,y)=({0:.2f}, {1:.2f}) <=> (Xc, Yc)=({2:.2f}, {3:.2f})".format(
                                                                                   x, y, Xc, Yc))

         return [psf_stamp, psf_actual_stamp_size, full_psf_stamp]

      except Exception as detail:

         if worker.logging_enabled():
            worker.logger.log_warning_p(           
               "{0} - Error while extracting galaxy postage in: "\
               "/{1}/image-{2:03d}-{3:1d} - {4} - [Centroids: {5}] - {6}".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                   request.image_file_type, (rel_xc, rel_yc), detail))

         return [None, 0, None]


   # -----------------------------------------------------------------------------------------------
   def _write_stamp_for_debugging(self, stamp, label, igal, (xc, yc), output_dir, job, worker):

      stamp_filename = "{0}_{1:03d}-{2:1d}_{3:04d}_{4:.2f}_{5:.2f}.fits".format(label,
                                                                         job.img_no, job.epoch, 
                                                                         igal, xc, yc)
      self.helper.write_as_fits(stamp, os.path.join(output_dir, stamp_filename))


   # -----------------------------------------------------------------------------------------------
   def _setup_stamp_mosaic(self, galaxy_image_data, 
                                 default_stamp_size, effective_stamp_size, max_common_stamp_size, 
                                 pixel_size, job, worker):

      input_mosaic_stamp_size = worker.config.get_as_int("MOSAIC_STAMP_SIZE", "DEBUGGING")
      mosaic_stamp_size = int(min(min(input_mosaic_stamp_size, max_common_stamp_size), 
                              effective_stamp_size))

      if input_mosaic_stamp_size > min(max_common_stamp_size, effective_stamp_size):
         if worker.logging_enabled():
            worker.logger.log_warning_p(           
               "{0} - /{1}/image-{2:03d}-{3:1d} - "\
               "Maximum possible mosaic galaxy postage size: {4} -> {5}".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch,
                   input_mosaic_stamp_size, mosaic_stamp_size))

      (mosaic_height, mosaic_width) = galaxy_image_data.shape
      nb_image_rows = mosaic_height / default_stamp_size
      nb_image_cols = mosaic_width / default_stamp_size
      
      mosaic_height = nb_image_rows * mosaic_stamp_size
      mosaic_width  = nb_image_cols * mosaic_stamp_size

      mosaic_data = numpy.zeros((mosaic_height, mosaic_width)).astype(pixel_size)

      return mosaic_data,  mosaic_stamp_size


   # -----------------------------------------------------------------------------------------------
   def _plot_stamp_for_debugging(self, stamp, label, igal, (xc, yc), output_dir, is_3D, 
                                       job, worker):

      plot_title = "{0} {1:03d}-{2:1d}-{3:04d} - ({4:.2f},{5:.2f})".format(label,
                                                                         job.img_no, job.epoch, 
                                                                         igal, xc, yc)
      if is_3D:
         plot_filename = "{0}_{1:03d}-{2:1d}_{3:04d}_{4:.2f}_{5:.2f}_3D.png".format(label,
                                                                            job.img_no, job.epoch, 
                                                                            igal, xc, yc)
      else:
         plot_filename = "{0}_{1:03d}-{2:1d}_{3:04d}_{4:.2f}_{5:.2f}.png".format(label,
                                                                            job.img_no, job.epoch, 
                                                                            igal, xc, yc)

      if is_3D:
         self.plotter.plot_stamp_3D(stamp, plot_title,
                                           output_dir=output_dir, output_file=plot_filename)
      else:
         self.plotter.plot_stamp(stamp, plot_title, color_bar=True,
                                        output_dir=output_dir, output_file=plot_filename)


   # -----------------------------------------------------------------------------------------------
   def _setup_galaxy_fitter(self, request, job, worker): 
      """!
         Create and configure the GalaxyFitter object
         @param request shape measurement request 
      """

      (image_filename, inage_fileext) = os.path.splitext(request.image_file_type)
      model_config_filename  = "{0}_{1}-{2:03d}-{3:1d}.cfg".format(request.galaxy_model_name,
                                                            image_filename, 
                                                            job.img_no, job.epoch)
      method_config_filename = "{0}_{1}-{2:03d}-{3:1d}.cfg".format(
                                                            request.galaxy_fitting_method_name, 
                                                            image_filename, 
                                                            job.img_no, job.epoch)

      # --- Create multifit ModelFitter
      fitter = ModelFitter(model_name=request.galaxy_model_name, 
                           method_name=request.galaxy_fitting_method_name, 
                           model_config_filename=model_config_filename, 
                           method_config_filename=method_config_filename,
                           model_config_dir=request.mfit_model_config_path, 
                           method_config_dir=request.mfit_method_config_path, 
                           fitting_config_dir=request.mfit_fitting_config_path, 
                           model_module_dir=request.mfit_model_module_path,
                           method_module_dir=request.mfit_method_module_path)

      # --- Use the standard Chi2 template objective function
      fitter.method.objective_func = fitter.method.get_chi2_objective_func() 

      # --- Dump the configuration files of the model and method      
      if worker.logging_enabled():
         
         self.helper.dump_config(fitter.model.config, worker.logger, "MODEL ")
         self.helper.dump_config(fitter.method.config, worker.logger, "METHOD ")
         if fitter.method.fitting_config is not None:
            self.helper.dump_config(fitter.method.fitting_config, worker.logger, "MODEL FITTING ")

      # --- Also make a copy of the configuration files to the log directory as a record
      target_model_config_dir, target_model_config_name = os.path.split(fitter.model.config.config_path) 
      target_model_config_path = os.path.join(worker.log_output_dir,
                                              target_model_config_name)

      if not worker.helper.file_exists(target_model_config_path):
         shutil.copy(fitter.model.config.config_path, worker.log_output_dir)  

      target_method_config_dir, target_method_config_name = os.path.split(fitter.method.config.config_path) 
      target_method_config_path = os.path.join(worker.log_output_dir, 
                                               target_method_config_name)

      if not worker.helper.file_exists(target_method_config_path):
         shutil.copy(fitter.method.config.config_path, worker.log_output_dir)  

      if not fitter.method.fitting_config is None:  
         target_fitting_config_dir, target_fitting_config_name = os.path.split(fitter.method.fitting_config.config_path) 
         target_fitting_config_path = os.path.join(worker.log_output_dir, 
                                                   target_fitting_config_name)
         if not worker.helper.file_exists(target_fitting_config_path):
            shutil.copy(fitter.method.fitting_config.config_path, worker.log_output_dir)  

      return fitter


#   # -----------------------------------------------------------------------------------------------
#   def _setup_galaxy_filter(self, request, job, worker):
#      """! Instantiate and initialize a ObjectFilter object """

#      return GalaxyFilter(request)

#   # -----------------------------------------------------------------------------------------------
#   def _setup_psf_filter(self, request, job, worker):
#      """! Instantiate and initialize a PSFFilter object """

#      pass

   
   # -----------------------------------------------------------------------------------------------
   def _setup_diag_options(self, fitter, igal, job, worker):   

      # Diagnostic options
      diag_log_file_name = "{0}_{1:03d}-{2:1d}_{3:04d}.log".format(
                                     fitter.method.name, job.img_no, job.epoch, igal)

      diag_options = fitter.method.diag_options
      if diag_options.diag_flags > 0:

         # --- Diagmostic options: directories and filenames
         diag_options.log_file_name = diag_log_file_name   
         fitter.method.diag_options = diag_options

         diag_log_directory = os.path.join(worker.log_output_dir, "multifit")
         if not os.path.isdir(diag_log_directory):
            self.helper.make_dir(diag_log_directory)

         #diag_options = fitter.method.diag_options
         diag_options.log_output_dir = diag_log_directory  

         diag_plot_directory = os.path.join(worker.plot_output_dir, job.get_branch_tree())
         diag_plot_directory = os.path.join(diag_plot_directory, "multifit")
         if not os.path.isdir(diag_plot_directory):
            self.helper.make_dir(diag_plot_directory)      
         diag_plot_directory = os.path.join(diag_plot_directory, "{0}_{1:03d}-{2:1d}_{3:04d}".format(
                                        fitter.method.name, job.img_no, job.epoch, igal))
         if not os.path.isdir(diag_plot_directory):
            self.helper.make_dir(diag_plot_directory)
         diag_options.plot_output_dir = diag_plot_directory

         fitter.method.diag_options = diag_options 


   # -----------------------------------------------------------------------------------------------
   def _create_galaxy_model(self, params, model, (x,y), galaxy_stamp, psf_stamp, request, 
                                  job, worker):
      """!
         Create a postage stamp representing a PSF-convolved galaxy from a given model.
         This method is invoked by the objective function of the underlying minimizer.   
      """   
      try:

         # --- Create a galaxy of the requested model
         if galaxy_stamp.shape[0] != 0:
            sersic_stamp = model.func(model, 
                                      galaxy_stamp.shape[0],       
                                      request.galaxy_subscale, *params)  

            if sersic_stamp is None:
               convolved_sersic_stamp = numpy.nan
               #raise Exception("Galaxy model coud not be created")
            else:
               # --- Convolve with the input PSF and returned the convolved sersic model profile
               # TODO: try shift in Fourier space instead

               # Adjust the PSF size to that of the galaxy and perform the cnvolution
               convolved_sersic_stamp = scipy.signal.fftconvolve(psf_stamp, 
                                                                 sersic_stamp,  mode='same')

               ########## TEMP #########               
#               guesses = model.get_best_guess_params( (x,y), galaxy_stamp, 
#                                                                request.se_galaxy_info_dico)
#               flux = model.get_param_by_name("flux").guess_value
#               flux = params[0]
##               flux = numpy.sum(galaxy_stamp)
#               print flux, params
#               convolved_sersic_stamp *= flux
               ########## TEMP #########               


               # --- The final galaxy stamp size may be larger than the original size because of the
               #     convolution with the PSF => Cut out the stamp size to match the original
               #     galaxy_actual_stamp_size (i.e. shape[0])   
               half_size = galaxy_stamp.shape[0] / 2.0
               convolved_sersic_stamp = self.helper.cut_stamp_around_centroid(
                                                                          convolved_sersic_stamp, 
                                                                          galaxy_stamp.shape[0], 
                                                                          (half_size, half_size))         
#            if __debug__:
#               output_dir =  os.path.join(worker.result_output_dir, job.get_branch_tree())
#               write_time = time.clock()
#               stamp_filename = "convolved_galaxy_{0:03d}-{1:1d}_{2}.fits".format(
#                                                                 job.img_no, job.epoch, write_time)
#               self.helper.write_as_fits(convolved_sersic_stamp, os.path.join(output_dir, stamp_filename))

#               #raise ZeroDivisionError("Numerical error during fitting")
         else:
            # --- malformed stamp
            convolved_sersic_stamp = numpy.nan

      except Exception as detail:

         # --- Exception during fitting
         if worker.logging_enabled():
            worker.logger.log_error_p(           
               "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - Fitting failure: {5}".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                                request.image_file_type, detail))

         convolved_sersic_stamp = numpy.nan

      return convolved_sersic_stamp


   # -----------------------------------------------------------------------------------------------
   def _find_unprocessed_ids(self, processed_ids, all_ids):

      return list(set(all_ids).difference(set(processed_ids)))

   # -----------------------------------------------------------------------------------------------
   def _find_unprocessed_values(self, processed_values, all_values):

      return list(set(all_values).difference(set(processed_values)))



# -------------------------------------------------------------------------------------------------
class ShapeMeasurementRequest(object):
   """! 
      Represent a request for shape measurement of the galaxies contained in an image made of
      galaxy postage stamps aranged in a grid   
   """
   def __init__(self, image_file_type, job, worker):
      """! 
         Construct a ShapeMeasurementRequest object 
         @param image_file_type galaxy image fie type (like image.fits or deep-image.fits)
         @param job GfitJob object
         @param worker Worker process object
      """         
      self._helper  = GfitHelper()                 # helper utility functions

      self._image_file_type = image_file_type      # should be a galaxy inage

      # --- Absolute file paths for the required images and catalogs
      self._galaxy_image_filepath = job.get_file_path(image_file_type)          
      self._psf_image_filepath = self._get_psf_image_filepath(image_file_type, job, worker)          
      self._galaxy_catalog_filepath = self._get_galaxy_catalog_filepath(job, worker)
      if self._galaxy_catalog_filepath is None:
         if worker.logging_enabled():
            msg = "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - "\
                   "Could not find the galaxy catalog for that image ...".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, image_file_type)
            worker.logger.log_error_p(msg)
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      self._galaxy_se_catalog_filepath = self._get_galaxy_se_catalog_filepath(image_file_type, 
                                                                              job, worker)  
      self._psf_se_catalog_filepath = self._get_psf_se_catalog_filepath(image_file_type, 
                                                                        job, worker)  

      # --- Galaxy image properties for shape measurement
      self._galaxy_image_hdu_no = worker.config.get_as_int("HDU_NO", 
                                                      "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")    
      self._galaxy_catalog_hdu_no = worker.config.get_as_int("HDU_NO", 
                                                      "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY")    

      self._galaxy_pixel_float_size = worker.config.get_as_string("PIXEL_FLOAT_SIZE", 
                                                      "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
      if  worker.config.has_key("IS_SEXTRACTOR_FORMAT", "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY"):  
         self._galaxy_catalog_is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                      "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY") 
      else:
         # Assumes SExtractor
         if worker.logging_enabled():
            worker.logger.log_warning_p(
                   "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - "\
                   "Assuming galaxy catalog has SExtractor format ...".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
         self._galaxy_catalog_is_sextractor = True

      self._galaxy_stamp_size = worker.config.get_as_int("ORIGINAL_STAMP_SIZE", 
                                                         "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
      self._galaxy_effective_stamp_size = worker.config.get_as_int("EFFECTIVE_STAMP_SIZE", 
                                                         "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")

      # --- Load relevant galaxy SExtractor catalog data per coordinates

      if  worker.config.has_key("IS_SEXTRACTOR_FORMAT", "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES"):  
         self._se_galaxy_catalog_is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                    "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES") 
      else:
         if worker.logging_enabled():
            worker.logger.log_warning_p(
                   "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - "\
                   "Assuming Sextractor galaxy catalog has indeed SExtractor format ...".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
         self._se_galaxy_catalog_is_sextractor = True

      self._se_galaxy_catalog_hdu_no = worker.config.get_as_int("HDU_NO", 
                                                    "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES") 

      self._se_galaxy_info_dico = self._get_se_galaxy_info_dico(job, worker)

      # --- Load relevant galaxy catalog data
      self._galaxy_info_dico = self._get_galaxy_info_dico(job, worker)

      # --- Load galaxy IDs per coordinates
      self._galaxy_id_dico = self._galaxy_info_dico["GAL_id"] 

#      # --- Load galaxy tile index per coordinate
#      self._galaxy_tile_index_dico = self._get_galaxy_tile_index_dico(self._galaxy_catalog_filepath, 
#                                                             self._galaxy_stamp_size, job, worker)

#      # --- Load galaxy tile position per coordinate
#      self._galaxy_tile_pos_dico = self._get_galaxy_tile_pos_dico(self._galaxy_catalog_filepath, 
#                                                             self._galaxy_stamp_size, job, worker)

#      # --- Load galaxy field position per coordinate
#      self._galaxy_field_pos_dico = self._get_galaxy_field_pos_dico(self._galaxy_catalog_filepath, 
#                                                             self._galaxy_stamp_size, 
#                                                             job, worker)

#      self._galaxy_field_angle_dico = self._get_galaxy_field_angle_dico(
#                                                                   self._galaxy_catalog_filepath, 
#                                                                   self._galaxy_stamp_size, 
#                                                                   job, worker)

      # --- PSF image properties for shape measurement
      self._psf_image_hdu_no = worker.config.get_as_int("HDU_NO", "PSF_DATASET.IMAGE_PROPERTIES")    
      self._psf_catalog_hdu_no = worker.config.get_as_int("HDU_NO", "PSF_DATASET.IMAGE_PROPERTIES")    
      self._psf_pixel_float_size = worker.config.get_as_string("PIXEL_FLOAT_SIZE", 
                                                               "PSF_DATASET.IMAGE_PROPERTIES")

      if  worker.config.has_key("IS_SEXTRACTOR_FORMAT", "PSF_DATASET.IMAGE_PROPERTIES"):  
         self._psf_catalog_is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                                   "PSF_DATASET.IMAGE_PROPERTIES") 
      else:
         # Assumes SExtractor
         self._psf_catalog_is_sextractor = True

      if  worker.config.has_key("IS_SEXTRACTOR_FORMAT", "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES"):  
         self._se_psf_catalog_is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                      "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES") 
      else:
         if worker.logging_enabled():
            worker.logger.log_warning_p(
                   "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - "\
                   "Assuming Sextractor star catalog has indeed SExtractor format ...".format(
                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))
         self._se_psf_catalog_is_sextractor = True

      self._se_psf_catalog_hdu_no = worker.config.get_as_int("HDU_NO", 
                                                     "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES") 

      self._psf_stamp_size = worker.config.get_as_int("ORIGINAL_STAMP_SIZE", 
                                                      "PSF_DATASET.IMAGE_PROPERTIES")
      self._psf_effective_stamp_size = worker.config.get_as_int(
                                                      "EFFECTIVE_STAMP_SIZE", 
                                                      "PSF_DATASET.IMAGE_PROPERTIES")

      # --- Load relevant PSF SExtractor catalog data per coordinates
      self._se_psf_info_dico = self._get_se_psf_info_dico(job, worker)


      # --- Galaxy Fitting
      self._galaxy_model_name = worker.config.get_as_string("GALAXY_MODEL_NAME", "GALAXY_FITTING")
      self._galaxy_fitting_method_name = worker.config.get_as_string("GALAXY_FITTING_METHOD", 
                                                                     "GALAXY_FITTING")
      self._galaxy_subscale = worker.config.get_as_int("GALAXY_SUBSCALE", "SHAPE_MEASUREMENT")

      self._normalize_galaxy = worker.config.get_as_boolean("NORMALIZE_GALAXY", "SHAPE_MEASUREMENT")

      self._normalize_psf = worker.config.get_as_boolean("NORMALIZE_PSF", "SHAPE_MEASUREMENT")

      self._substract_psf_sky = worker.config.get_as_boolean("SUBSTRACT_PSF_SKY", 
                                                            "SHAPE_MEASUREMENT")
      self._substract_galaxy_sky = worker.config.get_as_boolean("SUBSTRACT_GALAXY_SKY", 
                                                               "SHAPE_MEASUREMENT")
      self._failed_ellipticity_value = worker.config.get_as_float("FAILED_ELLIPTICITY_VALUE",  
                                                                  "SHAPE_MEASUREMENT")
      
      # --- Object filtering
      self._enable_galaxy_filtering = worker.config.get_as_boolean("ENABLE_FILTERING", 
                                                                   "OBJECT_FILTERING.GALAXY")
      self._enable_psf_filtering = worker.config.get_as_boolean("ENABLE_FILTERING", 
                                                                "OBJECT_FILTERING.PSF")
   
      self._galaxy_filtering_policy = worker.config.get_as_string("FILTERING_POLICY", 
                                                                  "OBJECT_FILTERING.GALAXY")
      self._psf_filtering_policy = worker.config.get_as_string("FILTERING_POLICY", 
                                                               "OBJECT_FILTERING.PSF")

      self._galaxy_filter_dico = worker.config.get_as_dict("FILTER", "OBJECT_FILTERING.GALAXY")
      self._psf_filter_dico = worker.config.get_as_dict("FILTER", "OBJECT_FILTERING.PSF")

      #self._filtered_ellipticity_value = worker.config.get_as_float("FILTERED_ELLIPTICITY_VALUE",  
      #                                                              "SHAPE_MEASUREMENT")

      # --- Multifit config
      self._mfit_model_config_path  = worker.config.get_as_string("MODEL_CONFIG_PATH", "MULTIFIT")   
      self._mfit_method_config_path = worker.config.get_as_string("METHOD_CONFIG_PATH", "MULTIFIT")   
      self._mfit_fitting_config_path = worker.config.get_as_string("FITTING_CONFIG_PATH", "MULTIFIT")   
      self._mfit_model_module_path = worker.config.get_as_string("MODEL_MODULE_PATH", "MULTIFIT")   
      self._mfit_method_module_path =  worker.config.get_as_string("METHOD_MODULE_PATH", "MULTIFIT")   

      # --- Regularized Least-squares      
      #self._regularization_lambda = worker.config.get_as_float("LEASTSQ_REG_LAMBDA", 
      #                                                         "GALAXY_FITTING")


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def helper(self):
      """! @return the GfitHelper instance. """
      return self._helper

   @property
   def image_file_type(self):
      """! image_file_type
         Get the file type of the galaxy image
         @return galaxy image file type """
      return self._image_file_type

#   @property
#   def catalog_file_extension(self):
#      """! 
#         Get the preferred catalog file extension (like .txt or .fits)
#         @return the preferred catalog file extension
#      """
#      return self._catalog_file_extension
      
   @property
   def galaxy_image_filepath(self):
      """! 
         Get the full path of the galaxy image
         @return galaxy image full path
      """
      return self._galaxy_image_filepath

   @property
   def psf_image_filepath(self):
      """! 
         Get the full path of the PSF image
         @return PSF image full path
      """
      return self._psf_image_filepath

   @property
   def galaxy_catalog_filepath(self):
      """! 
         Get the full path of the galaxy gatalog
         @return galaxy gatalog full path 
      """
      return self._galaxy_catalog_filepath

   @property
   def galaxy_se_catalog_filepath(self):
      """! 
         Get the full path of the SExtractor galaxy gatalog
         @return SExtractor galaxy gatalog full path 
      """
      return self._galaxy_se_catalog_filepath

   @property
   def psf_se_catalog_filepath(self):
      """! 
         Get the full path of the SExtractor PSF gatalog
         @return SExtractor PSF gatalog full path 
      """
      return self._psf_se_catalog_filepath

   @property
   def galaxy_info_dico(self):
      """! 
         Get the dictionary containing relevant galaxy catalog data
         @return dictionary containing  relevant galaxy catalog data
      """
      return self._galaxy_info_dico

   @property
   def se_galaxy_info_dico(self):
      """! 
         Get the dictionary containing SExtractor data on galaxies
         @return dictionary containing SExtractor data on galaxies
      """
      return self._se_galaxy_info_dico

   @property
   def se_psf_info_dico(self):
      """! 
         Get the dictionary containing SExtractor data on PSF kernels
         @return dictionary containing SExtractor data on PSF kernels
      """
      return self._se_psf_info_dico

   @property
   def galaxy_id_dico(self):
      """! 
         Get the dictionary containing galaxy IDs per coordinates
         @return dictionary containing galaxy IDs per coordinates
      """
      return self._galaxy_id_dico

#   @property
#   def galaxy_tile_index_dico(self):
#      """! 
#         Get the dictionary containing galaxy tile indice per coordinates
#         @return dictionary containing galaxy tile indice per coordinates
#      """
#      return self._galaxy_tile_index_dico

#   @property
#   def galaxy_tile_pos_dico(self):
#      """! 
#         Get the dictionary containing galaxy tile positions per coordinates
#         @return dictionary containing galaxy tile positions per coordinates
#      """
#      return self._galaxy_tile_pos_dico

#   @property
#   def galaxy_field_pos_dico(self):
#      """! 
#         Get the dictionary containing galaxy field positions per coordinates
#         @return dictionary containing galaxy field positions per coordinates
#      """
#      return self._galaxy_field_pos_dico

   @property
   def galaxy_field_angle_dico(self):
      """! 
         Get the dictionary containing galaxy field positions (alpha) per coordinates
         @return dictionary containing galaxy field positions (delta) per coordinates
      """
      return self._galaxy_field_angle_dico

   @property
   def galaxy_stamp_size(self):
      """! 
         Get the size of a full galaxy postage stamp
         @return size of a full galaxy postage stamp
      """
      return self._galaxy_stamp_size


   @property
   def psf_stamp_size(self):
      """! 
         Get the size of a full PSF postage stamp
         @return size of a full PF postage stamp
      """
      return self._psf_stamp_size

   @property
   def galaxy_effective_stamp_size(self):
      """! 
         Get the size of the effective galaxy stamp size to use for shape measurement
         @return size of the efffectve galaxy stamp size to use for shape measurement
      """
      return self._galaxy_effective_stamp_size


   @property
   def psf_effective_stamp_size(self):
      """! 
         Get the size of the effective PSF stamp size to use for shape measurement
         @return size of the effective PSF stamp size to use for shape measurement
      """
      return self._psf_effective_stamp_size

   @property
   def psf_actual_stamp_size(self):
      """! 
         Get the size of the actual PSF stamp size to use for shape measurement
         @return size of the actual PSF stamp size to use for shape measurement
      """
      return self._psf_actual_stamp_size

   @property
   def galaxy_image_hdu_no(self):
      """! 
         Get the the FITS HDU number for a galaxy image
         @return FITS HDU number for a galaxy image
      """
      return self._galaxy_image_hdu_no

   @property
   def galaxy_catalog_hdu_no(self):
      """! 
         Get the the FITS HDU number for a galaxy catalog
         @return FITS HDU number for a galaxy catalog
      """
      return self._galaxy_catalog_hdu_no

   @property
   def psf_image_hdu_no(self):
      """! 
         Get the the FITS HDU number for a PSF image
         @return FITS HDU number for a PSF image
      """
      return self._psf_image_hdu_no

   @property
   def psf_catalog_hdu_no(self):
      """! 
         Get the the FITS HDU number for a PSF catalog
         @return FITS HDU number for a PSF catalog
      """
      return self._psf_catalog_hdu_no

   @property
   def se_galaxy_catalog_hdu_no(self):
      """! 
         Get the the FITS HDU number for a galaxy SExtractor FITS catalog
         @return FITS HDU number for a galaxy SExtractor FITScatalog
      """
      return self._se_galaxy_catalog_hdu_no

   @property
   def se_psf_catalog_hdu_no(self):
      """! 
         Get the the FITS HDU number for a PSF SExtractor FITS catalog
         @return FITS HDU number for a PSF SExtractor FITScatalog
      """
      return self._se_psf_catalog_hdu_no

   @property
   def galaxy_pixel_float_size(self):
      """! 
         Get the float size used to represent galaxy image pixel values (like "float32", "float64")
         @return float size used to represent galaxy image pixel values
      """
      return self._galaxy_pixel_float_size

   @property
   def psf_pixel_float_size(self):
      """! 
         Get the float size used to represent PSF image pixel values (like "float32", "float64")
         @return float size used to represent PSF image pixel values
      """
      return self._psf_pixel_float_size

   @property
   def galaxy_catalog_is_sextractor(self):
      """! 
         Tell whether the input galaxy catalog, if in text format, has SExtractor format
         @return True if the catalog has SExtractor format, False otherwise
      """
      return self._galaxy_catalog_is_sextractor

   @property
   def psf_catalog_is_sextractor(self):
      """! 
         Tell whether the input star catalog, if in text format, has SExtractor format
         @return True if the input star catalog has SExtractor format, False otherwise
      """
      return self._psf_catalog_is_sextractor

   @property
   def se_galaxy_catalog_is_sextractor(self):
      """! 
         Tell whether the input galaxy SExtractor catalog, if in text format, has SExtractor format
         @return True if the catalog has SExtractor format, False otherwise
      """
      return self._se_galaxy_catalog_is_sextractor

   @property
   def se_psf_catalog_is_sextractor(self):
      """! 
         Tell whether the input star SExtractor catalog, if in text format, has SExtractor format
         @return True if the catalog has SExtractor format, False otherwise
      """
      return self._se_psf_catalog_is_sextractor

   @property
   def normalize_galaxy(self):
      """! 
         Tell whether the galaxy image should bbe nornmalized or not for shape measurement
         @return True if the galaxy image should bbe nornmalized or not for shape measurement
      """
      return self._normalize_galaxy

   @property
   def normalize_psf(self):
      """! 
         Tell whether the PSF image should bbe nornmalized or not for shape measurement
         @return True if the PSF image should bbe nornmalized or not for shape measurement
      """
      return self._normalize_psf

   @property
   def substract_galaxy_sky(self):
      """! 
         Tell whether the galaxy sky noise must be removed or not for shape measurement
         @return True if the galaxy noise must be removed or not for shape measurement
      """
      return self._substract_galaxy_sky

   @property
   def failed_ellipticity_value(self):
      """! 
         Return the failed ellipticity marker value
         @return failed ellipticity marker value
      """
      return self._failed_ellipticity_value

#   @property
#   def filtered_ellipticity_value(self):
#      """! 
#         Return the filtered object ellipticity marker value
#         @return filtered object ellipticity marker value
#      """
#      return self._filtered_ellipticity_value

   @property
   def substract_psf_sky(self):
      """! 
         Tell whether the PSF sky noise must be removed or not for shape measurement
         @return True if the sky noise must be removed or not for shape measurement
      """
      return self._substract_psf_sky

   @property
   def enable_galaxy_filtering(self):
      """! 
         Tell whether galaxy filtering should be applied during shear measurement
         @return True if galaxy filtering should be applied during shear measurement
      """
      return self._enable_galaxy_filtering

   @property
   def enable_psf_filtering(self):
      """! 
         Tell whether PSF filtering should be applied during shear measurement
         @return True if PSF filtering should be applied during shear measurement
      """
      return self._enable_psf_filtering

   @property
   def galaxy_filtering_policy(self):
      """! 
         Return the galaxy filtering policy, i.e. what to do with filtered data (remove, etc.)
         @return galaxy filtering policy, i.e. what to do with filtered data
      """
      return self._galaxy_filtering_policy

   @property
   def psf_filtering_policy(self):
      """! 
         Return the PSF filtering policy, i.e. what to do with filtered data (remove, etc.)
         @return PSF filtering policy, i.e. what to do with filtered data
      """
      return self._psf_filtering_policy

   @property
   def galaxy_filter_dico(self):
      """! 
         Get the galaxy filter dico containing the filtering directives
         @return galaxy filter dico containing the filtering directives
      """
      return self._galaxy_filter_dico

   @property
   def psf_filter_dico(self):
      """! 
         Get the PSF filter dico containing the filtering directives
         @return PSF filter dico containing the filtering directives
      """
      return self._psf_filter_dico

   @property
   def galaxy_model_name(self):
      """! 
         Get the galaxy_model_name used to fit the Galaxy
         @return galaxy model name used to fit the Galaxy
      """
      return self._galaxy_model_name

   @property
   def galaxy_subscale(self):
      """! 
         Get the galaxy_subscale used to fit the Galaxy
         @return galaxy subscale used to fit the Galaxy
      """
      return self._galaxy_subscale

   @property
   def galaxy_fitting_method_name(self):
      """! 
         Get the galaxy fitting method name used to fit the Galaxy
         @return galaxy fitting method name used to fit the Galaxy
      """
      return self._galaxy_fitting_method_name

   @property
   def mfit_model_config_path(self):
      """! 
         Get the multifit config path for models
         @return multifit config path for models
      """
      return self._mfit_model_config_path

   @property
   def mfit_method_config_path(self):
      """! 
         Get the multifit config path for methods
         @return multifit config path for methods
      """
      return self._mfit_method_config_path

   @property
   def mfit_fitting_config_path(self):
      """! 
         Get the multifit config path for parameter fitting
         @return multifit config path for parameter fitting
      """
      return self._mfit_fitting_config_path

   @property
   def mfit_model_module_path(self):
      """! 
         Get the multifit module path for models
         @return multifit module path for models
      """
      return self._mfit_model_module_path
   
   @property
   def mfit_method_module_path(self):
      """! 
         Get the multifit module path for methods
         @return multifit module path for methods
      """
      return self._mfit_method_module_path

#   @property
#   def regularization_lambda(self):
#      """! 
#         Get the lambda regularization term to apply to the least.squares chi2 
#         @return lambda regularization term
#      """
#      return self._regularization_lambda


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_info_dico(self, job, worker):
      """! 
         Load data from the input galaxy catalog available in the primary dataset
         @param job GfitJob object
         @param worker Worker process object
      """

      info_dico = OrderedDict()

      # --- Read SE column mapping information
      col_mapping_dico = worker.config.get_as_dict("COL_MAPPING", 
                                                   "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY") 

      # --- Open SE galaxy catalog
      galaxy_catalog = self._open_catalog(self.galaxy_catalog_filepath, 
                                          self.galaxy_catalog_is_sextractor,
                                          self.galaxy_catalog_hdu_no, job, worker)

      # --- Available SExtractor column names
      cat_colnames = galaxy_catalog.get_col_names()   


      # --- Coordinates of galaxies, to be converted to coordimates of stamps lower-left corners
      xc_colname = col_mapping_dico.get("GAL_x", None)    # x coordinate in mosaic
      yc_colname = col_mapping_dico.get("GAL_y", None)    # y coordinate in mosaic
      Xc_colname = col_mapping_dico.get("GAL_X", None)    # Physical X centroid
      Yc_colname = col_mapping_dico.get("GAL_Y", None)    # Physical Y centroid

      if xc_colname is None or yc_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "GAL_x", "GAL_y")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)    

      if Xc_colname is None or Yc_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "GAL_X", "GAL_Y")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)    

      # --- Mapped coordinate columns are required
      if not xc_colname in cat_colnames or not yc_colname in cat_colnames:
         msg = "{0} - The {1} and {2} columns must be available in the galaxy catalog "\
               "for processing to continue".format(worker.name, xc_colname, yc_colname)     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)
         
      if not Xc_colname in cat_colnames or not Yc_colname in cat_colnames:
         msg = "{0} - The {1} and {2} columns must be available in the galaxy catalog "\
               "for processing to continue".format(worker.name, Xc_colname, Yc_colname)     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      xc_array = galaxy_catalog.get_named_col_data(xc_colname)
      yc_array = galaxy_catalog.get_named_col_data(yc_colname)

      Xc_array = galaxy_catalog.get_named_col_data(Xc_colname)
      Yc_array = galaxy_catalog.get_named_col_data(Yc_colname)

      # --- Force conversion to stamp lower-left corner coordinates
      x_array = numpy.floor(xc_array / self.galaxy_stamp_size) * self.galaxy_stamp_size
      y_array = numpy.floor(yc_array / self.galaxy_stamp_size) * self.galaxy_stamp_size

      # --- Record galaxy coordinates
      galaxy_x_y_coords   = zip(x_array, y_array)
      galaxy_xc_yc_coords = zip(xc_array, yc_array)
      galaxy_Xc_Yc_coords = zip(Xc_array, Yc_array)

      info_dico["GAL_xc_yc"] = dict(zip(galaxy_x_y_coords, galaxy_xc_yc_coords))
      info_dico["GAL_Xc_Yc"] = dict(zip(galaxy_x_y_coords, galaxy_Xc_Yc_coords))

      # --- Object ID (or sequence number)
      id_colname = col_mapping_dico.get("GAL_id", None)  # object id
      if id_colname is not None:
         if not id_colname in cat_colnames:
            msg = "{0} - The {1} column must be available in the galaxy catalog "\
                  "for processing to continue".format(worker.name, id_colname)     
            raise GfitShapeEstimator.ShapeMeasurementError(msg)         
         id_array = galaxy_catalog.get_named_col_data(id_colname)
      else:
         id_array = numpy.arange(0, len(galaxy_catalog.get_col_data(xc_colname)))
   
      info_dico["GAL_id"] = dict(zip(galaxy_x_y_coords, id_array))

      # --- Object RA and DEC if available
      RA_colname  = col_mapping_dico.get("GAL_RA",  None)    # RA  physical coordinates
      DEC_colname = col_mapping_dico.get("GAL_DEC", None)    # DEC physical coordinates
      if RA_colname is not None or DEC_colname is not None:
         if not RA_colname in cat_colnames or not RA_colname in cat_colnames:
            msg = "{0} - The {1} or the {2} column is not available in the galaxy catalog - "\
                  "Ignored".format(worker.name, RA_colname, DEC_colname)     
         else:
            RA_array  = galaxy_catalog.get_named_col_data(RA_colname)
            DEC_array = galaxy_catalog.get_named_col_data(DEC_colname)
            galaxy_RA_DEC_coords = zip(RA_array, DEC_array)
            info_dico["GAL_RA_DEC"] = dict(zip(galaxy_x_y_coords, galaxy_RA_DEC_coords))

      required_colnames = info_dico.keys()
      excluded_columns = ["GAL_x", "GAL_y", "GAL_X", "GAL_Y", "GAL_RA", "GAL_DEC"]
      excluded_columns.extend(required_colnames)

      # --- Keep track of extra SE columns if specified
      for colname in col_mapping_dico.keys():
         if not colname in excluded_columns:
            mapped_col_name = col_mapping_dico[colname]
            if mapped_col_name in cat_colnames:
               coldata = galaxy_catalog.get_named_col_data(mapped_col_name)
               info_dico[colname] = dict(zip(galaxy_x_y_coords, coldata))

      return info_dico


   # -----------------------------------------------------------------------------------------------
   def _get_se_galaxy_info_dico(self, job, worker):
      """! 
         Load SExtractor relevant information that may be useful for shape measurment
         @param job GfitJob object
         @param worker Worker process object
      """

      info_dico = OrderedDict()

      # --- Read SE column mapping information
      se_col_mapping_dico = worker.config.get_as_dict(
                                                "COL_MAPPING", 
                                                "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES") 

      # --- Open SE galaxy catalog
      se_galaxy_catalog = self._open_catalog(self.galaxy_se_catalog_filepath, 
                                             self.galaxy_catalog_is_sextractor,
                                             self.se_galaxy_catalog_hdu_no, job, worker)

      # --- Available SExtractor column names
      se_cat_colnames = se_galaxy_catalog.get_col_names()   

      # --- Coordinates of galaxy centroids
      x_image_colname = se_col_mapping_dico.get("SE_GAL_X_IMAGE", None)
      y_image_colname = se_col_mapping_dico.get("SE_GAL_Y_IMAGE", None)
      if x_image_colname is None or y_image_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "SE_GAL_X_IMAGE", "SE_GAL_Y_IMAGE")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)    

      # --- The SE X_IMAGE and Y_IMAGE mapped columns are required
      if not x_image_colname in se_cat_colnames or not y_image_colname in se_cat_colnames:
         msg = "{0} - The {1} and {2} columns must be available in the SExtractor catalog "\
               "for processing to continue".format(worker.name, x_image_colname, y_image_colname)     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      xc_array = se_galaxy_catalog.get_col_data(x_image_colname) - 1.0   # SE: start from 1, not 0
      yc_array = se_galaxy_catalog.get_col_data(y_image_colname) - 1.0   # SE: start from 1, not 0

      # --- Coordinates of stamp lower-left corners
      x_array = numpy.floor(xc_array / self.galaxy_stamp_size) * self.galaxy_stamp_size
      y_array = numpy.floor(yc_array / self.galaxy_stamp_size) * self.galaxy_stamp_size

      x_array = x_array.astype(int)
      y_array = y_array.astype(int)

      # --- Relative galaxy centroids (within stamp)
      rel_xc_array = numpy.floor(xc_array - x_array + 0.5)
      rel_yc_array = numpy.floor(yc_array - y_array + 0.5)

      galaxy_half_size = self.helper.get_stamp_center(self.galaxy_stamp_size) + 1
      rel_xc_array = numpy.where(galaxy_half_size - rel_xc_array <= 1.0, 
                                 numpy.ceil(galaxy_half_size), rel_xc_array)
      rel_yc_array = numpy.where(galaxy_half_size - rel_yc_array <= 1.0, 
                                 numpy.ceil(galaxy_half_size), rel_yc_array)

      # --- Record galaxy coordinates
      galaxy_coords = zip(x_array, y_array) 
      galaxy_centroids = zip(xc_array, yc_array)
      galaxy_rel_centroids = zip(rel_xc_array, rel_yc_array)

      info_dico["centroids"] = dict(zip(galaxy_coords, galaxy_centroids))
      info_dico["rel_centroids"] = dict(zip(galaxy_coords, galaxy_rel_centroids))

#      info_dico["SE_GAL_rel_xc"] = dict(zip(galaxy_coords, rel_xc_array)) 
#      info_dico["SE_GAL_rel_yc"] = dict(zip(galaxy_coords, rel_yc_array)) 

      # --- Extra information, useful for shape measurement
      flux_colname = se_col_mapping_dico.get("SE_GAL_FLUX", None)
      hlr_colname  = se_col_mapping_dico.get("SE_GAL_FLUX_RADIUS", None)
      if flux_colname is None or hlr_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "SE_GAL_FLUX", "SE_GAL_FLUX_RADIUS")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)    


      # --- The SE FLUX and SE_HLR mapped columns are required
      if not flux_colname in se_cat_colnames or not hlr_colname in se_cat_colnames:
         msg = "{0} - The {1}, {2} columns must be available in the SExtractor catalog "\
               "for processing to continue".format(worker.name,\
                 flux_colname, hlr_colname)          
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      se_flux_data = se_galaxy_catalog.get_named_col_data(flux_colname)
      se_hlr_data  = se_galaxy_catalog.get_named_col_data(hlr_colname)

      # --- Keep track of SE column values per coordinates
      info_dico["SE_GAL_FLUX"] = dict(zip(galaxy_coords, se_flux_data))
      info_dico["SE_GAL_FLUX_RADIUS"]  = dict(zip(galaxy_coords, se_hlr_data))

      required_colnames = info_dico.keys()
      required_colnames = info_dico.keys()
      excluded_columns = ["SE_GAL_X_IMAGE", "SE_GAL_Y_IMAGE"]
      excluded_columns.extend(required_colnames)

      # --- Compute SNR data if possible
      if "SE_GAL_FLUXERR" in se_col_mapping_dico: 
         flux_err_colname = se_col_mapping_dico["SE_GAL_FLUXERR"]
         if flux_err_colname in se_cat_colnames:
            se_flux_err_data = se_galaxy_catalog.get_named_col_data(flux_err_colname)
            se_snr_data = se_flux_data / se_flux_err_data
            info_dico["SE_GAL_SNR"] = dict(zip(galaxy_coords, se_snr_data))
            required_colnames.append("SE_GAL_SNR")
      
      # --- Keep track of extra SE columns if specified
      for colname in se_col_mapping_dico.keys():
         if not colname in excluded_columns:
            mapped_se_col_name = se_col_mapping_dico[colname]
            if mapped_se_col_name in se_cat_colnames:
               coldata = se_galaxy_catalog.get_named_col_data(mapped_se_col_name)
               info_dico[colname] = dict(zip(galaxy_coords, coldata))

      # --- Close Sextractor catalog
      se_galaxy_catalog.close()

      return info_dico


   # -----------------------------------------------------------------------------------------------
   def _get_se_psf_info_dico(self, job, worker):
      """! 
         Load SExtractor relevant information that may be useful for shape measurment
         @param psf_se_catalog_filepath path to the psf SExtractor catalog
      """
      info_dico = OrderedDict()


      # --- Read SE column mapping information
      se_col_mapping_dico = worker.config.get_as_dict("COL_MAPPING", 
                                                      "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES") 

      # --- Read ID data from PSF catalog
      se_psf_catalog = self._open_catalog(self.psf_se_catalog_filepath, 
                                          self.psf_catalog_is_sextractor,
                                          self.se_psf_catalog_hdu_no, job, worker)

      # --- Available SExtractor column names
      se_cat_colnames = se_psf_catalog.get_col_names()   

      # --- Coordinates of PSF centroids
      x_image_colname = se_col_mapping_dico.get("SE_PSF_X_IMAGE", None)
      y_image_colname = se_col_mapping_dico.get("SE_PSF_Y_IMAGE", None)
      if x_image_colname is None or y_image_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "SE_GAL_X_IMAGE", "SE_GAL_Y_IMAGE")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      # --- The SE X_IMAGE and Y_IMAGE mapped columns are required
      if not x_image_colname in se_cat_colnames or not y_image_colname in se_cat_colnames:
         msg = "{0} - The {1} and {2} columns must be available in the SExtractor catalog "\
               "for processing to continue".format(worker.name, x_image_colname, y_image_colname)
         if worker.logging_enabled():
            worker.logger.log_error_p(msg)           
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      xc_array = se_psf_catalog.get_col_data(x_image_colname) - 1.0   # SE: start from 1, not 0
      yc_array = se_psf_catalog.get_col_data(y_image_colname) - 1.0   # SE: start from 1, not 0

      # --- Coordinates of stamp lower-left corners
      x_array = numpy.floor(xc_array / self.psf_stamp_size) * self.psf_stamp_size
      y_array = numpy.floor(yc_array / self.psf_stamp_size) * self.psf_stamp_size

      x_array = x_array.astype(int)
      y_array = y_array.astype(int)

      # --- Relative PSF centroids (within stamp)
      rel_xc_array = numpy.floor(xc_array - x_array + 0.5)
      rel_yc_array = numpy.floor(yc_array - y_array + 0.5)

      # --- Record PSF coordinates
      psf_coords = zip(x_array, y_array) 
      psf_centroids = zip(xc_array, yc_array)
      psf_rel_centroids = zip(rel_xc_array, rel_yc_array)

      info_dico["psf_centroids"] = dict(zip(psf_coords, psf_centroids))
      info_dico["psf_rel_centroids"] = dict(zip(psf_coords, psf_rel_centroids))

      # --- Extra information, useful for shape measurement
      flux_colname = se_col_mapping_dico.get("SE_PSF_FLUX", None)
      hlr_colname  = se_col_mapping_dico.get("SE_PSF_FLUX_RADIUS", None)
      if flux_colname is None or hlr_colname is None:
         msg = "{0} - {1} and/or {2} missing from the 'COL_MAPPING' key".format(
                                                  worker.name, "SE_PSF_FLUX", "SE_PSF_FLUX_RADIUS")     
         raise GfitShapeEstimator.ShapeMeasurementError(msg)    

      # --- The SE FLUX and SE_HLR mapped columns are required
      if not flux_colname in se_cat_colnames or not hlr_colname in se_cat_colnames:
         msg = "{0} - The {1}, {2} columns must be available in the SExtractor catalog "\
               "for processing to continue".format(worker.name,\
                 flux_colname, hlr_colname)
         raise GfitShapeEstimator.ShapeMeasurementError(msg)         

      se_flux_data = se_psf_catalog.get_named_col_data(flux_colname)
      se_hlr_data = se_psf_catalog.get_named_col_data(hlr_colname)

      # --- Keep track of SE column values per coordinates
      info_dico["SE_PSF_FLUX"] = dict(zip(psf_coords, se_flux_data))
      info_dico["SE_PSF_FLUX_RADIUS"] = dict(zip(psf_coords, se_hlr_data))

      info_dico["SE_PSF_rel_xc"] = dict(zip(psf_coords, rel_xc_array)) 
      info_dico["SE_PSF_rel_yc"] = dict(zip(psf_coords, rel_yc_array)) 

      required_colnames = info_dico.keys()
      required_colnames = info_dico.keys()
      excluded_columns = ["SE_PSF_X_IMAGE", "SE_PSF_Y_IMAGE"]

      excluded_columns.extend(required_colnames)
      # --- Compute SNR data if possible
      if "SE_PSF_FLUXERR" in se_col_mapping_dico: 
         flux_err_colname = se_col_mapping_dico["SE_PSF_FLUXERR"]
         if flux_err_colname in se_cat_colnames:
            se_flux_err_data = se_psf_catalog.get_named_col_data(flux_err_colname)
            se_snr_data = se_flux_data / se_flux_err_data
            info_dico["SE_PSF_SNR"] = dict(zip(psf_coords, se_snr_data))
            required_colnames.append("SE_PSF_SNR")
      
      # --- Keep track of extra SE columns if specified
      for colname in se_col_mapping_dico.keys():
         if not colname in excluded_columns:
            mapped_se_col_name = se_col_mapping_dico[colname]
            if mapped_se_col_name in se_cat_colnames:
               coldata = se_psf_catalog.get_named_col_data(mapped_se_col_name)
               info_dico[colname] = dict(zip(psf_coords, coldata))

      # --- Close Sextractor catalog
      se_psf_catalog.close()

      return info_dico

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_tile_index_dico(self, galaxy_catalog_filepath, galaxy_stamp_size, job, worker):

      # --- Read tile index from galaxy catalog
      if job.dataset.is_fits(worker, galaxy_catalog_filepath):
         hdu_no = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
         galaxy_catalog = FITSCatalog(galaxy_catalog_filepath, hdu_no=hdu_no, mem_map=True)
      else:
         if self.galaxy_catalog_is_sextractor:
            galaxy_catalog = SExCatalog(galaxy_catalog_filepath)
         else:
            galaxy_catalog = TextCatalog(galaxy_catalog_filepath)

      # Create tile index dictionary per coordinates
      galaxy_catalog.open()

      tile_index_dico = None
      if galaxy_catalog.col_exists("x_tile_index") and galaxy_catalog.col_exists("y_tile_index"):

         half_stamp_size = self.helper.get_stamp_center(galaxy_stamp_size)

         x_tile_index = galaxy_catalog.get_col_data("x_tile_index") 
         y_tile_index = galaxy_catalog.get_col_data("y_tile_index") 
         x_coords = galaxy_catalog.get_col_data("x") - half_stamp_size
         y_coords = galaxy_catalog.get_col_data("y") - half_stamp_size 

         tile_index_dico = dict(zip(zip(x_coords, y_coords), zip(x_tile_index, y_tile_index)))

      galaxy_catalog.close()

      return tile_index_dico

#   # -----------------------------------------------------------------------------------------------
#   def _get_galaxy_tile_pos_dico(self, galaxy_catalog_filepath,
#                                       galaxy_stamp_size, job, worker):

#      # --- Read tile index from galaxy catalog
#      if job.dataset.is_fits(worker, galaxy_catalog_filepath):
#         hdu_no = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
#         galaxy_catalog = FITSCatalog(galaxy_catalog_filepath, hdu_no=hdu_no, mem_map=True)
#      else:
#         if self.galaxy_catalog_is_sextractor:
#            galaxy_catalog = SExCatalog(galaxy_catalog_filepath)
#         else:
#            galaxy_catalog = TextCatalog(galaxy_catalog_filepath)

#      tile_pos_dico = None
#      galaxy_catalog.open()

#      if galaxy_catalog.col_exists("x_tile_position") and \
#         galaxy_catalog.col_exists("y_tile_position"):

#         # Create tile index dictionary per coordinates
#         half_stamp_size = self.helper.get_stamp_center(galaxy_stamp_size)
#         tile_x_pos_deg = galaxy_catalog.get_col_data("x_tile_position") 
#         tile_y_pos_deg = galaxy_catalog.get_col_data("y_tile_position") 
#         x_coords = galaxy_catalog.get_col_data("x") - half_stamp_size
#         y_coords = galaxy_catalog.get_col_data("y") - half_stamp_size 

#         tile_pos_dico = dict(zip(zip(x_coords, y_coords), zip(tile_x_pos_deg, tile_y_pos_deg)))

#      galaxy_catalog.close()

#      return tile_pos_dico

#   # -----------------------------------------------------------------------------------------------
#   def _get_galaxy_field_pos_dico(self, galaxy_catalog_filepath, galaxy_stamp_size, job, worker):

#      # --- Read field index from galaxy catalog
#      if job.dataset.is_fits(worker, galaxy_catalog_filepath):
#         hdu_no = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
#         galaxy_catalog = FITSCatalog(galaxy_catalog_filepath, hdu_no=hdu_no, mem_map=True)
#      else:
#         if self.galaxy_catalog_is_sextractor:
#            galaxy_catalog = SExCatalog(galaxy_catalog_filepath)
#         else:
#            galaxy_catalog = TextCatalog(galaxy_catalog_filepath)

#      # Create field index dictionary per coordinates
#      galaxy_catalog.open()

#      field_pos_dico = None
#      if galaxy_catalog.col_exists("GAL_Xc") and \
#         galaxy_catalog.col_exists("GAL_Yc"):

#         half_stamp_size = self.helper.get_stamp_center(galaxy_stamp_size)
#         x_field_true_deg = galaxy_catalog.get_col_data("GAL_Xc") 
#         y_field_true_deg = galaxy_catalog.get_col_data("GAL_Yc") 
#         x_coords = galaxy_catalog.get_col_data("x") - half_stamp_size
#         y_coords = galaxy_catalog.get_col_data("y") - half_stamp_size

#         field_pos_dico = dict(zip(zip(x_coords, y_coords), zip(x_field_true_deg, y_field_true_deg)))

#      galaxy_catalog.close()

#      return field_pos_dico

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_field_angle_dico(self, galaxy_catalog_filepath, galaxy_stamp_size, job, worker):

      # --- Read field index from galaxy catalog
      if job.dataset.is_fits(worker, galaxy_catalog_filepath):
         hdu_no = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
         galaxy_catalog = FITSCatalog(galaxy_catalog_filepath, hdu_no=hdu_no, mem_map=True)
      else:
         if self.galaxy_catalog_is_sextractor:
            galaxy_catalog = SExCatalog(galaxy_catalog_filepath)
         else:
            galaxy_catalog = TextCatalog(galaxy_catalog_filepath)

      # Create field index dictionary per coordinates
      galaxy_catalog.open()

      angle_pos_dico = None
      if galaxy_catalog.col_exists("ALPHA") and \
         galaxy_catalog.col_exists("DELTA"):

         half_stamp_size = self.helper.get_stamp_center(galaxy_stamp_size)
         x_field_true_deg = galaxy_catalog.get_col_data("ALPHA") 
         y_field_true_deg = galaxy_catalog.get_col_data("DELTA") 
         x_coords = galaxy_catalog.get_col_data("x") - half_stamp_size
         y_coords = galaxy_catalog.get_col_data("y") - half_stamp_size 

         angle_pos_dico = dict(zip(zip(x_coords, y_coords), zip(x_field_true_deg, y_field_true_deg)))

      galaxy_catalog.close()

      return angle_pos_dico

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_image_filepath(self, job, worker):
      """! 
         Get the full path of the image containing the galaxy postrage stamps
         @param job GfitJob object
         @param worker Worker process object
         @return full path of the image containing the galaxy postrage stamps
      """
      return job.get_file_path(self.image_file_type)

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_catalog_filepath(self, job, worker):
      """! 
         Get the full path of the galaxy catalog containing at least the coordinates of the galaxies
         @param job GfitJob object
         @param worker Worker process object
         @return full path of the catalog containing the galaxy coordinates
      """

      catalog_filepath = None
      path_dico = job.img_path_dico.get(job.branch_key, None)
      if not path_dico is None:
         file_types = [p for p in path_dico.keys() if job.dataset.is_galaxy_catalog(worker, p)] 
         if len(file_types) > 0:
            # Note: we assume there is only one galaxy catalog          
            catalog_filepath = path_dico.get(file_types[0], None)  

      return catalog_filepath

   # -----------------------------------------------------------------------------------------------
   def _get_galaxy_se_catalog_filepath(self, file_type, job, worker):
      """! 
         Get the full path of the galaxy catalog containing SExtractor information about the 
         galaxy objects
         @param job GfitJob object
         @param worker Worker process object
         @return full path of the SExtractor catalog
      """
      #print "*** galaxy_se_filepath_dico:", job.galaxy_se_filepath_dico

      if len(job.galaxy_se_filepath_dico) > 0 and \
         (job.img_no, job.epoch) in job.galaxy_se_filepath_dico:
         return job.galaxy_se_filepath_dico.get((job.img_no, job.epoch), None)[0] 
      else:
         return None  

#      [filetype_main, filetype_ext] = os.path.splitext(file_type)
#      if len(job.galaxy_se_filepath_dico) > 0 and \
#         (filetype_main, job.img_no, job.epoch) in job.galaxy_se_filepath_dico:
#         return job.galaxy_se_filepath_dico.get((filetype_main, job.img_no, job.epoch), None) 
#      else:
#         return None   

   # -----------------------------------------------------------------------------------------------
   def _get_psf_se_catalog_filepath(self, file_type, job, worker):
      """! 
         Get the full path of the PSF catalog containing SExtractor information about the 
         PSF objects
         @param job GfitJob object
         @param worker Worker process object
         @return full path of the SExtractor catalog
      """

      #print "*** star_se_filepath_dico:", job.star_se_filepath_dico

      if len(job.star_se_filepath_dico) > 0 and \
         (job.img_no, job.epoch) in job.star_se_filepath_dico:
         return job.star_se_filepath_dico.get((job.img_no, job.epoch), None)[0] 
      else:
         return None  

#      [filetype_main, filetype_ext] = os.path.splitext(file_type)
#      if len(job.star_se_filepath_dico) > 0 and \
#         (filetype_main, job.img_no, job.epoch) in job.star_se_filepath_dico:
#         return job.star_se_filepath_dico.get((filetype_main, job.img_no, job.epoch), None) 
#      else:
#         return None  

   # -----------------------------------------------------------------------------------------------
   def _get_psf_image_filepath(self, file_type, job, worker):
      """! 
         Get the full path of the PSF field image with postage stamps matching those of the 
         galaxy images
         @param job GfitJob object
         @param worker Worker process object
         @return full path of the PSF image
      """

#      print "*** psf_filepath_dico:", job.psf_filepath_dico

      if len(job.psf_filepath_dico) > 0 and (job.img_no, job.epoch) in job.psf_filepath_dico:
         return job.psf_filepath_dico[(job.img_no, job.epoch)][0]    # PSF image path
      else:
#         if logging_enabled:
#            worker.logger.log_warning_p(
#                   "{0} - Could not find PSF image for: /{1}/image-{2:03d}-{3:1d} - {4}"\
#                   " => Job will not be processed".format(
#                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, image_file_type))
         return None

#      [filetype_main, filetype_ext] = os.path.splitext(file_type)
#      print "====== PSF key:", (filetype_main, job.img_no, job.epoch)      
#      if len(job.psf_filepath_dico) > 0 and \
#         (filetype_main, job.img_no, job.epoch) in job.psf_filepath_dico:

#         return job.psf_filepath_dico.get((filetype_main, job.img_no, job.epoch))  # PSF image path
#      else:
##         if logging_enabled:
##            worker.logger.log_warning_p(
##                   "{0} - Could not find PSF image for: /{1}/image-{2:03d}-{3:1d} - {4}"\
##                   " => Job will not be processed".format(
##                      worker.name, job.get_branch_tree(), job.img_no, job.epoch, image_file_type))
#         return None

#   # -----------------------------------------------------------------------------------------------
#   def _get_galaxy_effective_stamp_size(self, job, worker):
#      """! 
#         Get the the size of the galaxy stamp to use for shape measurement
#         @param job GfitJob object
#         @param worker Worker process object
#         @return effective stamp size (necessarilly lower than full galaxy stamp size)
#      """         
#      stamp_size_dico = worker.config.get_as_dict("GALAXY_EFFECTIVE_STAMP_SIZE", "PRIMARY_DATASET")
#      return int(stamp_size_dico[job.branch][job.obs_type])

#   # -----------------------------------------------------------------------------------------------
#   def _get_psf_effective_stamp_size(self, job, worker):
#      """! 
#         Get the the size of the PSF stamp to use for shape measurement
#         @param job GfitJob object
#         @param worker Worker process object
#         @return effective stamp size (necessarilly lower than full PSF stamp size)
#      """         
#      stamp_size_dico = worker.config.get_as_dict("PSF_EFFECTIVE_STAMP_SIZE", "PRIMARY_DATASET")
#      return int(stamp_size_dico[job.branch][job.obs_type])


   # -----------------------------------------------------------------------------------------------
   def _open_catalog(self, filepath, is_sextractor, hdu_no, job, worker):

      # --- Read ID data from galaxy catalog
      if job.dataset.is_fits(worker, filepath):
         catalog = FITSCatalog(filepath, hdu_no=hdu_no, mem_map=True)
      else:
         if is_sextractor:
            catalog = SExCatalog(filepath)
         else:
            catalog = TextCatalog(filepath)

      catalog.open()

      return catalog



# -- EOF gfit_shape.py
