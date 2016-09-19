"""! 
   MultiFit - Bulge-Disk Sersic model using GalSim
"""

# -- Python imports
import math
import numpy

# --- External imports
import galsim
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models', config_filename="gbdsersic.cfg"):
      """! 
         Construct a @c gbdsersic FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="gbdsersic", 
                                      config_dir=config_dir, config_filename=config_filename)
      
      # --- Set parameter definitions from configuration
      self.set_params(self.get_params_from_config() )

      # --- Set module parametric function
      if self.config.has_section("MODEL.FUNC"):

         # --- Use the profile function specified in the configuration 
         func_name = self.name
         if self.config.has_key("func_name", "MODEL.FUNC"):
            func_name = self.config.get_as_string("func_name", "MODEL.FUNC")
         self.func = eval("self."+func_name)

      else:
         # --- Look for a profile function with same name as the model 

         self.func = self.gbdsersic

      # --- Read Galsim-specific configuration if available
      if self.config.has_section("MODEL.GALSIM"):

         # If True, abort the fit when an error occurs; otherwise tries to find another best fit
         self._abort_fit_uppon_error = self.config.get_as_boolean("abort_fit_uppon_error", 
                                                                  "MODEL.GALSIM")

         # Smallest allowed bulge effective radius; beyond, assumes single disk component galaxy
         self._min_bulge_radius = self.config.get_as_float("min_bulge_sersic_radius", 
                                                           "MODEL.GALSIM")
      else:
         self._abort_fit_uppon_error = True
   

#      # --- Read catalog column mapping iformation if available
#      self._col_mapping_dico = {}
#      if self.config.has_section("MODEL.GUESS"):

#         # Read catalog column mapping information
#         if self.config.has_key("guess_col_mapping", "MODEL.GUESS"):
#            self._col_mapping_dico = self.config.get_as_dict("guess_col_mapping", "MODEL.GUESS") 
            

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters     

   @property
   def out_of_bound_handler(self):
      """! 
         Return the callback handler function to notify a parameter value came outside the 
                 allowed range of varation specified by Param.guess_value
         @return the callback handler function
      """
      #print "***MFIT: gbdsersic, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: gbdsersic, getting  numerical_error_handler..."

      return self._numerical_error_handler


   @property
   def abort_fit_uppon_error(self):
      """! 
         Tells whether a fit should be aborted when an error occurs
         @return whethed a fit should be aborted when an error occurs
      """
      return self._abort_fit_uppon_error

   @property
   def min_bulge_radius(self):
      """! 
         Return the smallest allowed galaxy bulge effective radius
         @return the smallest allowed galaxy bulge effective radius
      """
      return self._min_bulge_radius


   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_model_func_from_config(self):
      """!
         Return a pointer to the model function for evaluating the models given parameters
         @return model function
         @note this method is called by multifit.BaseFittingModel and must be implemented
      """
      if self.config.has_section("MODEL.FUNC"):

         # --- Use the profile function specified in the configuration 
         model_func_name = self.name
         if self.config.has_key("func_name", "MODEL.FUNC"):
            model_func_name = self.config.get_as_string("func_name", "MODEL.FUNC")
         model_func = eval("self."+model_func_name)

      else:
         # --- Look for a profile function with same name as the model 

         model_func = self.gbdsersic

      return model_func 


   # -----------------------------------------------------------------------------------------------
   def get_params_from_config(self):
      """! 
         Obtain the list of parameter definitions from the configuration file. 
         In this particular case, the model is generic and parameters must be set programmatically 
      """

      return BaseFittingModel.get_params_from_config(self)


   # -----------------------------------------------------------------------------------------------
   def get_best_guess_params(self, coords, actual_galaxy_image, 
                                   galaxy_info_dico, se_galaxy_info_dico):
      """
         Attempt to provide estimates for parameter guess values based on input galaxy data
         @param coords coordinates of left-bottom corner of postage stamp (not in numpy convention)
         @param actual_galaxy_image galaxy postage stamp to fit
         @param se_galaxy_info_dico galaxy data dictionary (depends on dataset)
         @param se_galaxy_info_dico SExtractor data dictionary
         @return list of ModelParam objects with estimated guess values
         @todo possibly estimate (e1,e2) using quadrupole moments
      """   

      for param in self.params:

         if param.guess_value is None:

            if param.name == "flux":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_flux", 
                                                         se_galaxy_info_dico["SE_GAL_FLUX"])
            elif param.name == "dre":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_dre", 
                                                         se_galaxy_info_dico["SE_GAL_FLUX_RADIUS"])
            elif param.name == "bre":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_bre", 
                                                         se_galaxy_info_dico["SE_GAL_FLUX_RADIUS"])
            elif param.name == "xc": 
               if "guess_xc" in self.col_mapping_dico:
                  param.guess_value = self.get_guess_value(coords, 
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            param.name, "guess_xc",
                                                            (actual_galaxy_image.shape[1]-1.0)/2.0)
               else:
                  param.guess_value = (actual_galaxy_image.shape[1]-1.0)/2.0
            elif param.name == "yc": 
               if "guess_yc" in self.col_mapping_dico:
                  param.guess_value = self.get_guess_value(coords, 
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            param.name, "guess_yc", 
                                                            (actual_galaxy_image.shape[0]-1.0)/2.0)
               else:
                  param.guess_value = (actual_galaxy_image.shape[0]-1.0)/2.0

            elif param.name == "e1" or param.name == "e2":
               param.guess_value = 0.0

      return self.params

      
   # -----------------------------------------------------------------------------------------------
   def check_fitted_params(self, (x, y), actual_galaxy_image, 
                                         galaxy_info_dico, se_galaxy_info_dico, 
                                         guess_values, param_values):
      """!
         Check that the parameter values  @c param_values are within the allowed bounds
         @param coords coordinates of left-bottom corner of postage stamp (not in numpy convention)
         @param actual_galaxy_image galaxy postage stamp to fit
         @param se_galaxy_info_dico galaxy data dictionary (depends on dataset)
         @param se_galaxy_info_dico SExtractor data dictionaryds
         @return True if the valiity check is successful 
      """

      is_valid = BaseFittingModel.check_fitted_params(self, (x, y), actual_galaxy_image,
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            guess_values, param_values)

      if is_valid:
         # Check global distortion < 1
         e1 = param_values[self.get_param_names().index("e1")]
         e2 = param_values[self.get_param_names().index("e2")]
         is_valid = math.hypot(e1, e2) < 1.0
   
         #print "***MFIT: is_valid", is_valid, e1, e2, math.hypot(e1, e2)
      
      return is_valid

   # -----------------------------------------------------------------------------------------------
   def out_of_bound_handler(self, iiter, param, param_value, bound_index):
         
   #   if param.name == 'e1':
   #      status = 0  # continue fitting => may reach < +-1 later on
   #
   #   elif param.name == 'e2':
   #      status = 0  # continue fitting => may reach < +-1 later on

   #   elif param.name == 'dre' and bound_index == 0:
   #      status = -1
   #   elif param.name == 'dre' and bound_index == 1:
   #      status = -1
   #   elif param.name == 'bre' and bound_index == 0:
   #      status = -1
   #   elif param.name == 'bre' and bound_index == 1:
   #      status = -1
   #   else:
   #      status = 0  # continue fitting

      #print "***MFIT: in gbdsersic out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #print "***MFIT: in gbdsersic numerical_error_handler", iiter, param, value, value_name

      if self.abort_fit_uppon_error:
         return -1   # always abort fitting
      else:
         return 0    # attempt to fix 


  # -----------------------------------------------------------------------------------------------
   def gbdsersic(self, model, stamp_size, subscale, flux, dre, bre, df, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print flux, dn, dre, bre, df, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # --- Galaxy is modelled as a bulge + disk profile
         disk = galsim.Exponential(flux=df, half_light_radius=dre)
         disk.applyShear(e1=e1, e2=e2)

         # --- The flux of an object is the sum of the disk & bulge component fluxes
         bulge = galsim.DeVaucouleurs(flux=1.0-df, half_light_radius=bre)
         bulge.applyShear(e1=e1, e2=e2)
         gal = galsim.Add([disk, bulge])

         # --- Apply shear (distortion in this case)
         gal.applyShear(e1=e1, e2=e2)

         # --- The total flux of bulge+disk is set to flux
         gal.setFlux(flux)

         # --- Shift the object based on the estimated centroid
         half = stamp_size / 2.0
         gal.applyShift(dx=xc-half, dy=yc-half)

         # --- Render the image
         gal_image = galsim.ImageF(stamp_size, stamp_size, scale=1.0)
         gal.draw(gal_image, dx=1.0)

#         # Render the image
#         pixel_scale = 1.0
#         pix = galsim.Pixel(pixel_scale)
#         gal = galsim.Convolve([pix, gal])

#         gal_image = galsim.ImageF(stamp_size, stamp_size, pixel_scale)
#         gal.draw(gal_image, dx=pixel_scale)

         # --- Image data as a numpy array
         gal_array = gal_image.array
   
         #gal_image.write("bd.fits")


      except Exception as detail:

         gal_array = None
         #gal_array = numpy.nan  

         # --- Abort the fit
         print("***  Multifit: {0} *** An error occurred: {1}".format(model.name, detail))
         #self.helper.print_error("gbdsersic model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  gal_array



