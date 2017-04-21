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

   def __init__(self, config_dir='./config/models', config_filename="gbdshsersic.cfg"):
      """! 
         Construct a @c gbdshsersic FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="gbdshsersic", 
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

         self.func = self.gbdshsersic

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
      #print "***MFIT: gbdshsersic, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: gbdshsersic, getting  numerical_error_handler..."

      return self._numerical_error_handler

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

         model_func = self.gbdshsersic

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
   def check_fitted_params(self, (x, y), full_galaxy_image, actual_galaxy_image, 
                                         galaxy_info_dico, se_galaxy_info_dico, 
                                         guess_values, param_values):
      """!
         Check that the parameter values  @c param_values are within the allowed bounds
         @param coords coordinates of left-bottom corner of postage stamp (not in numpy convention)
         @param full_galaxy_image full-size galaxy postage stamp to fit
         @param actual_galaxy_image galaxy postage stamp to fit
         @param se_galaxy_info_dico galaxy data dictionary (depends on dataset)
         @param se_galaxy_info_dico SExtractor data dictionaryds
         @return True if the valiity check is successful 
      """

      is_valid = BaseFittingModel.check_fitted_params(self, (x, y), 
                                                            full_galaxy_image, actual_galaxy_image,
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

      #print "***MFIT: in gbdshsersic out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #print "***MFIT: in gbdshsersic numerical_error_handler", iiter, param, value, value_name

      if self.config.has_section("MODEL.GALSIM") and \
         self.config.get_as_boolean("abort_fit_uppon_error", "MODEL.GALSIM"):
         return -1   # always abort fitting
      else:
         return 0    # attempt to fix 


   # -----------------------------------------------------------------------------------------------
   def gbdshsersic(self, model, stamp_size, subscale, flux, dre, bre, df, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print flux, dn, dre, bre, df, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # --- Adjust parameters based on resolution and setup model pixel grid
         if subscale > 1:
            stamp_size *= subscale
            dre *= subscale
            bre *= subscale
            corr = (subscale - 1)/ (2.0 * subscale)    # account for centroid shift
            xc = (xc + corr) * subscale
            yc = (yc + corr) * subscale

         # Galaxy is modelled as a bulge + disk profile
         disk  = galsim.Exponential(flux=df, half_light_radius=dre)
         bulge = galsim.DeVaucouleurs(flux=1.0-df, half_light_radius=bre)

         # The flux of an Add object is the sum of the component fluxes
         gal = galsim.Add([disk, bulge])

         # Apply shear to both components
         gal.applyShear(g1=e1, g2=e2)

         # The total flux of bulge+disk is set to flux
         gal.setFlux(flux)

         # Shift the object based on the estimated centroid
         half = (stamp_size-1) / 2.0
         gal.applyShift(dx=xc-half, dy=yc-half)

         # Render the image
         gal_image = galsim.ImageF(stamp_size, stamp_size, scale=1.0)
         gal.draw(gal_image, dx=1.0)

         # Image data as a numpy array
         gal_array = gal_image.array
   
         #gal_image.write("bd.fits")

         # Subscaling if requested
         if subscale > 1:
            gal_array = self._create_low_resolution_stamp(gal_array, stamp_size)

      except:

         gal_array = None

         # --- Abort the fit
         print("*** gbdshsersic model *** An error occurred: {0}".format(sys.exc_info()))
         #self.helper.print_error("gbdshsersic model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  gal_array


