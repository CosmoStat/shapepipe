"""! 
   Various helper functions
""" 

# -- Python imports
import os, sys
import math
import numpy
import scipy
import pyfits
import shutil

# --- Python imports
try:
    from thread import get_ident as _get_ident
except ImportError:
    from dummy_thread import get_ident as _get_ident
from _abcoll import KeysView, ValuesView, ItemsView

# -- External import
from scatalog import *  # file catalog management             
from mpfg.mp_helper import *  # base Helper

# -- Module-specific imports  
from gfit_helper import *  # this module helper utility functions
from gfit_version import __version__  # version information

# -------------------------------------------------------------------------------------------------
class GfitHelper(Helper):

   """! Convenient utility functions """

   # -----------------------------------------------------------------------------------------------
   def get_version(self):
      """! 
          Get the version number of this pse code as text 
          @return version number of this pse code as text
      """
      return __version__ 

   # -----------------------------------------------------------------------------------------------
   def show_config_summary(self, master):
      
      if master.logging_enabled():

         # --- gfit version
         master.logger.log_info_p("\n*** gfit {0} ***\n".format(self.get_version()))

         # --- Python modules
         master.logger.log_info_p("Standard Python modules:")
         try:
            np = __import__('numpy', globals(), locals(), [], -1)
            master.logger.log_info_p("- numpy {0}\t\t{1}".format(np.__version__, np.__file__))        
            sp = __import__('scipy', globals(), locals(), [], -1)
            master.logger.log_info_p("- scipy {0}\t\t{1}".format(sp.__version__, sp.__file__))           
            pyfits = __import__('pyfits', globals(), locals(), [], -1)
            master.logger.log_info_p("- pyfits {0}\t\t{1}".format(pyfits.__version__, pyfits.__file__))     
            mpl = __import__('matplotlib', globals(), locals(), [], -1)
            master.logger.log_info_p("- matplotlib {0}\t{1}".format(mpl.__version__, mpl.__file__))  
            galsim = __import__('galsim', globals(), locals(), [], -1)
            master.logger.log_info_p("- galsim {0}\t\t{1}".format(galsim.version, galsim.__file__))  
            # asc = __import__('asciidata', globals(), locals(), [], -1)
            # master.logger.log_info_p("- asciidata {0}\t{1}".format(asc.__version__, asc.__file__))  
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be found: {0}\n".format(detail)) 

         master.logger.log_info_p("\nMPF Python modules:")
         try:
            mpfg = __import__('mpfg', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfg {0}\t\t{1}".format(mpfg.__version__, mpfg.__file__))  
            mpfx = __import__('mpfx', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfx {0}\t\t{1}".format(mpfx.__version__, mpfx.__file__))  
            slog = __import__('slogger', globals(), locals(), [], -1)
            master.logger.log_info_p("- slogger {0}\t\t{1}".format(slog.__version__, slog.__file__))  
            sconf = __import__('sconfig', globals(), locals(), [], -1)
            master.logger.log_info_p("- sconfig {0}\t\t{1}".format(sconf.__version__, sconf.__file__))  
            scat = __import__('scatalog', globals(), locals(), [], -1)
            master.logger.log_info_p("- scatalog {0}\t{1}".format(scat.__version__, scat.__file__)) 
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be found: {0}\n".format(detail)) 

         master.logger.log_info_p("\nFitting modules:")
         try:
            multi = __import__('multifit', globals(), locals(), [], -1)
            master.logger.log_info_p("- multifit {0}\t{1}".format(multi.__version__, multi.__file__))    
            scdm = __import__('scdm', globals(), locals(), [], -1)
            master.logger.log_info_p("- scdm {0}\t\t{1}".format(scdm.__version__, scdm.__file__))
            kapteyn = __import__('kapteyn.kmpfit', globals(), locals(), [], -1)
            master.logger.log_info_p("- kmpfit {0}\t\t{1}".format(kapteyn.__version__, kapteyn.__file__))
            
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be imported: {0}\n".format(detail)) 

#         try:
#            sklearn = __import__('sklearn', globals(), locals(), [], -1)
#            master.logger.log_info_p("- scikit learn {0}\n".format(sklearn.__version__))  
#         except:
#            pass

         # --- configuration
         gfit_config_filename = master.config_filename
         gfit_config_dir = master.config_dir
         gfit_config_path = os.path.join(gfit_config_dir, gfit_config_filename)

         model_config_path = master.config.get_as_string("MODEL_CONFIG_PATH", "MULTIFIT")
         method_config_path = master.config.get_as_string("METHOD_CONFIG_PATH", "MULTIFIT")
         fitting_config_path = master.config.get_as_string("FITTING_CONFIG_PATH", "MULTIFIT")

         master.logger.log_info_p("\nBase gfit configuration:")
         master.logger.log_info_p("- Gfit main config : {0}".format(gfit_config_path))
         master.logger.log_info_p("- Fitting methods  : {0}".format(method_config_path))
         master.logger.log_info_p("- Galaxy models    : {0}".format(model_config_path))
         master.logger.log_info_p("- Custom fitting   : {0}".format(fitting_config_path))

         # --- Multifit
         mfit_model_module_path = master.config.get_as_string("MODEL_MODULE_PATH", "MULTIFIT")   
         mfit_method_module_path = master.config.get_as_string("METHOD_MODULE_PATH", "MULTIFIT")   

         master.logger.log_info_p("\nMultifit configuration directories:")
         master.logger.log_info_p("- Models : {0}".format(mfit_model_module_path))
         master.logger.log_info_p("- Methods: {0}\n".format(mfit_method_module_path))

         # --- Selected minimizer & galaxy model
         galaxy_model = master.config.get_as_string("GALAXY_MODEL_NAME", "GALAXY_FITTING")
         fitting_method = master.config.get_as_string("GALAXY_FITTING_METHOD", "GALAXY_FITTING")
         master.logger.log_info_p("Selected galaxy model : {0}".format(galaxy_model))
         master.logger.log_info_p("Selected minimizer    : {0}".format(fitting_method))

         # --- Directories
         # master.logger.log_info_p("Base input directory   : {0}".format(master.base_input_dir))
         master.logger.log_info_p("\nBase output directory: {0}".format(master.base_output_dir))
         master.logger.log_info_p("- Run output directory   : {0}".format(master.run_output_dir))
         master.logger.log_info_p("- log output directory   : {0}".format(master.log_output_dir))
         master.logger.log_info_p("- plot output directory  : {0}".format(master.plot_output_dir))
         master.logger.log_info_p("- stats output directory : {0}".format(master.stat_output_dir))
         master.logger.log_info_p("- error output directory : {0}\n".format(master.error_output_dir))

         # --- Datasets   
         try:
            src_dataset_dir = master.config.get_as_string("BASE_DIR", "PRIMARY_DATASET")
            gal_se_dataset_dir = master.config.get_as_string("BASE_DIR",
                                                              "SEXTRACTOR_DATASET_GALAXY")
            star_se_dataset_dir = master.config.get_as_string("BASE_DIR",
                                                              "SEXTRACTOR_DATASET_PSF")
            psf_dataset_dir = master.config.get_as_string("BASE_DIR", "PSF_DATASET")

            master.logger.log_info_p("Primary dataset           : {0}".format(src_dataset_dir))
            master.logger.log_info_p("PSF dataset               : {0}".format(psf_dataset_dir))
            master.logger.log_info_p("SExtractor galaxy dataset : {0}".format(gal_se_dataset_dir))
            master.logger.log_info_p("SExtractor psf dataset    : {0}\n".format(star_se_dataset_dir))
         except Exception as detail:
            master.logger.log_error_p("Error while getting dataset information: {0}\n".format(detail)) 
   
   # -----------------------------------------------------------------------------------------------
   def read_fits_image(self, image_filepath, hdu=0, float_size="float32"):

      # --- Open catalog 
      if not os.path.exists(image_filepath):
         print 'Error: catalog %s not found' % (catalog_filepath)
         return False, None, None, {}

      fits_image = pyfits.open(image_filepath, mode='readonly', memmap=True)
      image_data = fits_image[hdu].data.astype(float_size)
      image_header = fits_image[hdu].header
      fits_image.close()

      return image_data, image_header

   # -----------------------------------------------------------------------------------------------
   def get_stamp_center_as_int(self, stamp_size):
      """ ! @return the postage stamp geometrical center pixel no, indexed from zero """      
      
      return int((stamp_size - 1.0) / 2.0)

   # -----------------------------------------------------------------------------------------------
   def cut_stamp_around_centroid(self, stamp, stamp_size, (row, col)):
      """!  Cut a postage stamp of a given size around a given centroid relative to the image """
      
      half_stamp_size = stamp_size / 2.0
      if stamp_size % 2 == 0:     
         # --- Even stamp size
         # print "EVEN stamp size:", stamp_size, "centroid:", (row, col), half_stamp_size, row-half_stamp_size
         return stamp[int(row - half_stamp_size):int(row + half_stamp_size),
                      int(col - half_stamp_size):int(col + half_stamp_size)]
      else:
         # --- Odd stamp size
         # print "ODD stamp size:", stamp_size, "centroid:", (row, col), half_stamp_size, row-half_stamp_size+0.5
         return stamp[int(row - half_stamp_size + 0.5):int(row + half_stamp_size + 0.5),
                      int(col - half_stamp_size + 0.5):int(col + half_stamp_size + 0.5)]


   # -----------------------------------------------------------------------------------------------
   def empty_stamp(self, stamp, border_size=1):
      """! 
         Empty a square stamp, leaving only the border. The original stamp is untouched
      """

      stamp = stamp.copy()
      stamp[border_size:stamp.shape[0] - border_size, border_size:stamp.shape[0] - border_size] = 0.0
      return stamp

   # -----------------------------------------------------------------------------------------------
   def get_stamp_border(self, stamp, border_size=1):
 
      """!
         get border from stamp
      """
      return stamp[border_size:stamp.shape[0] - border_size, border_size:stamp.shape[0] - border_size]   

#    # -----------------------------------------------------------------------------------------------
#    def get_stamp_border(self, stamp):
# 
#       """!
#          Extract border from stamp
#       """
# 
#       nn, mm = stamp.shape
#       a1 = stamp[0,:]
#       a2 = stamp[nn-1,:]
#       a3 = stamp[1:nn-1,0]
#       a4 = stamp[1:nn-1,mm-1]
#       border = numpy.concatenate((a1,a2))
#       border = numpy.concatenate((border,a3))
#       border = numpy.concatenate((border,a4))
# 
#       return border

#   # -----------------------------------------------------------------------------------------------
#   def get_stamp_border(self, stamp, border_size=1):
#      """!
#         Empy a square stamp, leaving only the border
#      """
#      e_stamp = self.empty_stamp(stamp, border_size)
#      return numpy.extract(numpy.not_equal(e_stamp, 0.0), e_stamp)

   # -----------------------------------------------------------------------------------------------
   def compute_stamp_sigma_back_noise(self, stamp):
      """!
         Get stamp background median noise evaluated from its border
      """
      border = self.get_stamp_border(stamp)
      return math.sqrt(numpy.mean(border ** 2))
#      return (numpy.sum((stamp_border - numpy.median(stamp_border.flat))**2) / len(stamp_border.flat))**0.5

   # -----------------------------------------------------------------------------------------------
   def compute_stamp_signal(self, stamp):
      """! Estimate the signal in a postage stamp image """
      return numpy.sum(stamp ** 2) / numpy.sum(stamp)

   # -----------------------------------------------------------------------------------------------
   def estimate_stamp_noise(self, request, stamp, section_name, job, worker):
      """
         ! Estimate the noise in a postage stamp image
         @param request ShapeMeasurementRequest object
         @param stamp object postage stamp (2D)
         @param section_name either 'GALAXY_NOISE' or 'PSF_NOISE' (see configuration file)
         @param job this job
         @param worker worker process   
      """

      noise_map_model_name = worker.config.get_as_string("NOISE_MODEL_NAME", section_name)
      if noise_map_model_name == "CONST_SIGMA":
         # --- Fixed sigma value   
         sigma_value = worker.config.get_as_float("SIGMA_VALUE",
                                          "{0}.{1}".format(section_name, noise_map_model_name))

      elif noise_map_model_name == "MATCHED_FILTER":  # kept / G3 compatibility: does not work well
         # ------------------------------------------------------------
         # Sigma value calculated as follows (see GREAT08, GREAT3):
         # S = sum W(x,y) I(x,y) / sum W(x,y)
         # N^2 = Var(S) = sum W(x,y)^2 Var(I(x,y)) / (sum W(x,y))^2
         # sigma = N = Sqrt(Var(S)) 
         # ------------------------------------------------------------
         sigma_value = scipy.std(stamp) * math.sqrt(numpy.sum(stamp ** 2)) / numpy.sum(stamp)         
      
      else:   
         # --- By default, noise estimated as the median value of the postage stamp border   
         border_size = worker.config.get_as_int("BORDER_SIZE",
                                          "{0}.{1}".format(section_name, noise_map_model_name)) 
         border = numpy.ravel(stamp[border_size:stamp.shape[0] - border_size, \
                                    border_size:stamp.shape[0] - border_size])
         sigma_value = math.sqrt(numpy.median(border ** 2))          

      # TODO: possibly add more sophisticated noise estimation models here...    
      
      return sigma_value
         

#    # -----------------------------------------------------------------------------------------------
#    def compute_stamp_noise(self, stamp):
#       """! Estimate the noise in a postage stamp image """
#       # ------------------------------------------------------------
#       # Sigma value calculated as follows (see GREAT08, GREAT3):
#       # S = sum W(x,y) I(x,y) / sum W(x,y)
#       # N^2 = Var(S) = sum W(x,y)^2 Var(I(x,y)) / (sum W(x,y))^2
#       # sigma = N = Sqrt(Vat(S)) 
#       # ------------------------------------------------------------
#       return scipy.std(stamp) * math.sqrt(numpy.sum(stamp**2)) / numpy.sum(stamp)

   # -----------------------------------------------------------------------------------------------
   def compute_stamp_SNR(self, stamp):
      """! Compute the SNR of a postage stamp image """
      return math.sqrt(numpy.sum(stamp ** 2)) / scipy.std(stamp)

   # -----------------------------------------------------------------------------------------------
   def compute_stamp_sky_background(self, stamp):
      """!
         Get stamp median sky background evaluated from the stamp border
      """   
      # return numpy.median(get_stamp_border(stamp).flat)
      return numpy.mean(numpy.ravel(self.get_stamp_border(stamp)))

#    # -----------------------------------------------------------------------------------------------
#    def estimate_sigma_noise(self, image):
#       """
#       Estimate noise sigma from image border 
#       """
#       nn, mm = image.shape
#       a1 = image[0,:]
#       a2 = image[nn-1,:]
#       a3 = image[1:nn-1,0]
#       a4 = image[1:nn-1,mm-1]
#       a1 = np.concatenate((a1,a2))
#       a1 = np.concatenate((a1,a3))
#       a1 = np.concatenate((a1,a4))
#       a4 = (a1**2).mean()
#       return math.sqrt(a4)

   # -----------------------------------------------------------------------------------------------
   def make_galaxy_sky_image(self, request, image, job, worker):

      sky_image = numpy.zeros(image.shape)
      sky_model_name = worker.config.get_as_string("SKY_MODEL_NAME", "SKY_MODEL_GALAXY")
      section_name = "SKY_MODEL_GALAXY.{0}".format(sky_model_name)
      
      if sky_model_name == "CONST_SKY":
         # --- Sky has the same value everywhere
         sky_value = worker.config.get_as_float("SKY_VALUE", section_name)
         sky_image = sky_value * numpy.ones(image.shape)
           
      elif sky_model_name == "GAUSSIAN_NOISE":     
         # --- Sky normal distribution with given mean and variance
         sky_mean = worker.config.get_as_float("SKY_MEAN", section_name) 
         sky_sigma = worker.config.get_as_float("SKY_SIGMA", section_name) 
         sky_image = numpy.random.randn(image.shape[0], image.shape[1]) * sky_sigma + sky_mean  
         
#       elif sky_model_name == "BORDER_MEDIAN":      
#          # --- Sky estimated as the median value of the image border   
#          border_size = worker.config.get_as_int("BORDER_SIZE", section_name)   
#          border = numpy.ravel(image[border_size:image.shape[0]-border_size,\
#                                     border_size:image.shape[0]-border_size])
#          sky_image = math.sqrt(numpy.median(border**2)) * numpy.ones(image.shape) 
            
      elif sky_model_name == "EXT_IMAGE":      
         # --- External image used as sky background
         sky_image, _ = self.read_fits_image(
                           request.galaxy_sky_image_filepath,
                           hdu=request.galaxy_sky_image_property_dico["HDU_NO"],
                           float_size=request.galaxy_sky_image_property_dico["PIXEL_FLOAT_SIZE"])
      else:                 
         if worker.logging_enabled():
            worker.logger.log_warning_p(
               "{0} - {1}/image-{2:03d}-{3:1d} - "\
               "unknown galaxy Sky model {4}, ignoring...".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch,
                   sky_model_name))

      return sky_image           

   # -----------------------------------------------------------------------------------------------
   def make_psf_sky_image(self, request, image, job, worker):

      sky_image = numpy.zeros(image.shape)
      sky_model_name = worker.config.get_as_string("SKY_MODEL_NAME", "SKY_MODEL_PSF")
      section_name = "SKY_MODEL_PSF.{0}".format(sky_model_name)
      
      if sky_model_name == "CONST_SKY":
         # --- Sky has the same value everywhere
         sky_value = worker.config.get_as_float("SKY_VALUE", section_name)
         sky_image = sky_value * numpy.ones(image.shape)
           
      elif sky_model_name == "GAUSSIAN_NOISE":     
         # --- Sky normal distribution with given mean and variance
         sky_mean = worker.config.get_as_float("SKY_MEAN", section_name) 
         sky_sigma = worker.config.get_as_float("SKY_SIGMA", section_name) 
         sky_image = numpy.random.randn(image_shape[0], image_shape[1]) * sky_sigma + sky_mean  
         
      elif sky_model_name == "EXT_IMAGE":      
         # --- External image used as sky background     
         sky_image, _ = self.read_fits_image(
                           request.psf_sky_image_filepath,
                           hdu=request.psf_sky_image_property_dico["HDU_NO"],
                           float_size=request.psf_sky_image_property_dico["PIXEL_FLOAT_SIZE"])
                       
      else:
         if worker.logging_enabled():
            worker.logger.log_warning_p(
               "{0} - {1}/image-{2:03d}-{3:1d} - "\
               "unknown PSF Sky model {4}, ignoring...".format(
                   worker.name, job.get_branch_tree(), job.img_no, job.epoch,
                   sky_model_name))
                        
      return sky_image      

   # -----------------------------------------------------------------------------------------------
   def get_galaxy_noise_map(self, request, image, job, worker):
      """! Create or Read a RMS noise map, whose values are the estimated standard deviation 
           (pixels represent sigma values) of the galaxy noise
      """     

      noise_map_image = None 
      
      noise_map_model_name = worker.config.get_as_string("NOISE_MODEL_NAME", "GALAXY_NOISE")
      if noise_map_model_name == "EXT_NOISE_MAP":
         # --- External image used as rms noise map
         noise_map_image, _ = self.read_fits_image(
                     request.galaxy_noise_map_filepath,
                     hdu=request.galaxy_noise_map_property_dico["HDU_NO"],
                     float_size=request.galaxy_noise_map_property_dico["PIXEL_FLOAT_SIZE"])
      
      return noise_map_image           

   # -----------------------------------------------------------------------------------------------
   def get_psf_noise_map(self, request, image, job, worker):
      """! Create or Read a RMS noise map, whose values are the estimated standard deviation 
           (pixels represent sigma values) of the PSF noise
      """     

      noise_map_image = None 
      
      noise_map_model_name = worker.config.get_as_string("NOISE_MODEL_NAME", "PSF_NOISE")
      if noise_map_model_name == "EXT_NOISE_MAP":
         # --- External image used as rms noise map
         noise_map_image, _ = self.read_fits_image(
                     request.psf_noise_map_filepath,
                     hdu=request.psf_noise_map_property_dico["HDU_NO"],
                     float_size=request.psf_noise_map_property_dico["PIXEL_FLOAT_SIZE"])
      
      return noise_map_image           

   # -----------------------------------------------------------------------------------------------
   def get_image_property_dico(self, section_name, job, worker):   
      """! Create or Read a RMS noise map, whose values are the estimated standard deviation 
           (pixels represent sigma values) of the PSF noise
      """     

      section_dico = {"DIR_INPUT_LIST":[], "DIR_RECURSE":True, \
                      "BASE_DIR":".", "FILE_PATTERNS":[], \
                      "HDU_NO":0, "PIXEL_FLOAT_SIZE":"float32"\
                     }
      if worker.config.has_section(section_name):
         if worker.config.has_key("FILE_PATTERNS", section_name):
            section_dico["FILE_PATTERNS"] = worker.config.get_as_list("FILE_PATTERNS",
                                                                      section_name)         
         if worker.config.has_key("DIR_INPUT_LIST", section_name):
            section_dico["DIR_INPUT_LIST"] = worker.config.get_as_list("DIR_INPUT_LIST",
                                                                       section_name)         
         if worker.config.has_key("BASE_DIR", section_name):
            section_dico["BASE_DIR"] = worker.config.get_as_string("BASE_DIR", section_name)
         if worker.config.has_key("FLOAT_SIZE", section_name):
            section_dico["FLOAT_SIZE"] = worker.config.get_as_string("FLOAT_SIZE", section_name)         
         if worker.config.has_key("DIR_RECURSE", section_name):
            section_dico["DIR_RECURSE"] = worker.config.get_as_boolean("DIR_RECURSE", section_name)
         if worker.config.has_key("HDU_NO", section_name):
            section_dico["HDU_NO"] = worker.config.get_as_int("HDU_NO", section_name)
         
      return section_dico   

   # -----------------------------------------------------------------------------------------------
   def get_catalog_property_dico(self, section_name, job, worker):   

      section_dico = {"DIR_INPUT_LIST":[], "DIR_RECURSE":True, \
                      "BASE_DIR":".", "FILE_PATTERNS":[], \
                      "HDU_NO":1, "IS_SEXTRACTOR_FORMAT":False\
                     }
      if worker.config.has_section(section_name):
         if worker.config.has_key("FILE_PATTERNS", section_name):
            section_dico["FILE_PATTERNS"] = worker.config.get_as_list("FILE_PATTERNS",
                                                                      section_name)         
         if worker.config.has_key("DIR_INPUT_LIST", section_name):
            section_dico["DIR_INPUT_LIST"] = worker.config.get_as_list("DIR_INPUT_LIST",
                                                                       section_name)         
         if worker.config.has_key("BASE_DIR", section_name):
            section_dico["BASE_DIR"] = worker.config.get_as_string("BASE_DIR", section_name)
         if worker.config.has_key("IS_SEXTRACTOR_FORMAT", section_name):
            section_dico["IS_SEXTRACTOR_FORMAT"] = worker.config.get_as_boolean(
                                                            "IS_SEXTRACTOR_FORMAT", section_name)         
         if worker.config.has_key("DIR_RECURSE", section_name):
            section_dico["DIR_RECURSE"] = worker.config.get_as_boolean("DIR_RECURSE", section_name)
         if worker.config.has_key("HDU_NO", section_name):
            section_dico["HDU_NO"] = worker.config.get_as_int("HDU_NO", section_name)

      return section_dico   


#   # -----------------------------------------------------------------------------------------------
#   def extract_stamp_around_centroid_even(self, xc, yc, half, image):
#      """! 
#         Cut out a square stamp around the specified centroid (@c xc, @c yc).
#         @param xc x coordinate of object centroid in postage stamp 
#         @param yc y coordinate of object centroid in postage stamp
#         @param half half of the postage stamp size (of a side)
#         @param image the field image from where the postage stamp has to be cut out

#         @note stamp size is supposed to be a even number of pixels.
#      """
#      return image[xc-half+1:xc+half+1, yc-half+1:yc+half+1]

#   # -----------------------------------------------------------------------------------------------
#   def extract_stamp_around_centroid_odd(self, xc, yc, half, image):
#      """! 
#         Cut out a square stamp around the specified centroid (@c xc, @c yc).
#         @param xc x coordinate of object centroid in postage stamp 
#         @param yc y coordinate of object centroid in postage stamp
#         @param half half of the postage stamp size (of a side)
#         @param image the field image from where the postage stamp has to be cut out

#         @note stamp size is supposed to be a even number of pixels.
#      """
#      return image[xc-half+1:xc+half+1, yc-half+1:yc+half+1]


#    # -----------------------------------------------------------------------------------------------
#    def extract_stamps(self, image_data, nb_stamps, stamp_size, worker.result_output_dir, subdir):  

#   def _extract_stamps(self, image_filepath, img_no, epoch, image_file_type, 
#                             stamp_size, output_dir, float_type=32, center_star=False):
#      
#      image_data = pyfits.getdata(image_filepath).astype(float_type)
#      (image_width, image_height) = image_data.shape
#      file_main, file_ext = os.path.splitext(image_file_type) 
#      stamp_filename_pattern = "{0}_stamp_{1:03d}-{2:1d}_{3}_{4}.fits"
#      ref_x_pixel = ref_y_pixel = (stamp_size-1)/2.0
# 
#      for x in xrange(0, image_width, stamp_size):
#         for y in xrange(0, image_height, stamp_size):
#            stamp = image_data[y:y+stamp_size, x:x+stamp_size].astype(float_type) # (row, col) np.
#            if center_star:
#               (yc, xc) = scipy.ndimage.center_of_mass(stamp)
#               corr = (ref_y_pixel - yc, ref_x_pixel - xc)

#               # TEMP 
#               corr = (ref_y_pixel - yc, ref_x_pixel - xc)

#               stamp = scipy.ndimage.interpolation.shift(stamp, corr)
#               (nyc, nxc) = scipy.ndimage.center_of_mass(stamp)
#               #print (x,y), (ref_x_pixel, ref_y_pixel), (xc, yc), "=>", (nxc, nyc), "corr:", corr

#            stamp_filename = stamp_filename_pattern.format(file_main, img_no, epoch, x, y)
#            self.write_as_fits(stamp, os.path.join(output_dir, stamp_filename))

   # -----------------------------------------------------------------------------------------------
   def write_as_fits(self, data, output_filepath, header=None):
      """! 
         Write a two-dimensional data numpy array as a .fits file to some given path
         @param data data to write
         @param output_filepath full path of the -FITS file
         @param header optional FITS header data 
      """

      if os.path.exists(output_filepath):
         os.remove(output_filepath)

      if data is not None and type(data) == numpy.ndarray:
         pyfits.writeto(output_filepath, data, header=header) 

   # -----------------------------------------------------------------------------------------------
   def _plot_stamps(self, plotter, image_filepath, img_no, epoch,
                                   image_file_type, stamp_size, output_dir):
      
      image_data = pyfits.getdata(image_filepath).astype("float32")
      (image_width, image_height) = image_data.shape
      file_main, file_ext = os.path.splitext(image_file_type) 
      stamp_filename_pattern = "{0}_stamp_{1:03d}-{2:1d}_{3}_{4}.png"
      stamp_filename_pattern_3D = "{0}_stamp_3D_{1:03d}-{2:1d}_{3}_{4}.png"
      plot_filename_pattern = "{0} PSF {1:03d}-{2:1d} at {3}_{4} (log scale)"
      plot_filename_pattern_3D = "{0} PSF 3D {1:03d}-{2:1d} at {3}_{4} (log scale)"
      for x in xrange(0, image_width, stamp_size):
         for y in xrange(0, image_height, stamp_size):
            stamp = image_data[y:y + stamp_size, x:x + stamp_size]  # (row, col) numpy format
            min_pix_value = numpy.min(stamp)
            stamp += 2 * math.fabs(min_pix_value)  # stamp values are slghtly inexact
            stamp_filename = stamp_filename_pattern.format(file_main, img_no, epoch, x, y)
            plot_title = plot_filename_pattern.format(file_main, img_no, epoch, x, y)
            
            plotter.plot_stamp(numpy.log10(stamp), plot_title=plot_title,
                                      cmap="jet", color_bar=True,
                                      output_dir=output_dir, output_file=stamp_filename,
                                      show=False, logger=None)

            stamp_filename = stamp_filename_pattern_3D.format(file_main, img_no, epoch, x, y)
            plot_title = plot_filename_pattern_3D.format(file_main, img_no, epoch, x, y)

            plotter.plot_stamp_3D(numpy.log10(stamp), plot_title=plot_title,
                                      cmap="jet", rstride=2, cstride=2,
                                      output_dir=output_dir, output_file=stamp_filename,
                                      show=False, logger=None)

#   # -----------------------------------------------------------------------------------------------
#   def get_nb_objects(self, image_filepath, stamp_size):
#      """! Find the actual number of objects (i.e. postage stamps) in the image """
#      return int(len(pyfits.getdata(image_filepath)) / stamp_size)**2
   
   # -----------------------------------------------------------------------------------------------
   def get_nb_objects(self, image_filepath, stamp_size):
      """! Find the actual number of objects (i.e. postage stamps) in the image """

      header = pyfits.getheader(image_filepath)
      nb_pixels = header.get("NAXIS1", 0) * header.get("NAXIS2", 0)
      return int(nb_pixels / stamp_size ** 2)   

   # -----------------------------------------------------------------------------------------------
   def remove_failed_fits_from_result_dico(self, result_dico, master):

      failed_marking_value = master.config.get_as_float("FAILED_ELLIPTICITY_VALUE",
                                                        "SHAPE_MEASUREMENT")

      if failed_marking_value in result_dico["e1"] or failed_marking_value in result_dico["e2"]:

         good_indice = numpy.logical_and(numpy.asarray(result_dico["e1"]) != failed_marking_value,
                                         numpy.asarray(result_dico["e2"]) != failed_marking_value)  

         for var in result_dico.keys():
            var_array = numpy.asarray(result_dico[var])   
            result_dico[var] = list(var_array[good_indice])

      return result_dico   
               
   # -----------------------------------------------------------------------------------------------
   def remove_nans_from_result_dico(self, result_dico, master):


      if "sum_norm_residuals" in result_dico:
         good_indice = numpy.logical_not(numpy.isnan(result_dico["sum_norm_residuals"]))

         for var in result_dico.keys():
            var_array = numpy.asarray(result_dico[var])   
            result_dico[var] = list(var_array[good_indice])

      return result_dico   


   # -----------------------------------------------------------------------------------------------
   def save_from_list_dico(self, data_dico, output_directory, output_filename, col_list=[],
                                 key_index_map={}, key_fmt_map={}, default_fmt="%.18e"):
      """! 
         Save a dictionary to disk as a catalog file. It is assumed a list of values is attached to
         each first-level key in the dictionary (a "list" dictionary)
         @param data_dico dictionary with the data
         @param output_directory directory where to create the file
         @param output_filename name of the file to create
         @param col_list list of column names. If empty, take all the keys of the dictionary 
         @param key_index_map if not empty, contains a map with the preferred order for some keys
         @param key_fmt_map if not empty, contains the preferred output format of some key values
         @param default_fmt default format if not explicitly specified
      """

      output_catalog = None  # output catalog

      try:

         # print data_dico

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

         # --- Sort column names according to their index
         if len(key_index_map) > 0:
            sorted_index_tuples = sorted(key_index_map.items(), key=itemgetter(1, 0))
            sorted_col_names = [ sorted_index_tuples[i][0] for (k, i) in sorted_index_tuples ]
            left_col_names = [c for c in col_names if not c in sorted_col_names]
            col_names = sorted_col_names
            col_names.extend(left_col_names)

         # --- Check the column names are indeed in the dictionary
         col_names_to_remove = []
         for col_name in col_names:
            if not col_name in data_dico:
               # self.print_warning( "column: {0} not found in the dictionary".format(col_name) )
               col_names_to_remove.append(col_name)

         for col_name in col_names_to_remove:
            # self.print_warning( "Removing column: {0}".format(col_name) )
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

            # ## print col_name, len(data_dico[col_name])
            col_data_list[col_no] = numpy.concatenate(
                                                    (col_data_list[col_no], data_dico[col_name])) 

         data_matrix = numpy.asmatrix(col_data_list).transpose().squeeze()

         numpy.savetxt(output_filepath, data_matrix, fmt=cat_col_fmt, header=cat_col_comments)

      except Exception as detail:
         self.print_error("could not create catalog from dictionary ({0})".format(detail))
   

   # -----------------------------------------------------------------------------------------------
   def create_from_list_dico(self, data_dico, output_directory, output_filename, \
                                   job, worker, col_list=[],
                                   key_index_map={}, key_fmt_map={}, default_fmt="%.6e",
                                   is_ascii=True, is_sextractor=False, hdu_no=1):
      """! 
         Save a dictionary to disk as a catalog file. It is assumed a list of values is attached to
         each first-level key in the dictionary (a "list" dictionary)
         @param data_dico dictionary with the data
         @param output_directory directory where to create the file
         @param output_filename name of the file to create
         @param job an object of class MksJob to process
         @param worker instance of the worker process  
         @param col_list list of column names. If empty, take all the keys of the dictionary 
         @param key_index_map if not empty, contains a map with the preferred order for some keys
         @param key_fmt_map if not empty, contains the preferred output format of some key values
         @param default_fmt default format if not explicitly specified
         @param is_sextractor if set to True and catalog has text format, create SExtractor file
                              with ASCII_HEAD format
      """

      try:

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

         # --- Sort column names according to their index
         if len(key_index_map) > 0:
            sorted_index_tuples = sorted(key_index_map.items(), key=itemgetter(1, 0))
            sorted_col_names = [ sorted_index_tuples[i][0] for (k, i) in sorted_index_tuples ]
            left_col_names = [c for c in col_names if not c in sorted_col_names]
            col_names = sorted_col_names
            col_names.extend(left_col_names)

         # --- Check the column names are indeed in the dictionary
         col_names_to_remove = []
         for col_name in col_names:
            if not col_name in data_dico:
               # self.print_warning( "column: {0} not found in the dictionary".format(col_name) )
               col_names_to_remove.append(col_name)

         for col_name in col_names_to_remove:
            # self.print_warning( "Removing column: {0}".format(col_name) )
            col_names.remove(col_name)

         for col_name in col_names:
            col_no = col_names.index(col_name)
            # ##print col_name, len(data_dico[col_name])

            if col_name in key_index_map.keys():
               cat_col_list.insert(key_index_map[col_name], col_name)

               if col_name in key_fmt_map:
                  cat_col_fmt.insert(key_index_map[col_name], key_fmt_map[col_name])
               else:
                  cat_col_fmt.insert(key_index_map[col_name], default_fmt)
            else:
               cat_col_list.append(col_name)
               cat_col_fmt.append(default_fmt)

         # print "col_list:", col_list

         # --- Insert the columns in the catalog
         col_data_list = []   
         for col_name in cat_col_list:
            # print col_name
            cat_col_comments += col_name + " " 
            col_data_list.append(numpy.asarray([]))

         for col_name in cat_col_list:
            # print col_name
            col_no = cat_col_list.index(col_name)
            col_data_list[col_no] = numpy.concatenate(
                                                    (col_data_list[col_no], data_dico[col_name])) 
         data_matrix = numpy.asmatrix(col_data_list).transpose().squeeze()

         # print "matrix:", data_matrix

         # --- Create catalog file of the correct format
         catalog = None
         if not is_ascii:
            # --- Assume FITS format
            catalog = FITSCatalog(output_filepath, hdu_no=hdu_no)
            catalog.create_from_numpy(data_matrix, cat_col_list, ext_name=None)
         else:   
            if is_sextractor:
               # --- ASCII_HEAD SExtractor format
               catalog = SExCatalog(output_filepath)
               catalog.create_from_numpy(data_matrix, cat_col_list, col_formats=cat_col_fmt)
            else:         
               # --- ASCII tabulated format
#                catalog = TextCatalog(output_filepath)
#                catalog.create_from_numpy(data_matrix, cat_col_list)
               numpy.savetxt(output_filepath, data_matrix,
                                              fmt=cat_col_fmt, header=cat_col_comments)


      except:
         self.print_error("could not create catalog from dictionary ({0})".format(sys.exc_info()))

      finally:
         if not catalog is None:
            catalog.close()   
            

#   # -----------------------------------------------------------------------------------------------
#   def RMSD(self, stamp):
#      return numpy.sqrt()

   # -----------------------------------------------------------------------------------------------
   def RMSD(self, stamp):
      return numpy.sqrt(numpy.mean(stamp))  

   # -----------------------------------------------------------------------------------------------
   def record_data_items(self, data_dico, qvars, qvals):

      qoper_tuples = [("sum", numpy.sum), ("min", numpy.min), \
                      ("max", numpy.max)] 
#      qoper_tuples = [("sum", numpy.sum), ("min", numpy.min),\
#                      ("max", numpy.max), ("mean", numpy.mean), ("rmsd", self.RMSD)] 

      for (qlabel, qoper) in qoper_tuples:
         for (qvar, qval) in zip(qvars, qvals):
            full_var = qlabel + "_" + qvar
            if not full_var in data_dico:
               data_dico[full_var] = []
            data_dico[full_var].append(qoper(qval))


# --------------------------------------------------------------------------------------------------
class OrderedDict(dict):
    'Dictionary that remembers insertion order'
    # An inherited dict maps keys to values.
    # The inherited dict provides __getitem__, __len__, __contains__, and get.
    # The remaining methods are order-aware.
    # Big-O running times for all methods are the same as for regular dictionaries.

    # The internal self.__map dictionary maps keys to links in a doubly linked list.
    # The circular doubly linked list starts and ends with a sentinel element.
    # The sentinel element never gets deleted (this simplifies the algorithm).
    # Each link is stored as a list of length three:  [PREV, NEXT, KEY].

    def __init__(self, *args, **kwds):
        '''Initialize an ordered dictionary.  Signature is the same as for
        regular dictionaries, but keyword arguments are not recommended
        because their insertion order is arbitrary.

        '''
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__root
        except AttributeError:
            self.__root = root = []  # sentinel node
            root[:] = [root, root, None]
            self.__map = {}
        self.__update(*args, **kwds)

    def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
        'od.__setitem__(i, y) <==> od[i]=y'
        # Setting a new item creates a new link which goes at the end of the linked
        # list, and the inherited dictionary is updated with the new key/value pair.
        if key not in self:
            root = self.__root
            last = root[0]
            last[1] = root[0] = self.__map[key] = [last, root, key]
        dict_setitem(self, key, value)

    def __delitem__(self, key, dict_delitem=dict.__delitem__):
        'od.__delitem__(y) <==> del od[y]'
        # Deleting an existing item uses self.__map to find the link which is
        # then removed by updating the links in the predecessor and successor nodes.
        dict_delitem(self, key)
        link_prev, link_next, key = self.__map.pop(key)
        link_prev[1] = link_next
        link_next[0] = link_prev

    def __iter__(self):
        'od.__iter__() <==> iter(od)'
        root = self.__root
        curr = root[1]
        while curr is not root:
            yield curr[2]
            curr = curr[1]

    def __reversed__(self):
        'od.__reversed__() <==> reversed(od)'
        root = self.__root
        curr = root[0]
        while curr is not root:
            yield curr[2]
            curr = curr[0]

    def clear(self):
        'od.clear() -> None.  Remove all items from od.'
        try:
            for node in self.__map.itervalues():
                del node[:]
            root = self.__root
            root[:] = [root, root, None]
            self.__map.clear()
        except AttributeError:
            pass
        dict.clear(self)

    def popitem(self, last=True):
        '''od.popitem() -> (k, v), return and remove a (key, value) pair.
        Pairs are returned in LIFO order if last is true or FIFO order if false.

        '''
        if not self:
            raise KeyError('dictionary is empty')
        root = self.__root
        if last:
            link = root[0]
            link_prev = link[0]
            link_prev[1] = root
            root[0] = link_prev
        else:
            link = root[1]
            link_next = link[1]
            root[1] = link_next
            link_next[0] = root
        key = link[2]
        del self.__map[key]
        value = dict.pop(self, key)
        return key, value

    # -- the following methods do not depend on the internal structure --

    def keys(self):
        'od.keys() -> list of keys in od'
        return list(self)

    def values(self):
        'od.values() -> list of values in od'
        return [self[key] for key in self]

    def items(self):
        'od.items() -> list of (key, value) pairs in od'
        return [(key, self[key]) for key in self]

    def iterkeys(self):
        'od.iterkeys() -> an iterator over the keys in od'
        return iter(self)

    def itervalues(self):
        'od.itervalues -> an iterator over the values in od'
        for k in self:
            yield self[k]

    def iteritems(self):
        'od.iteritems -> an iterator over the (key, value) items in od'
        for k in self:
            yield (k, self[k])

    def update(*args, **kwds):
        '''od.update(E, **F) -> None.  Update od from dict/iterable E and F.

        If E is a dict instance, does:           for k in E: od[k] = E[k]
        If E has a .keys() method, does:         for k in E.keys(): od[k] = E[k]
        Or if E is an iterable of items, does:   for k, v in E: od[k] = v
        In either case, this is followed by:     for k, v in F.items(): od[k] = v

        '''
        if len(args) > 2:
            raise TypeError('update() takes at most 2 positional '
                            'arguments (%d given)' % (len(args),))
        elif not args:
            raise TypeError('update() takes at least 1 argument (0 given)')
        self = args[0]
        # Make progressively weaker assumptions about "other"
        other = ()
        if len(args) == 2:
            other = args[1]
        if isinstance(other, dict):
            for key in other:
                self[key] = other[key]
        elif hasattr(other, 'keys'):
            for key in other.keys():
                self[key] = other[key]
        else:
            for key, value in other:
                self[key] = value
        for key, value in kwds.items():
            self[key] = value

    __update = update  # let subclasses override update without breaking __init__

    __marker = object()

    def pop(self, key, default=__marker):
        '''od.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised.

        '''
        if key in self:
            result = self[key]
            del self[key]
            return result
        if default is self.__marker:
            raise KeyError(key)
        return default

    def setdefault(self, key, default=None):
        'od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od'
        if key in self:
            return self[key]
        self[key] = default
        return default

    def __repr__(self, _repr_running={}):
        'od.__repr__() <==> repr(od)'
        call_key = id(self), _get_ident()
        if call_key in _repr_running:
            return '...'
        _repr_running[call_key] = 1
        try:
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())
        finally:
            del _repr_running[call_key]

    def __reduce__(self):
        'Return state information for pickling'
        items = [[k, self[k]] for k in self]
        inst_dict = vars(self).copy()
        for k in vars(OrderedDict()):
            inst_dict.pop(k, None)
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def copy(self):
        'od.copy() -> a shallow copy of od'
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        '''OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S
        and values equal to v (which defaults to None).

        '''
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        '''od.__eq__(y) <==> od==y.  Comparison to another OD is order-sensitive
        while comparison to a regular mapping is order-insensitive.

        '''
        if isinstance(other, OrderedDict):
            return len(self) == len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    # -- the following methods are only used in Python 2.7 --

    def viewkeys(self):
        "od.viewkeys() -> a set-like object providing a view on od's keys"
        return KeysView(self)

    def viewvalues(self):
        "od.viewvalues() -> an object providing a view on od's values"
        return ValuesView(self)

    def viewitems(self):
        "od.viewitems() -> a set-like object providing a view on od's items"
        return ItemsView(self)

# -- EOF gfit_helper.py

