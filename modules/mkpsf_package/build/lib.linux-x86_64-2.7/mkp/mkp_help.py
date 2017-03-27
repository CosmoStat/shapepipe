"""! 
   mkp_helper.py: various helper functions
""" 

# -- Python imports
import os, sys
import imp
import math
import numpy
import scipy.ndimage
import math
import pyfits

# -- External import
from scatalog import *              # catalog management
from mpfg.mp_helper import *         # base Helper

# -- Module-specific imports  



# -------------------------------------------------------------------------------------------------
class MkpHelper(Helper):

   """! Convenient utility functions """

   # -----------------------------------------------------------------------------------------------
   def __init(self):
      Helper.__init__(self)

   # -----------------------------------------------------------------------------------------------
   def extract_stamp_around_centroid(self, xc, yc, half, image):
      """! 
         Cut out a square stamp around the specified centroid (@c xc, @c yc).
         @param xc x coordinate of object centroid in postage stamp 
         @param yc y coordinate of object centroid in postage stamp
         @param half half of the postage stamp size (of a side)
         @param image the field image from where the postage stamp has to be cut out

         @note stamp size is supposed to be a even number of pixels.
      """

      return image[xc-half+1:xc+half+1, yc-half+1:yc+half+1]


   # -----------------------------------------------------------------------------------------------
   def _extract_stamps(self, image_filepath, img_no, epoch, image_file_type, 
                             stamp_size, output_dir, float_type=32, center_star=False):
      
      image_data = pyfits.getdata(image_filepath).astype(float_type)
      (image_width, image_height) = image_data.shape
      file_main, file_ext = os.path.splitext(image_file_type) 
      stamp_filename_pattern = "{0}_stamp_{1:03d}-{2:1d}_{3}_{4}.fits"
      ref_x_pixel = ref_y_pixel = (stamp_size-1)/2.0
 
      for x in xrange(0, image_width, stamp_size):
         for y in xrange(0, image_height, stamp_size):
            stamp = image_data[y:y+stamp_size, x:x+stamp_size].astype(float_type) # (row, col) np.
            if center_star:
               (yc, xc) = scipy.ndimage.center_of_mass(stamp)
               corr = (ref_y_pixel - yc, ref_x_pixel - xc)

               # TEMP 
               corr = (ref_y_pixel - yc, ref_x_pixel - xc)

               stamp = scipy.ndimage.interpolation.shift(stamp, corr)
               #(nyc, nxc) = scipy.ndimage.center_of_mass(stamp)
               #print (x,y), (ref_x_pixel, ref_y_pixel), (xc, yc), "=>", (nxc, nyc), "corr:", corr

            stamp_filename = stamp_filename_pattern.format(file_main, img_no, epoch, x, y)
            self.write_as_fits(stamp, os.path.join(output_dir, stamp_filename))

   # -----------------------------------------------------------------------------------------------
   def write_as_fits(self, data, output_filepath, header=None):
      """! 
         Write a two-dimensional data numpy array as a .fits file to some given path
         @param data data to write
         @param output_filepath full path of the -FITS file
      """

      if os.path.exists(output_filepath):
         os.remove(output_filepath)

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
            stamp = image_data[y:y+stamp_size, x:x+stamp_size] # (row, col) numpy format
            min_pix_value = numpy.min(stamp)
            stamp += 2*math.fabs(min_pix_value) # stamp values are slghtly inexact
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

   # -----------------------------------------------------------------------------------------------
   def get_nb_objects(self, image_filepath, stamp_size):
      """! Find the actual number of objects (i.e. postage stamps) in the image """
      return int(len(pyfits.getdata(image_filepath)) / stamp_size)**2

   # -----------------------------------------------------------------------------------------------
   def mark_fits_stamps_file(self, fits_path, coordinates, marking_value, stamp_size):
      """! Mark splitted or flagged objects in the check .fits files """

      check_image = pyfits.getdata(fits_path)   # SE-generated chech .fits image as numpy array
      for (row, col) in coordinates:
         stamp = check_image[row: row + stamp_size, col: col + stamp_size]
         stamp[:1] = marking_value
         stamp[stamp_size:stamp_size] = marking_value
         stamp[0:,0:1] = marking_value   
         stamp[0:,stamp_size:stamp_size] = marking_value

      self.write_as_fits(check_image, fits_path)

   # -----------------------------------------------------------------------------------------------   
   def mark_centroids(self, fits_path, centroids, marking_value):
      """! Mark centroid locations in the check .fits files """

      check_image = pyfits.getdata(fits_path)   # SE-generated chech .fits image as numpy array
      for (row, col) in centroids:
         check_image[math.floor(row + 0.5), math.floor(col + 0.5)] = marking_value        

      self.write_as_fits(check_image, fits_path)

   # -----------------------------------------------------------------------------------------------
   def make_stats(self, job_result, master):
   
      job = job_result.job                # associated job        
      object_per_type_dico = job_result.result

      file_types = object_per_type_dico.keys()   # all results for each image types

      for file_type in file_types:

         object_dico = object_per_type_dico[file_type]

         # --- Where plots will be stored
         branch_tree = job.get_branch_tree()
         stats_pathname = os.path.join(master.stat_output_dir, branch_tree)    

         if master.logging_enabled():
            master.logger.log_info_p(
              "{0} - Img {1:03d} - /{2} - Generating statistics for {3}...".format(
                                               master, job.img_no, branch_tree, file_type))
         # --- Stats catalog
         stats_filename = "stats_{0}_img {1:03}-{2:1d}_{3}.txt".format(
                          branch_tree, job.img_no, job.epoch, file_type)
         stats_filepath = os.path.join(stats_pathname, stats_filename)

         catalog = TextCatalog(stats_filepath)
         catalog.create()      

         # --- Produce stats catalogs...
         fd = open(stats_filepath, "w+")

         fd.write('\n')
         title = "STATISTICS FOR: Image {0:03d}-{1:1d} - /{2} - {3}\n".format(
                                 job.img_no, job.epoch, branch_tree, file_type)
         fd.write(title)
         fd.write("=" * len(title) + "\n\n")

         for col_name in job_result.stats_dico[file_type].keys():

            col_stats_dico = job_result.stats_dico[file_type][col_name]

            fd.write('--- Statistics on {0} ---\n'.format(col_name))
            oper_keys = col_stats_dico.keys()
            fd.write('{0}:\t'.format(col_name))
            for oper in oper_keys:
               fd.write('{0}={1:.6f}\t'.format(oper, col_stats_dico[oper]))
            fd.write('\n\n')
      
         fd.close()

   # -----------------------------------------------------------------------------------------------
   def import_module(self, module_name, module_dir):
      """! Find and load a module <module_name> from a directory <module_dir> """

      try:

         file_obj, filename, data = imp.find_module(module_name, [module_dir])
         imp.load_module(module_name, file_obj, filename, data)
         module = imp.load_module(module_name, file_obj, filename, data)
         return module

      except:
         self.print_error("Could not load module {0}.".format(module_name))

   # -----------------------------------------------------------------------------------------------
   def load_method_module(self, method_name, method_dir):
      """! Load the Python module to access the method's code """   

      return self._import_module(method_name, method_dir)



# -- EOF mkp_help.py
