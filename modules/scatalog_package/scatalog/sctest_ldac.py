#!/usr/bin/env python

"""! 
   @package scatalog.sctest Testing program for SCatalog
   @author Marc Gentile
   @file sctest_ldac.py
   Testing program for SCatalog
"""

import numpy as np
import numpy.random as random
import shutil

from scatalog import *

def test():

   shutil.copy("test/galaxy_catalog-001_orig.fits", "test/galaxy_catalog-001.fits")

   # --- Opening  a .FITS catalog
   fits_cat = LDACFITSCatalog("test/se_gal_image_LDAC.fits", frame_no=1)
   fits_cat.open()   
   
   print(fits_cat)
   print "Nb rows:", fits_cat.get_nb_rows(), "Nb cols:",  fits_cat.get_nb_cols()
   print "Col names:", fits_cat.get_col_names()
   print "Col exists FLUX_BEST:", fits_cat.col_exists("FLUX_BESt")    
   print "Col index('FLAGS')", fits_cat.get_col_index("FLAGS")
   print "Col formats:", fits_cat.get_col_formats()
   print "Col comments:", fits_cat.get_col_comments()
   print "Col data('FWHM_IMAGE')", fits_cat.get_named_col_data('FWHM_IMAGE')
   print "Col data(0)", fits_cat.get_col_data(0)
   print "get_data()", fits_cat.get_data()
   print "header:", fits_cat.get_header()
   fits_cat.flush()  
   fits_cat.close()   

 
   # --- Updating a .FITS LDAC Catalog
   print "*** UPDATE TEST ***"
   fits_cat = LDACFITSCatalog("test/se_gal_image_LDAC.fits", 
                              open_mode=LDACFITSCatalog.OpenMode.ReadWrite, frame_no=1)
   fits_cat.open()
   
   print fits_cat.get_col_names()
   print fits_cat.get_nb_rows()
   
   col_1_data = np.ones(10000)
   col_2_data = np.ones(10000)
 
   shutil.copy("test/se_gal_image_LDAC_orig.fits", "test/se_gal_image_LDAC.fits.fits")
 
   fits_cat.add_col(col_name="Col 1", col_comment="New column 1", col_data=col_1_data)
   fits_cat.add_col(col_name="Col 2", col_comment="New column 2", col_data=col_2_data)
  
   fits_cat.close()

   # --- Creating a FITS LDAC catalog
   fits_cat = LDACFITSCatalog("test/se_gal_image_LDAC.fits")
   fits_cat.open()     
      
   test_header = fits_cat.get_header()
   test_data = fits_cat.get_data()
   test_col_names = fits_cat.get_col_names()
   
   print "test_header:",test_header
   print "test_data:", test_data, test_data.shape
   print "test_col_names:", fits_cat.get_col_names()
   
   
   if os.path.exists("test/empty_np_ldac_cat.fits"):
      os.remove("test/empty_np_ldac_cat.fits")
      
   np_fits_cat = LDACFITSCatalog("test/empty_np_ldac_cat.fits")
   np_fits_cat.create()
   np_fits_cat.open()
   
   print np_fits_cat
   np_fits_cat.close()
 
   # --- Creating a FITS catalog from a numpy array
   print "Creating FITS catalog from a numpy array..."
   matrix = random.rand(100, 5)   
   np_fits_cat = LDACFITSCatalog("test/np_ldac_cat.fits")
   col_names = ["new_col 1", "new_col 2", "new_col 3", "new_col 4", "new_col 5"]
#    np_fits_cat.create_from_numpy(matrix, col_names, header=None, ext_name="Objects")

   np_fits_cat.create_from_numpy(test_data, col_names, header=test_header, ext_name="my objects")

   print "Nb rows:", np_fits_cat.get_nb_rows(), "Nb cols:",  np_fits_cat.get_nb_cols()
   print "Col names:", np_fits_cat.get_col_names()    
   np_fits_cat.close()
   

if __name__ == "__main__":
   test()
