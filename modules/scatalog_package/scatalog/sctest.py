#!/usr/bin/env python

"""! 
   @package scatalog.sctest Testing program for SCatalog
   @author Marc Gentile
   @file sctest.py
   Testing program for SCatalog
"""

import numpy as np
import numpy.random as random
import shutil

from scatalog import *

def test():

   shutil.copy("test/galaxy_catalog-001_orig.fits", "test/galaxy_catalog-001.fits")

   print "*** QUERY TEST ***"

   # --- Opening  a .FITS catalog
   fits_cat = FITSCatalog("test/galaxy_catalog-001.fits")
   fits_cat.open()
   
   #print "Nb rows:", fits_cat.get_nb_rows()   
   col_data = fits_cat.get_named_col_data('x')  
   print "Col data('x')", col_data

   
   #print(fits_cat)
   print "Col names:", fits_cat.get_col_names()
   #print "Nb rows:", fits_cat.get_nb_rows()
   print "Nb cols:",  fits_cat.get_nb_cols()
   print "Col exists ID:", fits_cat.col_exists("ID")    
   print "Col index('ID')", fits_cat.get_col_index("ID")
   print "Col formats:", fits_cat.get_col_formats()
   print "Comments:", fits_cat.get_comments()
   print "Col comments:", fits_cat.get_col_comments()
   print "Col data('x')", fits_cat.get_named_col_data("x")
   print "Col data(0)", fits_cat.get_col_data(0)
   print "get_data()", fits_cat.get_data()
   print "Header:", fits_cat.get_header()
   print "Headers:", fits_cat.get_headers()
   fits_cat.flush()  
   fits_cat.close()   
 
   # --- Updating a .FITS Catalog
   print "*** UPDATE TEST ***"
   fits_cat = FITSCatalog("test/galaxy_catalog-001.fits",
                           open_mode=FITSCatalog.OpenMode.ReadWrite)
   fits_cat.open()
   
   print fits_cat.get_col_names()
   #print fits_cat.get_nb_rows()
   
   col_1_data = np.ones(10000)
   col_2_data = np.ones(10000)

   fits_cat.add_col(hdu_no=1, col_name="Col 1", col_comment="New column 1", col_data=col_1_data)
   fits_cat.add_col(hdu_no=1, col_name="Col 2", col_comment="New column 2", col_data=col_2_data)
  
   fits_cat.close()   

   # --- Creating a FITS catalog
   fits_cat = FITSCatalog("test/galaxy_catalog-001.fits")
   fits_cat.open()     
      
   test_header = fits_cat.get_header()
   test_data = fits_cat.get_data()
   test_col_names = fits_cat.get_col_names()
   
   print "test_header:",test_header
   print "test_data:", test_data, test_data.shape
   print "test_col_names:", fits_cat.get_col_names()
   
   
   if os.path.exists("test/empty_np_cat.fits"):
      os.remove("test/empty_np_cat.fits")
      
   np_fits_cat = FITSCatalog("test/empty_np_cat.fits")
   np_fits_cat.create()
   np_fits_cat.open()
   print np_fits_cat
   print "Empty FITS, nb_rows=", np_fits_cat.get_nb_rows()
   np_fits_cat.close()
 
   # --- Creating a FITS catalog from a numpy array
   print "Creating FITS catalog from a numpy array..."
   matrix = random.rand(100, 5)   
   np_fits_cat = FITSCatalog("test/np_cat.fits")
   col_names = ["col 1", "col 2", "col 3", "col 4", "col 5"]
#    np_fits_cat.create_from_numpy(matrix, col_names, header=None, ext_name="Objects")

   np_fits_cat.create_from_numpy(test_data, test_col_names, header=test_header, ext_name="Objects")

   print "Nb rows:", np_fits_cat.get_nb_rows(), "Nb cols:",  np_fits_cat.get_nb_cols()
   print "Col names:", np_fits_cat.get_col_names()    
   np_fits_cat.close()

    # --- Creating a Text catalog
   new_txt_cat = TextCatalog("test/new_txt_cat.txt", null_char="*")
   new_txt_cat.create()
 
   col_data = list(np.arange(0,1000))   
   new_txt_cat.add_col("Col1", 
                       col_format="% 9.6f", col_comment="col comment", col_data=col_data)
   new_txt_cat.add_col("Col 2", 
                       col_format="% 9.6f", col_comment="col comment", col_data=col_data)
 
   print "Filename:", new_txt_cat.filename, "nb rows:", new_txt_cat.get_nb_rows(),\
                                            "nb_cols:", new_txt_cat.get_nb_cols()
   print "Null char:", new_txt_cat.null_char, "Comment char:", new_txt_cat.comment_char
   print "Top Header:", new_txt_cat.get_header(), "headers:", new_txt_cat.get_headers()
   print "column names:", new_txt_cat.get_col_names()   
   print "info:", new_txt_cat.get_info()     
 
   new_txt_cat.save()
   new_txt_cat.close()
 
   # --- Creating a Text catalog from a numpy array
   matrix = random.rand(100, 5)   
   np_txt_cat = TextCatalog("test/np_cat.txt", null_char="*")
   col_names = ["col1", "col2", "col3", "col4", "col5"]
   col_comments = [["#col1", " #col2", "#col3", "#col4", "#col5"]]      
   np_txt_cat.create_from_numpy(matrix, col_names)
 
   print "Filename:", np_txt_cat.filename, "nb rows:", np_txt_cat.get_nb_rows(),\
                                            "nb_cols:",np_txt_cat.get_nb_cols()
   print "Null char:", np_txt_cat.null_char, "Comment char:", np_txt_cat.comment_char
   print "Top Header:", np_txt_cat.get_header(), "headers:", np_txt_cat.get_headers()
   print "column names:", np_txt_cat.get_col_names()   
   print "info:", np_txt_cat.get_info()
 
   print "*** Col 1", np_txt_cat.get_named_col("col1")
   print "col 1 data:", np_txt_cat.get_named_col_data("col1")
 
   np_txt_cat.close()
 
   # --- Opening Text Catalog
   txt_cat = TextCatalog("test/epoch_dither-000-0.txt")
   txt_cat.open()
   print "Filename:", txt_cat.filename, "nb rows:", txt_cat.get_nb_rows(),\
                                        "nb_cols:", txt_cat.get_nb_cols()
   print "Top Header:", txt_cat.get_header(), "headers:", txt_cat.get_headers()
   print "column names:", txt_cat.get_col_names()   
   print "info:", txt_cat.get_info()   
    
   np_data = txt_cat.get_data()
   print "Numpy data:", np_data
 
   txt_cat.close()
 
   # --- Create an empty SE catalog from a numpy array
   matrix = []   
   np_SE_cat = SExCatalog("test/np_SE_cat.txt", null_char="*")
 
   col_names = ["col1", "col2", "col3", "col4", "col5"]
   col_comments = ["col1", " #col2", "#col3", "#col4", "#col5"] 
   col_formats = ["% 9.6f", "% 9.6f", "% 9.6f", "% 9.6f", "% 9.6f"] 
   np_SE_cat.create_from_numpy(matrix, col_names, col_comments, col_formats)
 
   np_SE_cat.save()
   np_SE_cat.close()
 
   # --- Creating a SE catalog from a numpy array
   matrix = random.rand(1000, 5)   
   np_SE_cat = SExCatalog("test/np_SE_cat.txt", null_char="*")
   col_names = ["col1", "col2", "col3", "col4", "col5"]
   col_comments = ["col1", " #col2", "#col3", "#col4", "#col5"] 
   col_formats = ["% 9.6f", "% 9.6f", "% 9.6f", "% 9.6f", "% 9.6f"] 
   np_SE_cat.create_from_numpy(matrix, col_names, col_comments, col_formats)
 
   col_data = list(np.arange(0,1000))   
   np_SE_cat.add_col("Col 6", col_format="% 9.6f", col_comment="col comment", col_data=col_data)
 
   print "Filename:", np_SE_cat.filename, "nb rows:", np_SE_cat.get_nb_rows(),\
                                          "nb_cols:",np_SE_cat.get_nb_cols()
   print "Null char:", np_SE_cat.null_char, "Comment char:", np_SE_cat.comment_char
   print "Top Header:", np_SE_cat.get_header(), "headers:", np_SE_cat.get_headers()
   print "column names:", np_SE_cat.get_col_names()   
   print "info:", np_SE_cat.get_info()      
   print "col comments:", np_SE_cat.get_col_comments()
   print "col formats:", np_SE_cat.get_col_formats()
 
   np_SE_cat.save()
   np_SE_cat.close()
 
   # --- Opening SE Catalog
   se_cat = SExCatalog("test/star_catalog-000.txt")
   se_cat.open()
   print "Filename:", se_cat.filename, "nb rows:", se_cat.get_nb_rows(),\
                                       "nb_cols:", se_cat.get_nb_cols()
   print "Top Header:", se_cat.get_header(), "headers:", se_cat.get_headers()
   print "column names:", se_cat.get_col_names()
   print "info:", se_cat.get_info()   
   np_data = se_cat.get_data()
   print "Numpy data:", np_data.shape
 
   print se_cat.get_item(0,3)
   se_cat.set_item(0,3, 4.0)
   print se_cat.get_item(0,3)
   se_cat.save()
 
   se_cat.close()
 
   # --- Opening SE Catalog again
   se_cat = SExCatalog("test/stars_000_S82m0m_predicted.cat")
   se_cat.open()
 
   np_data = se_cat.get_data()
 
   # Adding column
   print "Adding column..." 
   if se_cat.col_exists("Extra_col"):
      se_cat.remove_named_col("Extra_col")
   
   #col_data = np.arange(0,len(np_data)).astype(np.float)
   col_data = np.random.rand(len(np_data))
   se_cat.add_col("Extra col", col_format="% 5.4f", col_comment="extra col", col_data=col_data)
 
   print "info:", se_cat.get_info()  
 
   np_data = se_cat.get_data()
   print "Numpy data:", np_data
       
#   se_cat.set_data(np_data)   
 
   print "Filename:", se_cat.filename, "nb rows:", se_cat.get_nb_rows(),\
                                       "nb_cols:", se_cat.get_nb_cols()
   print "Top Header:", se_cat.get_header(), "headers:", se_cat.get_headers()
   print "column names:", se_cat.get_col_names()
   print "info:", se_cat.get_info()   
 
   se_cat.save(write_col_info=True, write_header=True)
 
   se_cat.close()
 
   # Creating and opening catalog
   print "Creating catalog..."
   se_cat_2 = SExCatalog("test/stars_000_S82m0m_predicted.cat")
   print "Opening catalog..."
   se_cat_2.open()
   print "Getting catalog data..."
   np_data = se_cat_2.get_data()
   print "Numpy data:", np_data.shape
 
   # Removing column
   print "Removing column..."
   if se_cat_2.col_exists("Extra_col"):
      se_cat_2.remove_named_col("Extra_col")
 
   print "info after col removed:", se_cat_2.get_info()   
 
   # Adding column
   print "Adding column..."   
   col_data = list(np.arange(0,len(np_data)))   
   se_cat_2.add_col("extra col", col_format="%9.6f", col_comment="extra col", col_data=col_data)
 
   print "Filename:", se_cat_2.filename, "nb rows:", se_cat_2.get_nb_rows(),\
                                         "nb_cols:", se_cat_2.get_nb_cols()
   print "Top Header:", se_cat_2.get_header(), "headers:", se_cat_2.get_headers()
   print "column names:", se_cat_2.get_col_names()
   print "info:", se_cat_2.get_info()   
 
   # Settng column data
   print "Setting column data..."   
   se_cat_2.set_named_col_data("extra col", np.random.rand(len(np_data)))
 
   # Sorting
   print "Sorting..."
   se_cat_2.sort("extra col", descending=True)
 
   print "info after col data set:", se_cat_2.get_info()   
 
   # Saving
   print "Saving..."   
   se_cat_2.save()   
   se_cat_2.close()

if __name__ == "__main__":
   test()
