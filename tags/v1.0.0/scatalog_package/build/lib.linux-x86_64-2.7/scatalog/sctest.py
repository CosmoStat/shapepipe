#!/usr/bin/env python

"""! 
   @package scatalog.sctest Testing program for SCatalog
   @author Marc Gentile
   @file sctest.py
   Testing program for SCatalog
"""

import sys
import numpy as np
import numpy.random as random

from scatalog import *

def test():

   # --- Opening  a .FITS catalog
   fits_cat = FITSCatalog("test/galaxy_catalog-001.fits", hdu_no=1, mem_map=True)
   fits_cat.open()   
   print(fits_cat)
   print "Nb rows:", fits_cat.get_nb_rows(), "Nb cols:",  fits_cat.get_nb_cols()
   print "Col names:", fits_cat.get_col_names()    
   print "Col index('ID')", fits_cat.get_col_index("ID")   
   print "Col data("'ID'")", fits_cat.get_col_data("ID")
   print "Col data(0)", fits_cat.get_col_data(0)
   print "get_data()", fits_cat.get_data()
   fits_cat.flush()  
   fits_cat.close()   

   # --- Creating a Text catalog
   new_txt_cat = TextCatalog("test/new_txt_cat.txt", null_char="*")
   new_txt_cat.create()

   col_data = list(np.arange(0,1000))   
   col_1 = TextCatalog.Column("col_1", format="% 9.6f", comment="col comment", data=col_data)
   print "col 1: name:", col_1.name, "Type:", col_1.get_type(), "Format:", col_1.format,\
                                                                "Comment:", col_1.comment
   print "Col 1: len(data)", len(col_1.data)
   print "Nb rows:", col_1.get_nb_rows()   
   print "Col 1: Info:", col_1.get_info()
   new_txt_cat.append_col(col_1)

   col_2 = TextCatalog.Column("col_2", format="% 9.6f", comment="col comment", data=col_data)
   print "Col 2: len(data)", len(col_2.data)
   print "Nb rows:", col_2.get_nb_rows()   
   print "Col 2: Info:", col_2.get_info()
   new_txt_cat.append_col(col_2)

   print "Filename:", new_txt_cat.filename, "nb rows:", new_txt_cat.get_nb_rows(),\
                                            "nb_cols:",new_txt_cat.get_nb_cols()
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
   col_6 = SExCatalog.Column("col_6", format="% 9.6f", comment="col comment", data=col_data)
   print "col 6: name:", col_6.name, "Type:", col_6.get_type(), "Format:", col_6.format,\
                                                                "Comment:", col_6.comment
   print "Col 6: len(data)", len(col_6.data)
   print "Nb rows:", col_6.get_nb_rows()   
   print "Col 6: Info:", col_6.get_info()
   np_SE_cat.append_col(col_6)

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
   new_col = SExCatalog.Column("Extra_col", format="% 5.4f", comment="extra col", data=col_data)
   se_cat.append_col(new_col)

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
   new_col = SExCatalog.Column("Extra_col", format="% 9.6f", comment="extra col", data=col_data)
   se_cat_2.append_col(new_col)

   print "Filename:", se_cat_2.filename, "nb rows:", se_cat_2.get_nb_rows(),\
                                         "nb_cols:", se_cat_2.get_nb_cols()
   print "Top Header:", se_cat_2.get_header(), "headers:", se_cat_2.get_headers()
   print "column names:", se_cat_2.get_col_names()
   print "info:", se_cat_2.get_info()   

   # Settng column data
   print "Setting column data..."   
   se_cat_2.set_named_col_data("Extra_col", np.random.rand(len(np_data)))

   # Sorting
   print "Sorting..."
   se_cat_2.sort("Extra_col", descending=True)

   print "info after col data set:", se_cat_2.get_info()   

   # Saving
   print "Saving..."   
   se_cat_2.save()   
   se_cat_2.close()

if __name__ == "__main__":
   test()
