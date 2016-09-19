"""! 
   @package scatalog.scatalog_helper Simple catalog management
   @author Marc Gentile
   @file scatalog_helper.py
   Helper functions for scatalog
"""

# -- Python imports
import os, sys
import numpy as np

# -------------------------------------------------------------------------------------------------
class Helper(object):

   """! Convenient utility functions that can be shared across sub-classes. """


#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Fix int64 and float64 types, not supported by asciidata
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   def fix_data_type(self, data):
#      """! Fix int64 and float64 types, not supported by asciidata! """
#      if not data is None and len(data) > 0:
#         if type(data) == list:
#            data = np.asarray(data)

#         print "***:", type(data[0])
#         if type(data[0]) == np.float64:
#            data = data.astype(np.float)
#         elif type(data[0]) == np.int64:
#            print "***: => numpy.int64"
#            data = data.astype(np.float)
#            print "***: => data type:", type(data), type(data[0])

#         print "### data type:", type(data[0])
#      return data

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Check if a file (not a directory) exists
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def file_exists(self, filepath):
      """! 
         Safe way of checking if a file (not a directory) does exist. 
         @param filepath file path to chech for existence
         @retval True if file exists
         @retval False otherwise
      """
      try:
        with open(filepath) as file:
            return True
      except IOError as e:
        return False

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Print Error, warning or information message
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def print_error(self, msg):
      """! 
         Convenience function to print error message. 
         @param msg error message to print   
      """
      print("{0} *** ERROR ***: {1}".format(self.get_main_basename(), msg))

   def print_warning(self, msg):
      """! 
         Convenience function to print warning message. 
         @param msg warning message to print   
      """
      print("{0} *** Warning ***: {1}".format(self.get_main_basename(), msg))

   def print_info(self, msg):
      """! 
         Convenience function to print information message. 
         @param msg information message to print   
      """
      print("{0} *** Info ***: {1}".format(self.get_main_basename(), msg))


# -- EOF scatalog_helper.py
