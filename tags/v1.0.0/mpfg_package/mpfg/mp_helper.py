"""!
   @package mpfg.mp_helper various helper functions
   @author Marc Gentile
   @file mp_helper.py

   mp_helper.py: various helper functions
""" 

# -- Python imports
import os, sys
import numpy

# -- External imports
from sconfig import *         # external configuration



# -------------------------------------------------------------------------------------------------
class Helper(object):

   """! Convenient utility functions that can be shared across sub-classes. """

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

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Get the base name of the main program 
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def get_main_basename(self):
      """! @return the base name of the main program """
      return os.path.basename(sys.argv[0])

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Read the external configuration file and return a SConfig object
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def read_config(self, option_dico):
      """! 
         Read the external configuration file and return a SConfig object.
         @param option_dico dictionary containing command-line options  
      """

      config_filepath = os.path.join(option_dico['-d'], option_dico['-c'])

      config = None
      if self.file_exists(config_filepath): 
         config = SConfig(config_filepath)
      else:
         self.print_error("Error reading external configuration file {0}".format(config_filepath))      
      return config

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Dump configuration data to a logger
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def dump_config(self, config, logger, label=""):
      """!
         Dump content of configuration file to the log
         @param config a SConfig object
         @param logger a FileLogger object
         @param label an optional text to insert as header & footer
      """
      logger.log_info("========= {0}BEGIN CONFIGURATION =========".format(label))
      config.dump_to_file(logger.filedescr)
      logger.log_info("========= {0}END CONFIGURATION  =========\n".format(label))
      logger.flush()

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Check if logging is enabled
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def is_logging(self, config):
      """!
         Tell whether logging is enabled
         @retval True enabled
         @retval False not enabled
      """
      return config.get_as_boolean("MASTER_LOGGING_ENABLED", "LOGGING")  

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Attempt to make a directory
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def make_dir(self, dir_name):
      """! 
         Attempt to make a directory 
         @param dir_name name of directory to create
      """  
      try:
         if not os.path.exists(dir_name):
            os.mkdir(dir_name)
            return True
      except(OSError):
         self.print_warning(sys.exc_info()[1])
         return False

#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # Save a "list" dictionary to a file
#   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   def save_from_list_dico(self, data_dico, output_directory, output_filename, col_list = [], 
#                                 key_index_map=[], key_fmt_map=[]):
#      """! 
#         Save a dictionary to disk as a catalog file. It is assumed a list of values is attached to
#         each first-level key in the dictionary (a "list" dicrionary)
#         @param data_dico dictionary with the data
#         @param output_directory directory where to create the file
#         @param output_filename name of the file to create
#         @param col_list list of column names. If empty, take all the keys of the dictionary 
#         @param key_index-map if not empty, contains a map with the preferred order for some keys
#         @param key_fmt_map if not empty, contains the preferred output format of some key values
#      """

#      output_catalog = None    # output catalog

#      try:
#         # --- Create output file
#         output_filepath = os.path.join(output_directory, output_filename)

#         # --- Build the list of columns in the required order
#         cat_col_list = []
#         cat_col_fmt = []
#         cat_col_comments = ""

#         if len(col_list) == 0:
#            col_names = data_dico.keys()
#         else:
#            col_names = col_list

#         # --- Check the column names are indeed in the dictionary
#         for col_name in col_names:
#            if not col_name in data_dico:
#               self.print_warning( "column: {0} not found in the dictionary".format(col_name) )
#               col_names.remove(col_name)

#         for col_name in col_names:
#            col_no = col_names.index(col_name)
#            if col_name in key_index_map:
#               cat_col_list.insert(key_index_map[col_name], col_name)
#               if col_name in key_fmt_map:
#                  cat_col_fmt.insert(key_index_map[col_name], key_fmt_map[col_name])
#            else:
#               cat_col_list.append(col_name)
#               cat_col_fmt.append("%.9f")

#         # --- Insert the columns in the catalog
#         col_data_list = []   
#         for col_name in cat_col_list:
#           cat_col_comments += col_name + " " 
#           col_data_list.append(numpy.asarray([]))

#         for col_name in cat_col_list:
#            col_no = cat_col_list.index(col_name)
#            col_data_list[col_no] = numpy.concatenate(  
#                                                    (col_data_list[col_no], data_dico[col_name] ) ) 

#         data_matrix = numpy.asmatrix(col_data_list).transpose().squeeze()
#         numpy.savetxt(output_filepath, data_matrix, fmt=cat_col_fmt, header=cat_col_comments)

#      except:
#         self.print_error("could not create catalog from dictionary ({0})".format(sys.exc_info()))



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


# -- EOF mp_helper.py
