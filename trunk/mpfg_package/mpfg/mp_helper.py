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
