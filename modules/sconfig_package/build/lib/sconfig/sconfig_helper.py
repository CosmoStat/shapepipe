"""! 
   @package sconfig.sconfig_helper Helper functions for sconfig
   @author Marc Gentile
   @file sconfig_helper.py
   Helper functions for sconfig
"""

import os, sys
import string


def filepath_exists(filepath):
   """! 
      Safe way of checking if a file (not a directory) does exist. 
      @param filepath file path to chech for existence
      @retval True if file exists
      @retval False otherwise
   """
   try:
      with open(filepath):  
         return True
   except IOError:
      return False

def expand_env_vars(var_string):
   """! 
      Locate and expand environment variables in a string 
      @param var_string Python string in which environment variables should be expanded 
      @return Python string with expanded environment variables
   """

   new_var_string = var_string
   try:
      new_var_string = string.Template(var_string).substitute(os.environ)
   except (KeyError, ValueError):
      print "Error substituting environment variables in: '%s' (%s)" %(var_string, sys.exc_info()[1])
      
   return new_var_string
