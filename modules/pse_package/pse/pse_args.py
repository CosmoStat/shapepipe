"""! 
   pse_args.py - Command line options handling.
""" 

# -- Python imports
import os, sys
import getopt

# -- External Modules
from mpfg.mp_args import *        # mpfg.Args: command line parsing and validation


# -------------------------------------------------------------------------------------------------
class PseArgs(Args):
   """! Management of command-line options. Based on parent class: mpf.Args. """
   
   def __init__(self):

      Args.__init__(self)
      
   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def default_config_dir(self):
      """! @return default configuration directory (overridden) """
      return "./config"

   @property
   def default_config_filename(self):
      """! @return default configuration filename (Overridden). """
      return "pse.cfg"  # pse

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Parse command-line input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def parse_input_arguments(self):
      """! Parse input arguments and options (Overridden) """

      success = Args.parse_input_arguments(self)   # here we do nothing more

      return success

   # ~~~~~~~~~~~~~~~~~~~~~~~~~
   # Validate input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~
   def validate_arguments(self):
      """! 
         Validate command-line input arguments (Overridden)
         @retval True if all arguments and options have been sucessfully validated
         @retval False if validation failed and in such case, also returns the corresponding error 
      """
      validated, err_msg = Args.validate_arguments(self)    # here we do nothing more

      return validated, err_msg

   # ~~~~~~~~~~~~~~~~~~~
   # Print Usage syntax
   # ~~~~~~~~~~~~~~~~~~~
   def print_usage(self, prog_name):
      """!
         Print Usage syntax (overridden)
         @param prog_name name of the main program
      """

      # Print usage information
      print('\nUsage: {0} [options]'.format(prog_name))
      print('\nHelp:')
      print('-h,  --help\tprint this help')

      # Optinal execution arguments
      print('\nOptions:')   
      print('-c,  --config-file\tconfiguration file name (default: pse.cfg)')    # pse
      print('-d,  --config-dir\tconfiguration directory (default: ./config)')
      print('\n')


# -- EOF pse_args.py
