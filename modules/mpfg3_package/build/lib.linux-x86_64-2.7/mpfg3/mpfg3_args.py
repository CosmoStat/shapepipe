"""! 
   @package mpfg3.mpfg3_args Command line option management
   @author Marc Gentile
   @file mpfg3_args.py
   Command line option management
""" 

# -- Python imports
import os, sys
import getopt

# -- External Modules
from mpfx.mpfx_args import *        # mpfx.Args: command line parsing and validation

# -------------------------------------------------------------------------------------------------
class Mpfg3Args(MpfxArgs):
   """! Management of command-line options """
   
   def __init__(self):
      """! Default class constructor """
      MpfxArgs.__init__(self)

      
   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def default_config_dir(self):
      """! @return default configuration directory """

      return "./config"

   @property
   def default_config_filename(self):
      """! @return default configuration filename """

      return "mpfg3.cfg"

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Parse command-line input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def parse_input_arguments(self):
      """! Parse input arguments and options (Override) """

      MpfxArgs.parse_input_arguments(self)


   # ~~~~~~~~~~~~~~~~~~~~~~~~~
   # Validate input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~
   def validate_arguments(self):
      """! 
         Validate command-line input arguments (Override)
         @retval True if all arguments and options have been sucessfully validated
         @retval False if validation failed and in such case, also returns the corresponding error 
      """
      validated, err_msg = MpfxArgs.validate_arguments(self)

      return validated, err_msg

   # ~~~~~~~~~~~~~~~~~~~
   # Print Usage syntax
   # ~~~~~~~~~~~~~~~~~~~
   def print_usage(self, prog_name):
      """!
         Print Usage syntax 
         @param prog_name name of the main program
      """

      # Print usage information
      #print('\n{0}'.format(prog_name))
      print('\nUsage: {0} [options]'.format(prog_name))
      print('\nHelp:')
      print('-h,  --help\tprint this help')

      # Optinal execution arguments
      print('\nOptions:')   
      print('-c,  --config-file\tconfiguration file name (default: mpfg3.cfg)')
      print('-d,  --config-dir\tconfiguration directory (default: ./config)')
      print('\n')


# -- EOF mpfg3_args.py
