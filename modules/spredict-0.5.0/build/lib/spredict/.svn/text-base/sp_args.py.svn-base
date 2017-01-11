"""! 
   sp_args.py - Command line options handling.
""" 

# -- Python imports
import os, sys
import getopt

# -- External Modules
from mpf.mp_args import *        # mpf.Args: command line parsing and validation


# -------------------------------------------------------------------------------------------------
class SpredictArgs(Args):
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
   def default_module_dir(self):
      """! @return default module directory"""
      return "./modules"

   @property
   def default_config_filename(self):
      """! @return default configuration filename (Overridden). """
      return "spredict.cfg"  # spredict

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Parse command-line input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def parse_input_arguments(self):
      """! Parse input arguments and options """

      # --- Parse command line options here (see getopt() documentation for syntax
      try:
         # --- Collect command-line arguments
         opts, self.args = getopt.getopt(sys.argv[1:], "c:d:m:qh", \
                                    ["config-file=", "config-dir=", "module-dir", "quiet", "help"])

         # --- Options: reference all key-value pairs in dictionary for later use
         for o, a in opts:
            if o in ("-h", "--help"):  # help
               
               self.add_option( ('-h','--help'), True)

            elif o in ("-c", "--config-file"):     # input configuration file
               self.add_option( ('-c','--config-file'), a)

            elif o in ("-d", "--config-dir"):      # input configuration directory
               self.add_option( ('-d','--config-dir'), a)

            elif o in ("-m", "--module-dir"):      # input module directory
               self.add_option( ('-m','--module-dir'), a)

            elif o in ("-q", "--quiet"):           # quiet mode
               self.add_option( ('-q','--quiet'), True)

      except getopt.GetoptError, err:

         # Print error, help information and exit:
         self.print_usage(self.helper.get_main_basename())

         raise Args.ParsingError(str(err))

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
      if validated:
         # Check module directory
         if '-m' in self.options:   
            if not os.path.isdir(self.options['-m']):
               err_msg = '{0} is not a valid directory'.format(self.options['-m'])
               validated = False
         else:
            # Set current directory as default config direcory
            self.add_option( ('-m','--module-dir'), self.default_module_dir)

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

      # Optional execution arguments
      print('\nOptions:')   
      print('-c,  --config-file\tconfiguration file name (default: spredict.cfg)')    # spredict
      print('-d,  --config-dir\tconfiguration directory (default: ./config)')
      print('-m,  --module-dir\tmodule directory (default: ./module)')
      print('\n')


# -- EOF sp_args.py
