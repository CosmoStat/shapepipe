"""! 
   @package mpfg.mp_args Command line options management for mp_MPI.py or mp_SMP.py
   @author Marc Gentile
   @file mp_args.py
   Command line options management for mp_MPI.py or mp_SMP.py
""" 

# -- Python imports
import os, sys
import getopt

# -- Module-specific imports
from mp_helper import *       # miscellaneous utility functions


# -------------------------------------------------------------------------------------------------
class Args(object):
   """! Management of command-line options """
   
   def __init__(self):
      """! Default class constructor """

      self._args = []                  # optional command-line arguments
      self._options = {}               # command-line options      
      self._helper = Helper()          # helper utility functions
   
      # --- Parse and validate command-line args and options 
      self.parse_command_line()


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def args(self):
      """! @return command-line arguments as a list """
      return self._args

   @property
   def options(self):
      """! @return command-line options as a dictionary """
      return self._options

   @property
   def default_config_dir(self):
      """! @return default configuration directory (e.g. ../config) """
      return "./config"

   @property
   def default_config_filename(self):
      """! @return default configuration filename (e.g. multiproc.cfg) """
      return "mpfg.cfg"

   # --- Setters

   @args.setter 
   def args(self, args):
      """!
         set command-line arguments 
         @param args list of arguments   
      """
      self._args = args 

   @options.setter 
   def options(self, options):
      """!
         set command-line options
         @param options list of options   
      """
      self._options = options 

   @default_config_dir.setter 
   def default_config_dir(self, default_config_dir):
      """! 
         Set default configuration directory (e.g. ../config) 
         @param default_config_dir default configuration directory
      """
      self._default_config_dir = default_config_dir 

   @default_config_filename.setter 
   def default_config_filename(self, default_config_filename):
      """! 
         Set default configuration filename (e.g. multiproc.cfg) 
         @param default_config_filename default configuration filename
      """
      self._default_config_filename = default_config_filename 


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


   # ~~~~~~~~~~~~~~~~~~~~~~~
   # Get configuration info 
   # ~~~~~~~~~~~~~~~~~~~~~~~
   def parse_command_line(self):
      """! Parse and validate command-line input arguments """

      # Parse input arguments
      self.parse_input_arguments()

      # Check input arguments
      success, err_msg = self.validate_arguments()
      if not success:

         # Some error occurred
         # MK new: print usage string
         self.print_usage(self.helper.get_main_basename())
         raise Args.ValidationError(err_msg)

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Parse command-line input arguments
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def parse_input_arguments(self):
      """! Parse input arguments and options """

      # --- Parse command line options here (see getopt() documentation for syntax
      try:
         # --- Collect command-line arguments
         opts, self.args = getopt.getopt(sys.argv[1:], "c:d:qh", ["config-file=", "config-dir=", "quiet", "help"])

         # --- Options: reference all key-value pairs in dictionary for later use
         for o, a in opts:
            if o in ("-h", "--help"):  # help
               
               self.add_option( ('-h','--help'), True)

            elif o in ("-c", "--config-file"):     # input configuration file
               self.add_option( ('-c','--config-file'), a)

            elif o in ("-d", "--config-dir"):      # input configuration directory
               self.add_option( ('-d','--config-dir'), a)

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
         Validate command-line input arguments (Override)
         @retval True if all arguments and options have been sucessfully validated
         @retval False if validation failed and in such case, also returns the corresponding error 
      """

      validated = True
      err_msg = ''
      
      # --- Check presence of mandatory parameters
   #   if len(args) == 0:
   #      err_msg = 'one or more command line arguments must be specified'
   #      validated = False

      if validated:
         # Check config filename
         if '-d' in self.options:   
            if not os.path.isdir(self.options['-d']):
               err_msg = '{0} is not a valid directory'.format(self.options['-d'])
               validated = False
         else:
            # Set current directory as default config direcory
            self.add_option( ('-d','--config-dir'), self.default_config_dir)

      if validated:
         if not '-c' in self.options:
            # Set current directory as default config direcory
            self.add_option( ('-c','--config-file'), self.default_config_filename)

         # Check config directory
         config_filepath = os.path.join(self.options['-d'], self.options['-c'])
         if not self.helper.file_exists(config_filepath):
            err_msg = '{0} does not exist or is not a valid file'.format(config_filepath)
            validated = False

      if validated:
         # Check help
         if not '-h' in self.options:   
            self.add_option( ('-h','--help'), False)

         # Check quiet mode
         if not '-q' in self.options:   
            self.add_option( ('-q','--quiet'), False)

      return validated, err_msg


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Register a command-line option
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def add_option(self, (short_opt, long_opt), value):
      """!
         Record input command line arguments 
         @param short_opt String: short option (like "-q")
         @param long_opt String: long option (like "--quiet")
         @param value String: option value
      """

      if not short_opt in self.options and not long_opt in self.options:
         self.options[short_opt] = self.options[long_opt] = value
      else:
         self.helper.print_warning('Option {0} or {1} is redundant'.format(short_opt, long_opt)) 

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
      #print('(C) 2013 astro craft')
      print('\nUsage: {0} [options]'.format(prog_name))
      print('\nHelp:')
      print('-h,  --help\tprint this help')

      # Optinal execution arguments
      print('\nOptions:')   
      print('-c,  --config-file\tconfiguration file')
      print('-d,  --config-dir\tconfiguration directory')
      print('\n')

   # ~~~~~~~~~~~~~~~~~~~~~~~
   # Get Main program name 
   # ~~~~~~~~~~~~~~~~~~~~~~~
   def get_prog_name(self):
      """! Return the name of the main program """
      return self.helper.get_main_basename()


   # ------------------------------------------------------------------------------
   class ParsingError(Exception):

      """! Exception thrown in case of error detected while parsing the command-line """
      def __init__(self, msg):
         """!
            Class constructor
            @param msg message to print upon exception
         """
         self._msg = msg   

      def __str__(self):
         """! 
             Formatting for display 
             @return string representation of the exception  
         """
         return "Args *** ERROR ***: {0}".format(self._msg)

   # ------------------------------------------------------------------------------
   class ValidationError(Exception):

      """! Exception thrown in case the validating of the arguments and options failed """
      def __init__(self, msg):
         """!
            Class constructor
            @param msg message to print upon exception
         """
         self._msg = msg   

      def __str__(self):
         """!
             Formatting for display 
             @return string representation of the exception  
         """
         return "Args *** ERROR ***: {0}".format(self._msg)

# -- EOF mp_args.py
