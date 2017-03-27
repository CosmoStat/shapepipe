"""! 
   mkp_pse.py - PSF Reconstruction
"""

# -- Python imports
import os, sys
import imp
import time
import numpy
from operator import itemgetter, attrgetter

# -- External imports
from mpfg3.mpfg3_job import *    # base job processing for GREAT3
from mpfg3.mpfg3_data import *   # data access to GREAT3

# --- Module-specific imports
from mkp_help import *           # helper utility functions


## -------------------------------------------------------------------------------------------------
#class PSFBuilder(object):
#   
#   """! 
#      PSFBuilder: Reconstruct a PSF field.   
#   """

#   def __init__(self, method_name=None, method_dir="./methods"):
#      """! Construct a PSFBuilder class """

#      self._build_method_name = method_name
#      self._build_method_dir  = method_dir
#      self._helper = MkpHelper()    # helper utility functions
#      
#   # ~~~~~~~~~~
#   # Properties 
#   # ~~~~~~~~~~

#   @property
#   def build_method_name(self):
#      """! @return the BuildMethod name. """
#      return self._build_method_name

#   @property
#   def build_method_dir(self):
#      """! @return the BuildMethod directory. """
#      return self._build_method_dir

#   @property
#   def helper(self):
#      """! @return the MkpHelper instance. """
#      return self._helper


#   # ~~~~~~~~~~~~~~~
#   # Public methods 
#   # ~~~~~~~~~~~~~~~

#   def build_PSF(self, job, worker):

#      try:

#         # --- Get relevant configuration info
#         module = self._load_method_module(self.build_method_name, self.build_method_dir)   # method-specific code



##         # --- Get relevant configuration info
##         build_method_name = worker.config.get_as_string("BUILD_METHOD", "METHODS")
##         base_method_dir = worker.config.get_as_string("BASE_METHOD_DIR", "METHODS")

##         # --- Setup the PSF build method
##         build_method = PSFBuildMethod(build_method_name, base_method_dir)

#         # --- Apply the method to reconstruct the PSF    
#         pass

#      except(PSFBuilder.BuildingError):

#         if worker.logging_enabled():
#            worker.logger.log_error_p(
#                        "{0} - Some error occurred while processing job: {1} ({2})".format(
#                                                                                worker.name, 
#                                                                                job,
#                                                                                sys.exc_info()[1]))


#   # ~~~~~~~~~~~~~~~
#   # Private methods 
#   # ~~~~~~~~~~~~~~~

##   def _load_method_config(self, base_method_dir):
##      """! Load the method configuration file """
##      
##      # Config name is derived from the method name 
##      return os.path.join(base_method_dir, self.name) + '.cfg'  
#  

#   # -----------------------------------------------------------------------------------------------
#   def _load_method_module(self, method_name, method_dir):
#      """! Load the Python module to access the method's code """   

#      module_code_name = method_name
#      #module_path = os.path.join(base_method_dir, module_code_name)

##      return self._import_module("numpy", "/usr/local/lib/python2.6/dist-packages")
#      return self._import_module(module_code_name, method_dir)

#   # -----------------------------------------------------------------------------------------------
#   def _import_module(self, module_name, module_dir):
#      """! Find and load a module <module_name> from a directory <module_dir> """

#      try:


#         file_obj, filename, data = imp.find_module(module_name, [module_dir])

#         print "LOADED:", file_obj, filename, data
#         imp.load_module(module_name, file_obj, filename, data)

#         module = imp.load_module(module_name, file_obj, filename, data)

#         print "MODULE:", module

#         return module

#      except:
#         raise(PSFBuilder.BuildingError("Could not load module {0}.".format(module_name)))

#   # ~~~~~~~~~~
#   # Exceptions 
#   # ~~~~~~~~~~ 

##   # -----------------------------------------------------------------------------------------------
#   class BuildingError(Exception):
#      """! 
#         Exception thrown when a PSF reconstruction error has occurred
#      """
#      def __init__(self, msg):
#         """!
#            @param msg error message 
#         """
#         self._msg = msg
#         self._exc_info = sys.exc_info()[1]

#         print self._msg, self._exc_info

#      def __str__(self):
#         """! String representation of the BuildingError class """
#         if self._exc_info is not None:
#            return "MKPSF *** ERROR ***: {0} ({1})".format(self._msg, self._exc_info)
#         else:
#            return "MKPSF *** ERROR ***: {0}".format(self._msg)


## -------------------------------------------------------------------------------------------------
#class PSFBuildMethod(object):
#   
#   """! 
#      PSFBuildMethod: method for reconstructing a PSF field
#   """

#   def __init__(self, method_name, method_dir):
#      """! Construct a PSFBuildMethod class """
#   
#      self._name = method_name   # reconstruction method name -> xxx.cfg, xxx.py 

#      # --- Get relevant configuration info
#      self._config = self._load_method_config(method_dir)   # method-specific configuration
#      self._module = self._load_method_module(method_dir)   # method-specific code

##      print self.module

#   # --- Getters

#   @property
#   def name(self):
#      """! @return the name of the build method """   
#      return self._name

#   @property
#   def config(self):
#      """! @return the configuration object of the build method """   
#      return self._config

#   @property
#   def module(self):
#      """! @return the instantiated Python module of the build method """   
#      return self._module


#   # ~~~~~~~~~~~~~~~
#   # Public methods 
#   # ~~~~~~~~~~~~~~~

#   # ~~~~~~~~~~~~~~~
#   # Private methods 
#   # ~~~~~~~~~~~~~~~


#   def _load_method_config(self, base_method_dir):
#      """! Load the method configuration file """
#      
#      # Config name is derived from the method name 
#      return os.path.join(base_method_dir, self.name) + '.cfg'  
#  

#   def _load_method_module(self, base_method_dir):
#      """! Load the Python module to access the method's code """   

#      module_code_name = self.name
#      #module_path = os.path.join(base_method_dir, module_code_name)

#      return self._import_module("numpy", "/usr/local/lib/python2.6/dist-packages")
##      return self._import_module(module_code_name, base_method_dir)

#   def _import_module(self, module_name, module_dir):
#      """! Find and load a module <module_name> from a directory <module_dir> """

#      try:

#         print "MODULE:", module_name, module_dir

#         file_obj, filename, data = imp.find_module(module_name, [module_dir])

#         print "LOADED:", file_obj, filename, data
#         imp.load_module(module_name, file_obj, filename, data)
#         return imp.load_module(module_name, file_obj, filename, data)

#      except:
#         raise(PSFBuildMethod.BuildMethodError("Could not load module {0}.".format(module_name)))



#   # -----------------------------------------------------------------------------------------------
#   class BuildMethodError(Exception):
#      """! 
#         Exception thrown related to the PSF build method
#      """
#      def __init__(self, msg):
#         """!
#            @param msg error message 
#         """
#         self._msg = msg
#         self._exc_info = sys.exc_info()[1]

#      def __str__(self):
#         """! String representation of the BuildMethodError class """
#         if self._exc_info is not None:
#            return "MKPSF *** ERROR ***: {0} ({1})".format(self._msg, self._exc_info)
#         else:
#            return "MKPSF *** ERROR ***: {0}".format(self._msg)





# -- EOF mkp_psf.py
