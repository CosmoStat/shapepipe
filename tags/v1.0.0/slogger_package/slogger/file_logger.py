"""! 
   @package slogger.file_logger Simple and lightweight file logging facility
   @author Marc Gentile
   @file file_logger.py
   File logging management
"""

"""! 
   @package mpfx.mpfx_data - Dataset management
   @author Marc Gentile
   @file mpfx_data.py
   Dataset management
""" 

import os, sys 
import socket
import string
import time

# ~~~~~~~~~~~~~~~~~~~~~~~
# Base file Logger class
# ~~~~~~~~~~~~~~~~~~~~~~~
class BaseFileLogger(object):

   """!
      Simple file logging management - Base class 
   """

   def __init__(self, logger='', directory='.', filename='default.log', open_at_creation=True):
      """! 
         Construct a BaseFileLoggger object
         @param logger a prefix for each line inserted in the log
         @param directory output directory where the log is to be created
         @param filename name of the log file
         @param open_at_creation if @c True, the log file is also opened once created
      """

      self._logger_name   = logger
      self._log_directory = directory
      self._log_filename  = filename
      self._log_fullpath  = os.path.join(directory, filename)
      self._log_filedescr = None
      self._host_name = socket.gethostname()

      self._open_at_creation = open_at_creation
      if self._open_at_creation:
         self.open()

   @property
   def logname(self):
      """!
         Return the name of the log prefix 
         @return the name of the log prefix 
      """
      return self._logger_name

   @property
   def directory(self):
      """! Return the directory where the log resid
         @return the directory where the log reside 
      """
      return self._log_directory

   @property
   def filename(self):
      """!
         Return the name of the log file 
         @return the name of the log file 
      """
      return self._log_filename

   @property
   def fullpath(self):
      """! 
         Return the full path of the log file 
         @return the full path of the log file 
      """
      return self._log_fullpath

   @property
   def hostname(self):
      """! 
         Return the name of the machine where the log reside
         @return the name of the machine where the log reside 
      """
      return self._host_name

   @property
   def filedescr(self):
      """!
         Return the file descriptor of the open log file  
         @return the file descriptor of the open log file 
      """
      return self._log_filedescr

   @property
   def open_at_creation(self):
      """!
         Tell if the log file has been open at creation
         @retval True if the log is has been open at creation
         @retval False otherwise """
      return self._open_at_creation

   def open(self):
      """! Open the log file """
      try:

         # Open file
         if not os.path.exists(self._log_directory):
            os.mkdir(self._log_directory)
         self._log_filedescr = open(self._log_fullpath, "w")

         # Write header information 
         self._log_filedescr.write('=== Host: %s - Log file %s created: %s\n\n' %(self._host_name, self._log_fullpath, time.strftime('%d.%m.%y_%H:%M:%S', time.localtime())))
         os.fsync(self._log_filedescr.fileno())

      except:
         print("SFileLogger: file {0} could not be created ({1})".format(self._log_fullpath, sys.exc_info()[1]))
         if not self._log_filedescr is None:
            self._log_filedescr.close()

   def close(self):
      """! Close the log file """
      try:
         self._log_filedescr.write('\n=== Host: %s - Log file closed: %s\n\n' %(self._host_name, time.strftime('%d.%m.%y %H:%M:%S', time.localtime())))
         os.fsync(self._log_filedescr.fileno())
         self._log_filedescr.close()
      except:
         print("SFileLogger: file {0} could not be closed ({1})".format(self._log_fullpath, sys.exc_info()[1]))

   def is_open(self):
      """! 
          Tells if the log file is currently open
          @retval True if opened   
          @retval False if closed   
      """
      return not (self._log_filedescr is None)

   def flush(self):
      """! Force a flush of the log file buffer to disk """
      self._log_filedescr.flush()
      os.fsync(self._log_filedescr.fileno())


# ~~~~~~~~~~~~~~~~~~~~~
# A Simple file Logger
# ~~~~~~~~~~~~~~~~~~~~~
class FileLogger(BaseFileLogger):

   """!
      Simple file logging management: File Logger
   """

   def __init__(self, logger='', directory='.', filename='default.log', open_at_creation=True):
      """! 
         Construct a FileLoggger object
         @param logger a prefix for each line inserted in the log
         @param directory output directory where the log is to be created
         @param filename name of the log file
         @param open_at_creation if @c True, the log file is also opened once created
      """
      BaseFileLogger.__init__(self, logger, directory, filename, open_at_creation)
      
   def log_info(self, message):
      """!
         Log an information message
         @param message message to log
         @see log_info_p()
      """      
      self._log_filedescr.write('{0}{1} {2}\n'.format(self._logger_name, time.strftime("%H:%M:%S", time.localtime()), message))

   def log_info_p(self, message):
      """!
         Log an information message and display the message to the console
         @param message message to log and display
         @see log_info()
      """      
      self.log_info(message)
      print('{0}'.format(message))

   def log_warning(self, message):
      """!
         Log an warning message
         @param message message to log
         @see log_warning_p()
      """ 
      self._log_filedescr.write('{0}{1} *** Warning *** {2}\n'.format(self._logger_name, time.strftime("%H:%M:%S", time.localtime()), message))

   def log_warning_p(self, message):
      """!
         Log a warning message and display the message to the console
         @param message message to log and display
         @see log_warning()
      """ 
      self.log_warning(message)
      print('*** Warning *** {0}'.format(message))

   def log_error(self, message):
      """!
         Log an error message
         @param message message to log
         @see log_error_p()
      """ 
      self._log_filedescr.write('{0}{1} *** ERROR *** {2}\n'.format(self._logger_name, time.strftime("%H:%M:%S", time.localtime()), message))

   def log_error_p(self, message):
      """!
         Log an error message and display the message to the console
         @param message message to log and display
         @see log_error()
      """ 
      self.log_error(message)
      print('*** ERROR *** {0}'.format(message))

   def log_blank_line(self):
      """!
         Insert a blank line
      """ 
      self._log_filedescr.write('\n')

   def is_null(self):
      """!
         Tells is this is an instance of a NullFileLogger
         @retval False
         @note always return @c False for this class
         @see NullFileLogger
      """ 
      return False

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# A do-nothing file Logger
# ~~~~~~~~~~~~~~~~~~~~~~~~~
class NullFileLogger(BaseFileLogger):

   """!
      Simple file logging management - "do-nothing" file logger
   """

   def __init__(self, logger='', directory='.', filename='default.log', open_at_creation=False):
      """! 
         Construct a NullFileLoggger object
         @param logger a prefix for each line inserted in the log
         @param directory output directory where the log is to be created
         @param filename name of the log file
         @param open_at_creation if @c True, the log file is also opened once created
      """
      pass

   def open(self):
      """! Simulate the opening of a log file 
         @note does nothing for this class
         @retval None
      """
      return None

   def close(self):
      """! Simulate the closing of a log file 
         @note does nothing for this class
      """
      pass

   def is_open(self):
      """! 
         Tells if the log file is currently open
         @note always returns False for this class
         @retval False   
      """
      return False

   def is_null(self):
      """!
         Tells is this is an instance of a NullFileLogger
         @note always return @c True for this class
         @retval True
      """ 
      return True

   def log_info(self, message):
      """!
         Log an information message and display the message to the console
         @param message message to log and display
         @note does nothing for this class
         @see log_info_p()
      """
      pass

   def log_info_p(self, message):
      """!
         Log an information message and display the message to the console
         @param message message to log and display
         @note does nothing for this class
         @see log_info()
      """ 
      pass

   def log_warning(self, message):
      """!
         Log a warning message and display the message to the console
         @param message message to log and display
         @note does nothing for this class
         @see log_warning_p()
      """
      pass

   def log_warning_p(self, message):
      """!
         Log a warning message and display the message to the console
         @param message message to log and display
         @note does nothing for this class
         @see log_warning()
      """ 
      pass

   def log_error(self, message):
      """!
         Log an error message and display the message to the console
         @param message message to log and display
         @see log_error_p()
         @note does nothing for this class
      """      
      pass

   def log_error_p(self, message):
      """!
         Log an error message and display the message to the console
         @param message message to log and display
         @note does nothing for this class
         @see log_error()
      """ 
      pass

   def log_blank_line(self):
      """!
         Insert a blank line
         @note does nothing for this class
      """ 
      pass

   def flush(self):
      """! 
         Force a flush of the log file buffer to disk 
         @note does nothing for this class
      """
      pass


# ~~~~~~~~~~~~~
# Unit Testing
# ~~~~~~~~~~~~~
def main():
   """! Performs a few tests """

   log_time = time.strftime('%d.%m.%y_%H.%M.%S', time.localtime())
   log_file_name = ('img_{0:03d}_{1}.log'.format(1,  log_time))
   logger = FileLogger(logger='', directory='.', filename=log_file_name)

   print 'logname:', logger.logname
   print 'hostname:', logger.hostname
   print 'directory:', logger.directory
   print 'filename:', logger.filename
   print 'fullpath:', logger.fullpath
   print 'filedescr:', logger.filedescr
   print 'open_at_creation:', logger.open_at_creation

   logger.log_info_p('%s - Image %03d - Processing data ...' %("string", 1))
   logger.log_info_p('{0} - Image {1:03d} - Processing data ...'.format("message", 2))
   logger.log_error('%s - Image %03d - Error while processing data ...' %("string", 1))

   must_log = False

   # --- Just Log
   start = time.clock()
   for i in xrange (1000000):
      logger.log_info('{0} - Image {1:03d} - Processing data ...'.format("message", 2))
   elapsed = time.clock() - start
   print("Just log:", elapsed)

   # --- If ... Log
   start = time.clock()
   for i in xrange (1000000):
      if must_log:
         logger.log_info('{0} - Processing data ...'.format("message", 2))
   elapsed = time.clock() - start
   print("If log:", elapsed)

   # --- No Log
   start = time.clock()
   for i in xrange (1000000):
      pass
   elapsed = time.clock() - start
   print("Nothing:", elapsed)

   logger.flush()
   logger.close()

if __name__ == "__main__":
    main()

