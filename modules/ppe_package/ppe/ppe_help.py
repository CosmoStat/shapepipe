"""!
   @package ppe.ppe_helper helper class
   @author Marc Gentile
   @file ppe_help.py
   Helper class
"""

# -- Python imports
import os
import numpy as np
import pyfits
import glob

# -- External import
from scatalog import *
from mpfg.mp_helper import *
from ppe_version import __version__


# ----------------------------------------------------------------------------
class PpeHelper(Helper):

    """!
    Convenient utility functions

    """

    # ------------------------------------------------------------------------
    def __init(self):

        Helper.__init__(self)

    # ------------------------------------------------------------------------
    def get_version(self):
        """!
        Get the version number of this ppe code as text
        @return version number of this ppe code as text

        """

        return __version__

    # ------------------------------------------------------------------------
    def show_config_summary(self, master):

        if master.logging_enabled():

            # --- PPE version
            master.logger.log_info_p('\n*** PPE v{0} ***\n'.format(
                                     self.get_version()))

            # --- Python modules
            master.logger.log_info_p('Standard Python modules:')
            master.logger.log_info_p('- numpy {0}\t\t{1}'.format(
                                     np.__version__, np.__file__))
            master.logger.log_info_p('- pyfits {0}\t\t{1}'.format(
                                     pyfits.__version__, pyfits.__file__))
            try:
                mpl = __import__('matplotlib', globals(), locals(), [], -1)
                master.logger.log_info_p('- matplotlib {0}\t{1}'.format(
                                         mpl.__version__, mpl.__file__))
                asc = __import__('asciidata', globals(), locals(), [], -1)
                master.logger.log_info_p('- asciidata {0}\t{1}'.format(
                                         asc.__version__, asc.__file__))
            except Exception as detail:
                master.logger.log_error_p('- some modules could not be '
                                          'imported: {0}\n'.format(detail))

            master.logger.log_info_p('\nMPF Python modules:')

            try:
                mpfg = __import__('mpfg', globals(), locals(), [], -1)
                master.logger.log_info_p('- mpfg {0}\t\t{1}'.format(
                                         mpfg.__version__, mpfg.__file__))
                mpfx = __import__('mpfx', globals(), locals(), [], -1)
                master.logger.log_info_p('- mpfx {0}\t\t{1}'.format(
                                         mpfx.__version__, mpfx.__file__))
                slog = __import__('slogger', globals(), locals(), [], -1)
                master.logger.log_info_p('- slogger {0}\t\t{1}'.format(
                                         slog.__version__, slog.__file__))
                sconf = __import__('sconfig', globals(), locals(), [], -1)
                master.logger.log_info_p('- sconfig {0}\t\t{1}'.format(
                                         sconf.__version__, sconf.__file__))
                scat = __import__('scatalog', globals(), locals(), [], -1)
                master.logger.log_info_p('- scatalog {0}\t{1}'.format(
                                         scat.__version__, scat.__file__))
            except Exception as detail:
                master.logger.log_error_p('- some modules could not be '
                                          'imported: {0}\n'.format(detail))

            master.logger.log_info_p('\n')

    # ------------------------------------------------------------------------
    def _open_catalog(self, catalog_filepath, job, worker, hdu_no=1):

        catalog = None

        catalog = FITSCatalog(catalog_filepath, hdu_no=hdu_no,
                              open_mode=FITSCatalog.OpenMode.ReadWrite)

        catalog.open()

        return catalog

    # ------------------------------------------------------------------------
    def locate_files(self, master, pattern_list, directory, sort=True,
                     recurse=True, err_check=True):
        """!
        Locate files matching a search pattern list

        @param master master object instance
        @param pattern_list Unix-style file search pattern list (e.g.
        [*.fits, ".txt] )
        @param directory base directory from where to search for matching files
        @param sort [optional]  tell whether to sort the output file paths
        (default True)
        @param recurse [optional] tell whether to walk down directories
        (default True)
        @param err_check [optional] tell whether to the validity of check
        directories

        @return list of absolute paths of the files matching the search
        criteria
        @note the search is through the entire directory tree, not only the
        top nodes.

        """

        if err_check and not os.path.isdir(directory):
            self.helper.print_warning('{0} could not be found or is not a '
                                      'directory'.format(directory))
            return []

        filepaths = []

        if recurse:

            for pattern in pattern_list:
                # --- Recursively search the whole directory tree
                for filepath in self._walk_directory(master, pattern,
                                                     directory):
                    if not os.path.isdir(filepath):
                        filepaths.append(filepath)
        else:
            # --- Search only the top nodes of the directory
            for pattern in pattern_list:
                filepaths.extend([os.path.join(directory, f) for f in
                                 os.listdir(directory) if not
                                 os.path.isdir(os.path.join(directory, f)) and
                                 self.match_file(master, directory, f,
                                 pattern)])

        if sort:
            filepaths = sorted(filepaths)

        return list(set(filepaths))

    # ------------------------------------------------------------------------
    def match_file(self, master, directory, filename, pattern):
        """!
        File matching predicate method. Must return True in case of matching,
        False otherwise.
        May be overriden by subclasses to set additional criteria.
        @param master master object instance
        @param directory directory of filename
        @param filename file name
        @param pattern Unix-like file pattern
        @return True of a match is found, False otherwise

        """

        return glob.fnmatch.fnmatch(filename, pattern)

    # ~~~~~~~~~~~~~~~
    # Private methods
    # ~~~~~~~~~~~~~~~

    # -------------------------------------------------------------------------
    def _walk_directory(self, master, pattern, directory):
        """!
        Recursively locate files matching pattern in a given directory
        @param master master object instance
        @param pattern Unix-style file search pattern (e.g. *.fits)
        @param directory: base directory from where to search for matching
        files

        @return list of files matching the search criteria

        """

        for path, dirs, files in os.walk(directory):
            for filename in [os.path.abspath(os.path.join(path, filename)) for
                             filename in files if self.match_file(master,
                             path, filename, pattern)]:
                yield filename

# -- EOF ppe_helper.py
