# -*- coding: utf-8 -*-

"""HELPER MODULE

This module contain a class with 'helper' methods.

:Authors: Samuel Farrens and Marc Gentile

:Date: 31/10/2017 (Happy Halloween!)

"""

# -- Python imports
import os
import numpy as np
import glob

# -- External import
from mpfg.mp_helper import Helper
from info import __version__, __whoami__


class PackageHelper(Helper):

    """Package Helper

    This class contains convenient utility methods.

    """

    def __init__(self):

        Helper.__init__(self)
        self.version = __version__
        self.name = __whoami__

    def show_config_summary(self, master):

        """Show Summary of Configuration

        Parameters
        ----------
        master : class
            Master class instance

        """

        if master.logging_enabled():

            # --- Package name and version
            master.logger.log_info_p('\n*** {0} v{1} ***\n'.format(self.name,
                                     self.version))

            try:

                # --- Python modules
                master.logger.log_info_p('Standard Python modules:')

                master.logger.log_info_p('- numpy {0}\t\t{1}'.format(
                                         np.__version__, np.__file__))

                mpl = __import__('matplotlib', globals(), locals(), [], -1)
                master.logger.log_info_p('- matplotlib {0}\t{1}'.format(
                                         mpl.__version__, mpl.__file__))

                master.logger.log_info_p('\nMPF Python modules:')

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

    def locate_files(self, pattern_list, directory, sort=True,
                     recurse=True, err_check=True):

        """Locate Files

        Locate files matching a search pattern list.

        Parameters
        ----------
        pattern_list : list
            Unix-style file search pattern list (e.g.[*.fits, ".txt] )
        directory : str
            Base directory from where to search for matching files
        sort : bool, optional
            Tell whether to sort the output file paths (default True)
        recurse : bool, optional
            Tell whether to walk down directories (default True)
        err_check : bool, optional
            Tell whether to the validity of check directories

        Returns
        -------
        list of absolute paths of the files matching the search criteria

        Notes
        -----
        The search is through the entire directory tree, not only the top
        nodes.

        """

        if err_check and not os.path.isdir(directory):
            self.helper.print_warning('{0} could not be found or is not a '
                                      'directory'.format(directory))
            return []

        filepaths = []

        if recurse:
            for pattern in pattern_list:
                # --- Recursively search the whole directory tree
                for filepath in self._walk_directory(pattern, directory):
                    if not os.path.isdir(filepath):
                        filepaths.append(filepath)
        else:
            # --- Search only the top nodes of the directory
            for pattern in pattern_list:
                filepaths.extend([os.path.join(directory, f) for f in
                                 os.listdir(directory) if not
                                 os.path.isdir(os.path.join(directory, f)) and
                                 self._match_file(f, pattern)])

        if sort:
            filepaths = sorted(filepaths)

        return list(set(filepaths))

    @staticmethod
    def _match_file(filename, pattern):

        """Match File

        This method checks if the input filename matches the required Unix
        style pattern.

        Parameters
        ----------
        filename : str
            Input file name
        pattern : str
            Unix style file pattern

        Returns
        -------
        bool True if a match is found, False otherwise

        """

        return glob.fnmatch.fnmatch(filename, pattern)

    def _walk_directory(self, pattern, directory):

        """Walk Directory

        Recursively locate files matching pattern in a given directory.

        Parameters
        ----------
        pattern : str
            Unix style file pattern
        directory : str
            Base directory from where to search for matching files

        Returns
        -------
        list of files matching the search criteria

        """

        for path, dirs, files in os.walk(directory):
            for filename in [os.path.abspath(os.path.join(path, filename)) for
                             filename in files if self._match_file(filename,
                             pattern)]:
                yield filename
