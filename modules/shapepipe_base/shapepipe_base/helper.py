# -*- coding: utf-8 -*-

"""HELPER MODULE

This module contain a class with 'helper' methods.

:Authors: Samuel Farrens and Marc Gentile

:Date: 09/04/2018

"""

# -- Python imports
import os
import glob
import subprocess

# -- External import
from mpfg.mp_helper import Helper


def is_executable(exe_name):
    """Check if Input is Executable

    This methid checks if the input executable exists.

    Parameters
    ----------
    exe_name : str
        Executable name

    Returns
    -------
    Bool result of test

    Raises
    ------
    TypeError
        For invalid input type

    """

    if not isinstance(exe_name, str):

        raise TypeError('Executable name must be a string.')

    def is_exe(fpath):

        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(exe_name)

    if not fpath:

        res = any([is_exe(os.path.join(path, exe_name)) for path in
                   os.environ["PATH"].split(os.pathsep)])

    else:

        res = is_exe(exe_name)

    if not res:
        raise IOError('{} does not appear to be a valid executable on this '
                      'system.'.format(exe_name))


class PackageHelper(Helper):

    """Package Helper

    This class contains convenient utility methods.

    """

    def __init__(self, version, name, pydepend, sysdepend):

        Helper.__init__(self)
        self.version = version
        self.name = name
        self.pydepend = pydepend
        self.sysdepend = sysdepend

    def _check_python_dependencies(self, master):

        master.logger.log_info_p('- Python dependencies.')

        try:

            for module_name in self.pydepend:

                module = __import__(module_name, globals(), locals(), [],
                                    -1)

                master.logger.log_info_p('-- {0} {1}\t{2}'.format(
                                         module_name, module.__version__,
                                         module.__file__))

        except Exception as detail:
            master.logger.log_error_p('- some modules could not be '
                                      'imported: {0}\n'.format(detail))

    def _check_system_dependencies(self, master):

        master.logger.log_info_p('- System dependencies.')

        try:

            for exe_name in self.sysdepend:

                is_executable(exe_name)

                exe_path, err = (subprocess.Popen('which {0}'.format(exe_name),
                                 shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate())

                master.logger.log_info_p('-- {0} {1}'.format(exe_name,
                                         exe_path.rstrip()))

        except Exception as detail:
            master.logger.log_error_p('- some executables could not be '
                                      'found: {0}\n'.format(detail))

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

            master.logger.log_info_p('Checking package dependencies.')
            self._check_python_dependencies(master)
            self._check_system_dependencies(master)
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
