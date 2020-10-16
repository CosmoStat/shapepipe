# -*- coding: utf-8 -*-

"""FILE SYSTEM TOOLS

This module defines methods for managing actions on the file system.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os


class FileSystemError(Exception):
   """ File System Error

   Generic error that is raised by the file system.

   """

   pass


def check_dir(dir_name):
    """Check Directory

    Check if directory exists.

    Parameters
    ----------
    dir_name : str
        Directory name

    """

    return os.path.isdir(dir_name)


def mkdir(dir_name, check_created=True, exist_ok=True):
    """ Make Directory

    This method creates a directory in the specified path.

    Parameters
    ----------
    dir_name : str
        Directory name with full path
    check_created : bool
        Check if directory is properly created or already exists, raise error
        if not found (default is True)
    exist_ok : bool
        If False raise an error if the directory alredy exists (default is
        True)

    Raises
    ------
    FileSystemError
        If directory already exists.
    FileSystemError
        If directory not properly created.

    """

    os.makedirs(dir_name, exist_ok=exist_ok)

    if check_created and not check_dir(dir_name):
        raise FileSystemError('Directory \"{}\" not found after mkdir command.'
                              ''.format(dir_name))
