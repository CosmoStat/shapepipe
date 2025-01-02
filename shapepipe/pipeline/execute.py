"""EXECUTE.

This module defines methods for running the pipeline modules.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
import subprocess as sp


def execute(command_line):
    """Execute.

    This method executes a given command line.

    Parameters
    ----------
    command_line : str
        The command line to be executed

    Returns
    -------
    tuple
        Stdout and stderr (both type str)

    Raises
    ------
    TypeError
        For invalid input type

    """
    if not isinstance(command_line, str):
        raise TypeError("Command line must be a string.")

    command = command_line.split()
    check_executable(command[0])

    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = process.communicate()

    return stdout.decode("utf-8"), stderr.decode("utf-8")


def check_executable(exe_name):
    """Check if Input is Executable.

    This methid checks if the input executable exists.

    Parameters
    ----------
    exe_name : str
        Executable name

    Raises
    ------
    TypeError
        For invalid input type
    OSError
        For non-existent executable

    """
    if not isinstance(exe_name, str):
        raise TypeError("Executable name must be a string.")

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(exe_name)

    if not fpath:
        res = any(
            [
                is_exe(os.path.join(path, exe_name))
                for path in os.environ["PATH"].split(os.pathsep)
            ]
        )

    else:
        res = is_exe(exe_name)

    if not res:
        raise OSError(
            f"{exe_name} does not appear to be a valid executable on this "
            + "system."
        )
