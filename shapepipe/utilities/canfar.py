# -*- coding: utf-8 -*-

"""CANFAR TOOLS

This module defines methods for managing CANFAR specific actions.

:Author: Samuel Farrens <samuel.farrens@cea.fr>
         Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import os
import sys
from io import StringIO
from contextlib import redirect_stdout

try:
    import vos.commands as vosc
except ImportError:  # pragma: no cover
    import_fail = True
else:
    import_fail = False


class vosError(Exception):
    """VOS Error

    Generic error that is raised by the vosHandler.

    """
    pass


class vosHandler:
    """VOS Handler

    This class manages the use of VOS commands.

    Parameters
    ----------
    command : str
        VOS command name

    """

    def __init__(self, command):

        self._check_vos_install()
        self._avail_commands = tuple(vosc.__all__)
        self.command = command

    @staticmethod
    def _check_vos_install():
        """Check VOS Install

        Check if VOS is correctly installed.

        Raises
        ------
        ImportError
            if vos package cannot be imported
        """
        if import_fail:
            raise ImportError(
                'vos package not found, re-install ShapePipe '
                + 'with \'./install_shapepipe --vos\''
            )

    @property
    def command(self):
        """Command

        This method sets the VOS command property.

        Raises
        ------
        ValueError
            if value is not valid vos command
        """
        return self._command

    @command.setter
    def command(self, value):

        if value not in self._avail_commands:
            raise ValueError(
                f'vos command must be one of {self._avail_commands}'
            )

        self._command = getattr(vosc, value)

    def __call__(self, *args, **kwargs):
        """Call Method

        This method allows class instances to be called as functions.

        Raises
        ------
        vosError
            if error in vos command occurs
        """
        try:
            self._command()

        except:
            raise vosError(
                f'Error in VOs command: {self._command.__name__}'
            )


def download(source, target, verbose=False):
    """Download file from vos.

    Parameters
    ----------
    source : str
        source path on vos
    target : str
        target path
    verbose : bool, optional, default=False
        verbose output if True

    Returns
    -------
    status : bool
        status, True/False or success/failure
    """

    cmd = 'vcp'

    if not os.path.exists(target):
        sys.argv = [cmd, source, target]
        if verbose:
            print(f'Downloading file {source} to {target}...')
        vcp = vosHandler(cmd)

        vcp()
        if verbose:
            print('Download finished.')
    else:
        if verbose:
            print(f'Target file {target} exists, skipping download.')


def dir_list(path, verbose=False):
    """list

    List content of path on vos

    Parameters
    ----------
    path : str
        path on vos, starts with 'vos:cfis/...'
    verbose : bool, optional, default=False
        verbose output if True

    Raises
    ------
    HTTPError, KeyError
        if error occurs during vos command

    Returns
    -------
    vls_out : array of str
        file or directory at path
    """

    cmd = 'vls'
    sys.argv = [cmd, path]
    vls = vosHandler(cmd)

    if verbose:
        print('Getting vos directory content from vls...')

    f = StringIO()

    try:
        with redirect_stdout(f):
            vls()
    except:
        print('Error during vls command')
        raise

    vls_out = f.getvalue()

    return vls_out.split('\n')
