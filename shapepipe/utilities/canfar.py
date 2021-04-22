# -*- coding: utf-8 -*-

"""CANFAR TOOLS

This module defines methods for managing CANFAR specific actions.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

try:
    import vos.commands as vosc
except ImportError:  # pragma: no cover
    import_fail = True
else:
    import_fail = False


class vosError(Exception):
   """ VOS Error

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

        """

        if import_fail:
            raise ImportError('vos package not found, re-install ShapePipe '
                              'with \'./install_shapepipe --vos\'')

    @property
    def command(self):
        """Command

        This method sets the VOS command property.

        """

        return self._command

    @command.setter
    def command(self, value):

        if value not in self._avail_commands:
            raise ValueError('vos command must be one of {}'
                             ''.format(self._avail_commands))

        self._command = getattr(vosc, value)

    def __call__(self, *args, **kwargs):
        """Call Method

        This method allows class instances to be called as functions.

        """

        try:
            self._command()

        except:
            raise vosError('Error in VOs command: {}'
                           ''.format(self._command.__name__))
