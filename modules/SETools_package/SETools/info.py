# -*- coding: utf-8 -*-

"""INFO

Set the package information

:Authors: Samuel Farrens

:Date: 16/04/18

"""

__version__ = "2.0.1"
__whoami__ = "SETools"

# Python dependencies
pipe_depend = ['mpfg', 'mpfx', 'scatalog', 'sconfig', 'shapepipe_base',
               'slogger']
external_depend = ['numpy', 'matplotlib']

__python_depend__ = pipe_depend + external_depend

# System dependencies

__system_depend__ = ['diff', 'ls']
