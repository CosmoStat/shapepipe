# -*- coding: utf-8 -*-

"""INFO

Set the package information

:Authors: Samuel Farrens

:Date: 16/04/18
"""

__version__ = "1.0.2"
__whoami__ = "PSFExRun"

# Python dependencies
pipe_depend = ['mpfg', 'mpfx', 'scatalog', 'sconfig', 'shapepipe_base',
               'slogger']
external_depend = ['numpy']

__python_depend__ = pipe_depend + external_depend

# System dependencies

__system_depend__ = ['psfex']
