# -*- coding: utf-8 -*-

"""INFO

Set the package information

:Authors: Samuel Farrens

:Date: 09/04/2018

"""

__version__ = "1.0.0"
__whoami__ = "shapepipe_base"

# Python dependencies
pipe_depend = ['mpfg', 'mpfx', 'scatalog', 'sconfig', 'slogger']
external_depend = ['numpy']

__python_depend__ = pipe_depend + external_depend

# System dependencies

__system_depend__ = ['diff']
