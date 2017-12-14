# -*- coding: utf-8 -*-

"""PACKAGE NAME

!!! Don't forget to update the package name !!!

@mainpage Multiprocessing Framework - Parallel Package Name execution

Author: Marc Gentile (marc.gentile@epfl.ch) and Samuel Farrens

@section Description

@section code Source code

Default configuration file:
- @c config.cfg

Executable files:

Modules:

The Python program files are located in the ./package_name/ directory and the
default configuration file
@c config.cfg reside in the ./config directory.

Prerequisites

- Python 2.6 or 2.7 on Unix or Mac OS
- mpfg v1.x: base multiprocessing framework module
- mpfx v1.x: extension layer of multiprocessing framework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module
- scatalog v1.x: catalog management

- An installation of the PSFEx tool
(http://www.astromatic.net/software/psfex) is also required.

Installation

- Download and unpack the file package_name-1.0.0.tar.gz in the directory of
  your choice
- package_name is a pure Python module and can then be installed using the
  standard procedure:

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where package_name-1.0.0.tar.gz was unpacked </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install
        --home=target_directory @endcode to install the module to
        target_directory. This directory must be included in the @c
        PYTHONPATH environment variable </li>
</ol>

@section Execution

The general syntax is:

  SMP version: package_name_SMP.py [options]
  MPI version: package_name_MPI.py [options]

by default, package_name will look for a configuration file named config.cfg in
the ./package_name/config directory. The location of this file can be changed
using the -c and -d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file </li>
</ul>

To run package_name_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N package_name_MPI.py [options]</code>

where @c N is the number of processors (or nodes) to use. MPI must has been
installed and configured appropriately.

Since the manager process of Quadg3 uses one processor for itself, workers will
share N-1 processors.

For example: 'mpirun -np 6 package_name_MPI -d mydir -c myconfig.cfg' will run
Quadg3 on 6 processors, 1 for the manager and 5 for the workers. The manager
and each of the workers will look for a configuration file:
./mydir/myconfig.cfg.

Notes
-----
Before running package_name, edit the configuration file to set the
BASE_INPUT_DIR value, which should point to the directory where the input files
reside.

In the configuration file, any environment variable prefixed with (such as
$HOME) will be expanded when read.

"""

from . import *
from .info import __version__, __whoami__
