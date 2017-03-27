"""!
@mainpage Multiprocessing framework - Adaptation layer

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

MPFX is an extension of the base Mpf multiprocessing framework. 


@section code Source code

Default configuration file:
- mpfx.cfg: configuration file

Executable files:
- mpfx_SMP.py: main program for the SMP implementation
- mpfx_MPI.py: main program for the MPI implementation

Modules:
- mpfx_args: command-line argument and options
- mpfx_job: parallel job processing
- mpfx_data: @c dataset management

The Python program files are located in the ./mpfx/ directory and the default configuration file
@c mpfx.cfg reside in the ./mpfx/config directory.

@section  Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- mpfg v1.0 or greater: base multiprocessing framework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

@section  Installation

- Download and unpack the file mpfx-1.0.4.tar.gz in the directory of your
  choice
- mpfx is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where mpfx-1.0.4.tar.gz was detared </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
     @endcode to install the module to target_directory. This directory must be
     included in the @c PYTHONPATH environment variable </li>
</ol> 

@section Execution

The general syntax is: 

  - SMP version: <code>mpfx_SMP.py [options]</code>  
  - MPI version: <code>mpfx_MPI.py [options]</code>  

by default, MPFX will look for a configuration file named mpfx.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file (default: mpfx.cfg) </li> 
</ul> 

To run mpfx_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N mpfx_MPI.py [options]</code>

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of MPFX uses one processor for itself, workers will 
share N-1 processors. 

For example: <code>mpirun -np 6 mpfx_MPI -d mydir -c myconfig.cfg</code> will run mpfx
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@note Before running MPFX, edit the configuration file to set the @c BASE_INPUT_DIR value, which
      should point to the directory where the @c input files reside. 
@note In the configuration file, any environment variable prefixed with @c $ such as @c $HOME will
      be expanded when read.

@section notes Release notes

- v1.0.0: initial stable version 
- v1.0.2: improvement of dataset methods
- v1.0.4: now allows '-' characters in image and catalog filenames

"""

from mpfx import *
from mpfx.mpfx_version import __version__
