"""!
@mainpage Multiprocessing framework - Extension layer for the CFHT Stripe 82 Survey (CS82)

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

Mpfcs82 is an extension of the base Mpf multiprocessing framework. 


@section code Source code

Default configuration file:
- mpfcs82.cfg: configuration file

Executable files:
- mpfcs82_SMP.py: main program for the SMP implementation
- mpfcs82_MPI.py: main program for the MPI implementation

Modules:
- mpfcs82_args: command-line argument and options
- mpfcs82_job: parallel job processing
- mpfcs82_data: @c GREAT3 dataset management

The Python program files are located in the ./mpfcs82/ directory and the default configuration file
@c mpfcs82.cfg reside in the ./mpfcs82/config directory.

@section  Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- mpfg v1.x: base multiprocessing framwework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

@section  Installation

- Download and unpack the file mpfcs82-2.0.0.tar.gz in the directory of your
  choice
- mpfcs82 is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where mpfcs82-2.0.0.tar.gz was detared </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
     @endcode to install the module to target_directory. This directory must be
     included in the @c PYTHONPATH environment variable </li>
</ol> 

@section Execution

The general syntax is: 

  - SMP version: <code>mpfcs82_SMP.py [options]</code>  
  - MPI version: <code>mpfcs82_MPI.py [options]</code>  

by default, Mpfcs82 will look for a configuration file named mpfcs82.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file (default: mpfcs82.cfg) </li> 
</ul> 

To run mpfcs82_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N mpfcs82_MPI.py [options]</code>

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of Mpfcs82 uses one processor for itself, workers will 
share N-1 processors. 

For example: <code>mpirun -np 6 mpfcs82_MPI -d mydir -c myconfig.cfg</code> will run mpfcs82
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@note Before running Mpfcs82, edit the configuration file to set the @c BASE_INPUT_DIR value, which
      should point to the directory where the @c GREAT3 files reside. 
@note In the configuration file, any environment variable prefixed with @c $ such as @c $HOME will
      be expanded when read.

@section notes Release notes

- v2.0.0: initial stable version

"""

from mpfcs82 import *
from mpfcs82.mpfcs82_version import __version__
