"""!
@mainpage Multiprocessing framework - Extension layer

Author: Marc Gentile (marc.gentile@epfl.ch), Martin Kilbinger (martin.kilbinger@cea.fr)

@section Description

Mpfx is an extension of the base Mpf multiprocessing framework. 


@section code Source code

Default configuration file:
- mpfcfhtlens.cfg: configuration file

Executable files:
- mpfcfhtlens_SMP.py: main program for the SMP implementation
- mpfcfhtlens_MPI.py: main program for the MPI implementation

Modules:
- mpfcfhtlens_args: command-line argument and options
- mpfcfhtlens_job: parallel job processing
- mpfcfhtlens_data: @c GREAT3 dataset management

The Python program files are located in the ./mpfcfhtlens/ directory and the default configuration file
@c mpfcfhtlens.cfg reside in the ./mpfcfhtlens/config directory.

@section  Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- mpfg v1.x: base multiprocessing framwework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

@section  Installation

- Download and unpack the file mpfcfhtlens-0.0.1.tar.gz in the directory of your
  choice
- mpfcfhtlens is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where mpfcfhtlens-0.0.1.tar.gz was detared </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
     @endcode to install the module to target_directory. This directory must be
     included in the @c PYTHONPATH environment variable </li>
</ol> 

@section Execution

The general syntax is: 

  - SMP version: <code>mpfcfhtlens_SMP.py [options]</code>  
  - MPI version: <code>mpfcfhtlens_MPI.py [options]</code>  

by default, Mpfx will look for a configuration file named mpfcfhtlens.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file (default: mpfcfhtlens.cfg) </li> 
</ul> 

To run mpfcfhtlens_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N mpfcfhtlens_MPI.py [options]</code>

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of Mpfx uses one processor for itself, workers will 
share N-1 processors. 

For example: <code>mpirun -np 6 mpfcfhtlens_MPI -d mydir -c myconfig.cfg</code> will run mpfcfhtlens
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@note Before running Mpfx, edit the configuration file to set the @c BASE_INPUT_DIR value, which
      should point to the directory where the @c GREAT3 files reside. 
@note In the configuration file, any environment variable prefixed with @c $ such as @c $HOME will
      be expanded when read.

@section notes Release notes

- v0.0.1: initial version

"""

from mpfcfhtlens import *
from mpfcfhtlens.mpfcfhtlens_version import __version__
