"""!
@mainpage Multiprocessing Framework - Parallel SEXtractor execution

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

@section code Source code

Default configuration file:
- @c pse.cfg

Executable files:

Modules:

The Python program files are located in the ./pse/ directory and the default configuration file
@c pse.cfg reside in the ./pse-1.0.0/config directory.

@section Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- pyfits to manage files in .FITS format from Python 
  (http://www.stsci.edu/institute/software_hardware/pyfits)
- mpfg v1.x: base multiprocessing framework module
- mpfx v1.x: extension layer of multiprocessing framework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module
- scatalog v1.x: catalog management

- An installation of the SEXtractor tool (http://www.astromatic.net/software/sextractor) is also
  required.

@section Installation

- Download and unpack the file pse-1.0.0.tar.gz in the directory of your
  choice
- pse is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where pse-1.0.0.tar.gz was unpacked </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable </li>
</ol> 

@section Execution

The general syntax is: 

  SMP version: pse_SMP.py [options]  
  MPI version: pse_MPI.py [options]  

by default, pse will look for a configuration file named pse.cfg in the
./pse/config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file </li> 
</ul> 

To run pse_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N pse_MPI.py [options]</code>

where @c N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of Quadg3 uses one processor for itself, workers will 
share N-1 processors. 

For example: 'mpirun -np 6 pse_MPI -d mydir -c myconfig.cfg' will run Quadg3
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@note Before running pse, edit the configuration file to set the @c BASE_INPUT_DIR value, which
      should point to the directory where the @c input files reside. 
@note In the configuration file, any environment variable prefixed with @c $ (such as @c $HOME) will
      be expanded when read.
"""

from pse import *
from pse.pse_version import __version__
