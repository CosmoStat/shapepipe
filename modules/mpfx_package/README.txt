mpfx v1.0.4
===========

Multiprocessing framework adaptation layer

Author: Marc Gentile (marc.gentile@epfl.ch)

MPFX is an extension of the base Mpf multiprocessing framework.
 

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.

Source code
-----------

The MPFX code is made of the following files:
- mpfx.cfg: configuration file
- mpfx_SMP.py: main program for the SMP implementation
- mpfx_MPI.py: main program for the MPI implementation
- mpfx_args: command-line argument and options
- mpfx_jobs: parallel job processing
- mpfx_data: dataset management

Prerequisites 
-------------

- python 2.6 or 2.7 on Unix or Mac OS

- mpf v1.0 or greater: base multiprocessing framework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

Installation:
-------------

- Download and unpack the file mpfx-1.0.4.tar.gz in the directory of your
  choice
- MPFX is a pure Python module and can then be installed using the standard
  procedure: 

     1. create a console or terminal session
     2. cd to the directory where mpfx-1.0.4.tar.gz was detared
     3. enter the command: sudo python setup.py install

  If no root access, use: python setup install --home=<target directory> to 
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable. 

Execution
---------

The general syntax is: 

  SMP version: mpfx_SMP.py [options]  
  MPI version: mpfx_MPI.py [options]  

by default, MPFX will look for a configuration file named mpfx.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported 
                      options.   

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file (default: mpfx.cfg)  

To run mpfx_MPI.py, one has to use the mpirun executable as:

   mpirun -np N mpfx_MPI.py [options]

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of MPFX uses one processor for itself, workers will 
share N-1 processors. 

For example: 'mpirun -np 6 mpfx_MPI -d mydir -c myconfig.cfg' will run MPFX
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Release notes
-------------
- v1.0.0: initial stable version 
- v1.0.2: improvement of dataset methods
- v1.0.4: now allows '-' characters in image and catalog filenames

