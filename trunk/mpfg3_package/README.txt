Mpfg3 v2.0.0
============

Multiprocessing framework to be used as a pipeline on a GREAT3 Challenge dataset.

Author: Marc Gentile (marc.gentile@epfl.ch)

Mpfg3 is an extension of the base Mpf multiprocessing framework.
 

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.

Source code
-----------

The mpfg3 code is made of the following files:
- mpfg3.cfg: configuration file
- mpfg3_SMP.py: main program for the SMP implementation
- mpfg3_MPI.py: main program for the MPI implementation
- mpfg3_args: command-line argument and options
- mpfg3_jobs: parallel job processing
- mpfg3_data: GREAT3 dataset management

Prerequisites 
-------------

- python 2.6 or 2.7 on Unix or Mac OS

- mpfg v1.x: base multiprocessing framwework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

Installation:
-------------

- Download and unpack the file mpfg3-0.5.0.tar.gz in the directory of your
  choice
- mpfg3 is a pure Python module and can then be installed using the standard
  procedure: 

     1. create a console or terminal session
     2. cd to the directory where mpfg3-0.5.0.tar.gz was detared
     3. enter the command: python setup.py install

  If no root access, use: sudo python setup install --home=<target directory> to 
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable. 

Execution
---------

The general syntax is: 

  SMP version: mpfg3_SMP.py [options]  
  MPI version: mpfg3_MPI.py [options]  

by default, mpfg3 will look for a configuration file named mpfg3.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported 
                      options.   

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file (default: mpfg3.cfg)  

To run mpfg3_MPI.py, one has to use the mpirun executable as:

   mpirun -np N mpfg3_MPI.py [options]

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of mpfg3 uses one processor for itself, workers will 
share N-1 processors. 

For example: 'mpirun -np 6 mpfg3_MPI -d mydir -c myconfig.cfg' will run mpfg3
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Release notes
-------------



