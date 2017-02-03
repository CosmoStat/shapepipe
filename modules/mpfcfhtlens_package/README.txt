Mpfcfhtlens v2.0.0
===================

Multiprocessing framework to be used as a pipeline on the CFHTLenS data.

Author: Marc Gentile (marc.gentile@epfl.ch), Martin Kilbinger (martin.kilbinger@cea.fr)

Mpfcfhtlens is an extension of the base Mpf multiprocessing framework.
 

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.

Source code
-----------

The mpfcfhtlens code is made of the following files:
- mpfgcfhtlens.cfg: configuration file
- mpfcfhtlens.py: main program for the SMP implementation
- mpfcfhtlens.py: main program for the MPI implementation
- mpfcfhtlens: command-line argument and options
- mpfcfhtlens: parallel job processing
- mpfcfhtlens: CFHTLenS dataset management

Prerequisites 
-------------

- python 2.6 or 2.7 on Unix or Mac OS

- mpfg v1.x: base multiprocessing framwework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

Installation:
-------------

- Download and unpack the file mpfcfhtlens-0.5.0.tar.gz in the directory of your
  choice
- mpfcfhtlens is a pure Python module and can then be installed using the standard
  procedure: 

     1. create a console or terminal session
     2. cd to the directory where mpfcfhtlens-0.5.0.tar.gz was detared
     3. enter the command: python setup.py install

  If no root access, use: sudo python setup install --home=<target directory> to 
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable. 

Execution
---------

The general syntax is: 

  SMP version: mpfcfhtlens.py [options]  
  MPI version: mpfcfhtlens.py [options]  

by default, mpfcfhtlens will look for a configuration file named mpfcfhtlens.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported 
                      options.   

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file (default: mpfcfhtlens.cfg)  

To run mpfcfhtlens.py, one has to use the mpirun executable as:

   mpirun -np N mpfcfhtlens.py [options]

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of mpfcfhtlens uses one processor for itself, workers will 
share N-1 processors. 

For example: 'mpirun -np 6 mpfcfhtlens -d mydir -c myconfig.cfg' will run mpfcfhtlens
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Release notes
-------------



