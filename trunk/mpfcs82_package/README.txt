mpfcs82 v1.0.0
==============

Multiprocessing framework to be used for the CFHT Stripe 82 Survey (CS82)

Author: Marc Gentile (marc.gentile@epfl.ch)

mpfcs82 is an extension of the base Mpf multiprocessing framework.
 

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.

Source code
-----------

The mpfcs82 code is made of the following files:
- mpfcs82.cfg: configuration file
- mpfcs82_SMP.py: main program for the SMP implementation
- mpfcs82_MPI.py: main program for the MPI implementation
- mpfcs82_args: command-line argument and options
- mpfcs82_jobs: parallel job processing
- mpfcs82_data: GREAT3 dataset management

Prerequisites 
-------------

- python 2.6 or 2.7 on Unix or Mac OS

- mpfg v1.x: base multiprocessing framwework module
- sconfig v1.x or greater: configuration file parser module
- slogger v1.x or greater: logging module

Installation:
-------------

- Download and unpack the file mpfcs82-2.0.0.tar.gz in the directory of your
  choice
- mpfcs82 is a pure Python module and can then be installed using the standard
  procedure: 

     1. create a console or terminal session
     2. cd to the directory where mpfcs82-2.0.0.tar.gz was detared
     3. enter the command: sudo python setup.py install

  If no root access, use: sudo python setup install --home=<target directory> to 
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable. 

Execution
---------

The general syntax is: 

  SMP version: mpfcs82_SMP.py [options]  
  MPI version: mpfcs82_MPI.py [options]  

by default, mpfcs82 will look for a configuration file named mpfcs82.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported 
                      options.   

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file (default: mpfcs82.cfg)  

To run mpfcs82_MPI.py, one has to use the mpirun executable as:

   mpirun -np N mpfcs82_MPI.py [options]

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of mpfcs82 uses one processor for itself, workers will 
share N-1 processors. 

For example: 'mpirun -np 6 mpfcs82_MPI -d mydir -c myconfig.cfg' will run mpfcs82
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Before running mpfcs82, edit the configuration file to set the BASE_INPUT_DIR value, which
should point to the directory where the input files reside. 
In the configuration file, any environment variable prefixed with '$' (such as $HOME) will
be expanded when read.

