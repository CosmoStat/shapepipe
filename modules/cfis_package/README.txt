cfis v1.0.0
===================

!!! [Parallel] execution of CFIS field selection, book keeping, coordinate management, etc.

Documentation
-------------

A Doxygen reference documentation is available in the 'doc' directory.

Prerequisites
------------

- Python 2.6 or 2.7 on Unix or Mac OS
- mpfg v1.x: base multiprocessing framework module
- mpfx v1.x: extension layer of multiprocessing framework module
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module
- scatalog v1.x: catalog management

- numpy
- pylab
- astropy


Installation:
-------------

- Download and unpack the file cfis-1.0.0.tar.gz in the directory of your
  choice
- mpfcs82 is a pure Python module and can then be installed using the standard
  procedure:

     1. create a console or terminal session
     2. cd to the directory where cfis-1.0.0.tar.gz was unpacked
     3. enter the command: sudo python setup.py install

  If no root access, use: sudo python setup install --home=<target directory> to
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable.

Execution:
----------

The general syntax is:

  SMP version: cfis_SMP.py [options]
  MPI version: cfis_MPI.py [options]

by default, cfis will look for a configuration file named cfis.cfg in the
./cfis_package/config directory. The location of this file can be changed using the -c and
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported
                      options.

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file

To run cfis_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N cfis_MPI.py [options]</code>

where @c N is the number of processors (or nodes) to use. MPI must has been
installed and configured appropriately.

Since the manager process of Quadg3 uses one processor for itself, workers will
share N-1 processors.

For example: 'mpirun -np 6 cfis_MPI -d mydir -c myconfig.cfg' will run Quadg3
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Before running cfis, edit the configuration file to set the BASE_INPUT_DIR value, which
should point to the directory where the input files reside.
In the configuration file, any environment variable prefixed with '$' (such as $HOME) will
be expanded when read.
