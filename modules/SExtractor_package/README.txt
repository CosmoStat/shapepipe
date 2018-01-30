SExtractor v1.0.0
===================

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
- scatalog v2.x: catalog management


Installation:
-------------

- Download and unpack the file SExtractor-1.0.0.tar.gz in the directory of your
  choice
- mpfcs82 is a pure Python module and can then be installed using the standard
  procedure:

     1. create a console or terminal session
     2. cd to the directory where SExtractor-1.0.0.tar.gz was unpacked
     3. enter the command: sudo python setup.py install

  If no root access, use: sudo python setup install --home=<target directory> to
  install the module to <target directory>. This directory must be included in
  the PYTHONPATH environment variable.

Execution:
----------

The general syntax is:

  SMP version: SExtractor_SMP.py [options]
  MPI version: SExtractor_MPI.py [options]

by default, SExtractor will look for a configuration file named SExtractor.cfg in the
./SExtractor/config directory. The location of this file can be changed using the -c and
-d options (see below).

The supported options are:

  -h, --help          Display the usage syntax and the list of supported
                      options.

  -d, --config-dir    Directory where to find the configuration file

  -c, --config-file   Name of the configuration file

To run SExtractor_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N SExtractor_MPI.py [options]</code>

where @c N is the number of processors (or nodes) to use. MPI must has been
installed and configured appropriately.

Since the manager process of Quadg3 uses one processor for itself, workers will
share N-1 processors.

For example: 'mpirun -np 6 SExtractor_MPI -d mydir -c myconfig.cfg' will run Quadg3
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

Before running SExtractor, edit the configuration file to set the BASE_INPUT_DIR value, which
should point to the directory where the input files reside.
In the configuration file, any environment variable prefixed with '$' (such as $HOME) will
be expanded when read.
