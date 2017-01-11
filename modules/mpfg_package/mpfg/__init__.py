"""!
@mainpage MPF - Multiprocessing framework

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

MPF is a simple multiprocessing frawework for executing code on both SMP 
(Symmetric Multiprocessing) or on serveurs supporting the Message Passing Interface (MPI).

The implemented multiprocessing scheme follows a Producer-Consumer model 
whereby a "manager" process is responsible for producing and dispatching "jobs"
to a pool of "worker" processes. The task of each worker is to process the 
job assigned to it and once done, supply the result to the manager. The manager
can process that result and provide the worker with a new job and so on until 
no work is left.

This scheme works well when the same processing must be applied to a great 
number of files, e.g. thousands of astronomical images to be processed.

Every time it is run, MPF will create a new time-stamped run_xxx directory in
the base output directory specified in the configuration file [DIR.OUTPUT] 
section.  

In MPF, a Job object is created for each file found to process in a base
input directory specified by the user. Such a file may be for example an image
or a catalog to process. All the files under the base input directory 
constitute a 'Dataset' in the sense of MPF. A Dataset can be extended to 
represent a more complex directory structure.

As it is, MPF is only a skeleton that implements a default behavior. One has
to extend its classes and overrides some of its methods in order to implement 
concrete tasks, such as deconvolution or denoising in astronomy. 

One thus has to create a python module based on MPF, i.e. provide classes that 
inherits from relevant MPF classes and overrides some of their methods.

Assuming such a module has been created, it will have to provide:

  1. A subclass of class Args in order to tell where to find the 
   default configuration file and optionnally support additional command-line 
   arguments or options

  2. A subclass of class JobProcessor where the following methods must be 
     overridden to implement custom behavior.
  
   - process_job(): where the specific code of the task to perform must be
                    invoked from.

   - process_job_result(): which defines what to do with the results of the 
                           job, like creating catalogs or making plots. 

   - create_job(): factory method for creating a job object if the base 
                   Job class must be extended. In such a case, the derived
                   Job class must be defined and its constructor invoked there.

   - Optionally, it may be useful to override the following methods:
     - preprocess_job() : processing to be done before process_job() is called   
     - postprocess_job(): processing to be done afer process_job() is called   
     - all_jobs_processed(): invoked once all the jobs have been processed.

Modifying other methods has to be made with great care.   

@section code Source code

The MPF code is made of the following files:
- mpf.cfg: configuration file
- mp_SMP.py: main program for the SMP implementation
- mp_MPI.py: main program for the MPI implementation
- mp_calc_SMP.py: SMP multiprocessing classes
- mp_calc_MPI.py: MPI multiprocessing classes
- mp_calc.py: common code for mp_SMP.py and mp_MPI.py
- mp_args: command-line argument and options
- mp_job: parallel job processing
- mp_data: data management in the form of datasets
- mp_helper: utility functions

@section Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- sconfig v1.x: configuration file parser module
- slogger v1.x: logging module

@section Installation

Download and unpack the file mpfg-1.0.2.tar.gz in the directory of your
choice. MPF is a pure Python module and can then be installed using the standard
procedure: 

<ol>
<li> create a console or terminal session </li>
<li> cd to the directory where mpfg-1.0.2.tar.gz was detared </li>
<li> enter the command: command: @code sudo python setup.py install @endcode </li>

<li> if no root access, use: @code python setup install --home=target_directory 
  @endcode to install the module to target_directory. This directory must be
  included in the @c PYTHONPATH environment variable </li>
</ol> 

@section Execution

The general syntax is: 

  - SMP version: <code>mp_SMP.py [options] </code>  
  - MPI version: <code>mp_MPI.py [options]  </code> 

by default, MPF will look for a configuration file named mpfg.cfg in the 
./config directory. The location of this file can be changed using the -c and 
-d options (see below).

The supported options are:

<ul>
  <li> -h, --help     Display the usage syntax and the list of supported
                      options.  </li>

  <li> -d, --config-dir    Directory where to find the configuration file </li>

  <li> -c, --config-file   Name of the configuration file (default: mpfg.cfg) </li> 
</ul>

To run mp_MPI.py, one has to use the mpirun executable as:

   - <code>mpirun -np N mp_MPI.py [options] </code>

where N is the number of processors (or nodes) to use. MPI must has been 
installed and configured appropriately.
   
Since the manager process of MPF uses one processor for itself, workers will 
share N-1 processors. 

For example: <code>mpirun -np 6 mp_MPI -d mydir -c myconfig.cfg</code> will run mpf
on 6 processors, 1 for the manager and 5 for the workers. The manager and each
of the workers will look for a configuration file: ./mydir/myconfig.cfg.

@section notes Release notes

- v1.0.0: initial version

"""

from mpfg import *
from mpfg.mp_version import __version__ 
