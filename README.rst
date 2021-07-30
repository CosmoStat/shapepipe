ShapePipe
=========

|CI| |python38|

.. |CI| image:: https://github.com/CosmoStat/shapepipe/workflows/CI/badge.svg
  :target: https://github.com/CosmoStat/shapepipe/actions?query=workflow%3ACI

.. |python38| image:: https://img.shields.io/badge/python-3.8-green.svg
  :target: https://www.python.org/

:Version: 0.0.4

:Date: 15/07/2021

ShapePipe is a galaxy shape measurement pipeline developed within the
CosmoStat lab at CEA Paris-Saclay.

Contents
========

1. `Dependencies`_

   a. `Core Dependencies`_
   b. `Module Dependencies`_

2. `Installation`_

   a. `Installing ShapePipe`_
   b. `Installing the ShapePipe Library Only`_
   c. `Installing the Module Python Dependencies`_

3. `Execution`_

   a. `Running the Pipeline with SMP`_
   b. `Running the Pipeline with MPI`_
   c. `Configuration`_

4. `Development`_

   a. `Modules`_
   b. `Examples`_

5. `Additional Documentation <docs/wiki/shapepipe.md>`_

Dependencies
============

Core Dependencies
-----------------
- joblib>=0.13.1
- modopt>=1.2.0
- numpy>=1.16.0

Module Dependencies
-------------------

**Python Dependencies**

- astropy>=3.1.1
- GalSim>=2.1.4
- matplotlib>=3.0.2
- numbda==0.40.0
- ngmix>=1.2
- pandas>=1.0
- sf-tools>=2.0.0
- sqlitedict>=1.6.0

**Executable Dependencies**

- CDSclient
- PSFEx
- SExtractor
- WeightWatcher

Installation
============

Installing ShapePipe
--------------------

The entire ShapePipe package, include dependencies, can be built as follows:

.. code-block:: bash

  $ git clone https://github.com/CosmoStat/shapepipe
  $ cd ShapePipe
  $ ./install_shapepipe

The ``install_shapepipe`` script will create the recommended Conda environment
along with all of the required core and module dependencies.

A list installation options can be seen using the ``--help`` option:

.. code-block:: bash

  $ ./install_shapepipe --help


Installing the ShapePipe Library Only
-------------------------------------

The ShapePipe library, *i.e.* the core package not including module dependencies,
can be installed in the following ways.

**Using Git**

The ShapePipe package can be downloaded from the repository
and built as follows:

.. code-block:: bash

  $ git clone https://github.com/CosmoStat/shapepipe
  $ cd shapepipe
  $ python setup.py install


This method is recommend for development.

**Using Pip**

Alternatively, the ShapePipe library can be directly installed from the
repository as follows:

.. code-block:: bash

  $ pip install git+https://github.com/CosmoStat/shapepipe

Note, this method will not include any executable scripts or examples.

Installing the Module Python Dependencies
-----------------------------------------

Module Python dependencies can be installed in the following ways.

**Using Conda**

The ShapePipe Conda environment can be built and activated by running:

.. code-block:: bash

  $ conda env create -f environment.yml
  $ source activate shapepipe

**Using Pip**

Module Python dependencies can also be installed using ``pip`` as follows:

.. code-block:: bash

  $ pip install -r requirements.txt
  $ pip install -r requirements_git.txt

Execution
=========

Running the Pipeline with SMP
-----------------------------

The pipeline can be run with SMP as follows:

.. code-block:: bash

  $ shapepipe_run

A list of command line arguments can be displayed using the ``--help``
option:

.. code-block:: bash

  $ shapepipe_run --help

Running the Pipeline with MPI
-----------------------------

The pipeline can be run with MPI as follows:

.. code-block:: bash

  $ mpiexec -n <number_of_cores> shapepipe_run

where ``<number_of_cores>`` is the number of cores to allocate to the run.

Configuration
-------------

The pipeline requires a configuration file (by default called ``conifg.ini``)
in order to be run. An example configuration file is provided in the
``example`` directory.

The configuration parameters for the pipeline are:

**Default Options**

1. ``VERBOSE`` : (``bool``) Set the verbosity level. Default value is ``True``.
2. ``RUN_NAME`` : (``str``) The pipeline run name. Default value is
   ``shapepipe_run``.
3. ``RUN_DATETIME`` : (``bool``) Option to add date and time to ``RUN_NAME``.
   Default value is ``True``.

**Execution Options**

1. ``MODULE`` : (``str`` or ``list``) A valid module runner name (or a comma
   separated list of names).
2. ``MODE`` : (``str``) The pipeline execution mode. Options are ``smp`` or
   ``mpi``. Default value is ``smp``.

**File Options**

1. ``LOG_NAME`` : (``str``) Current run log file name. Default value is
   ``shapepipe``.
2. ``RUN_LOG_NAME`` : (``str``) Run history log file name. Default value is
   ``shapepipe_runs``.
3. ``INPUT_DIR`` : (``str`` or ``list``) A valid directory containing input
   files for the first module or a comma separated list of directories. This
   parameter also recognises the following special strings:

   a. ``last:MODULE`` : This will point to the output directory of the last run
      of the specified module.
   b. ``all:MODULE`` : This will point to all the output directories in which
      the specified module was run.
   c. ``PATTERN:MODULE`` : This will point to the output directory of a
      specified module from a run matching the specified pattern.

   In all cases the module name can be succeded by the run number (*e.g.*
   ``MODULE_run_2``)

4. ``OUTPUT_DIR`` : (``str``) A valid directory to write the pipeline output
   files.
5. ``FILE_PATTERN`` : (``str`` or ``list``) A list of string patterns to
   identify input files for the first module.
6. ``FILE_EXT`` : (``str`` or ``list``) A list of file extensions to identify
   input files for the first module.
7. ``NUMBERING_SCHEME`` : (``str``) A string indicating the expected numbering
   system for the input files (*e.g.* ``000-0``). Single digits indicate
   integer values without limit, multiples of digits indicate integers with a
   maximum value. Standard characters can be placed around digits (*e.g.*
   ``.``, ``-``, ``:``, *etc.*). Optionally a regular expression can also be
   passed if it is preceded by ``RE:`` (*e.g.* ``RE:-\d{9}``).

**Job Options**

1. ``SMP_BATCH_SIZE`` : (``int``) Number of SMP jobs to run in parallel.
   Default value is ``1``.
2. ``TIMEOUT`` : (``int``) Timeout limit in seconds for a given job.

**Module Options**

Additional module options can be added using the following structure:

.. code-block:: bash

  [MODULE_NAME]
  PARAMETER = PARAMETER VALUE

This mechanism can also be used to modify module decorator properties or append
additional values to list properties as follows:

.. code-block:: bash

  [MODULE_NAME]
  ADD_PARAMETER = PARAMETER VALUE

If a given module is run more than once, run specific parameter values can be
specified as follows:

.. code-block:: bash

  [MODULE_NAME_RUN_X]
  PARAMETER = PARAMETER VALUE

Where ``X`` is an integer greater than or equal to ``1``.


Development
===========

Modules
-------

New modules can be implemented in the pipeline by simply writing a
*module runner*.

The basic requirement for a new module runner is a single function decorated
with the ``module_runner`` wrapper that outputs the module ``stdout`` and
``stderr``. *e.g.*:

.. code-block:: python

  @module_runner(version=0.1)
  def example_module(*args)

    # DO SOMETHING

    return stdout, stderr

The module runner decorator takes the following keyword arguments:

1. ``input_module`` :  (``str`` or ``list``) The name of a preceding module(s)
   whose output provide(s) the input to this module. Default value is ``None``.
2. ``version`` : (``str``) The module version. Default value is ``'0.0'``. The
   module version should always be explicitly declared and be greater than the
   default value.
3. ``file_pattern`` : (``str`` or ``list``) The input file pattern(s) to look
   for. Default value is ``''``.
4. ``file_ext`` : (``str`` or ``list``) The input file extensions(s) to look
   for. Default value is ``''``.
5. ``depends`` : (``str`` or ``list``) The Python package(s) the module depends
   on. Default value is ``[]``.
6. ``executes`` : (``str`` or ``list``) The system executable(s) the module
   implements. Default value is ``[]``.
7. ``numbering_scheme`` : (``str``) The numbering scheme implemented by the
   module to find input files.

The arguments passed to the module runner are the following:

1. ``input_file_list`` : The list of input files.
2. ``run_dirs`` : The dictionary containing module run paths (*e.g.* the output
   path for a given module run is ``run_dirs['output']``).
3. ``file_number_string`` : The number pattern corresponding to the current
   process.
4. ``config`` : The config parser instance, which provides access to the
   configuration file parameter values. Module specific parameters can be passed
   using the following structure:

   .. code-block:: python

     parameter_value = config.get('MODULE_NAME', 'PARAMETER')

5  ``module_config_sec`` : A string specifying the configuration file section
   to be read.
6. ``w_log`` : The worker log instance, which can be used to record additional
   messages in the module output logs using the following structure:

   .. code-block:: python

      w_log.info('MESSAGE')

Examples
--------

The following example module runners are provided in ``shapepipe.modules``.

**Python Example**

In this example a Python script using a ``PythonExample()`` class is implemented.
This module does not read inputs from any preceding module, but looks for files
in the ``INPUT_DIR`` that match the file patterns ``'numbers'`` and
``'letters'`` with file extension ``'.txt'``. This module depends on ``numpy``.

As this module does not implement any system executable, it is not
necessary to return a ``stderr``. Instead any output content that should be
recorded in the log can be returned, otherwise the module runner should simply
return ``None, None``.

.. code-block:: python

  @module_runner(
    version='1.1',
    file_pattern=['numbers', 'letters'],
    file_ext='.txt',
    depends='numpy',
  )
  def python_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
  ):

    output_file_name = (
        f'{run_dirs["output"]}/pyex_output{file_number_string}.cat'
    )
    message = config.get(module_config_sec, 'MESSAGE')

    inst = pe.PythonExample(0)
    inst.read_files(*input_file_list)
    inst.write_file(output_file_name, message)

    return inst.content, None

**Executable Example**

In this example the module runner call the system executable ``head``. This
module read input files from the ``python_example`` module output that match
the file pattern ``'process'`` with file extension ``'.cat'``.

.. code-block:: python

  @module_runner(
    input_module='python_example',
    version='1.0',
    file_pattern='pyex_output',
    file_ext='.cat',
    executes='head',
  )
  def execute_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
  ):

    command_line = f'head {input_file_list[0]}'
    output_file_name = (
        f'{run_dirs["output"]}/head_output{file_number_string}.txt'
    )

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, 'w')
    text_file.write(stdout)

    return stdout, stderr
