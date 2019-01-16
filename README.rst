ShapePipe
=========

|gitlab-ci| |python35| |python36|

.. |gitlab-ci| image:: https://drf-gitlab.cea.fr/cosmostat/ShapePipe/badges/shapepipe2_dev/pipeline.svg
  :target: https://drf-gitlab.cea.fr/cosmostat/ShapePipe/tree/shapepipe2_dev

.. |python35| image:: https://img.shields.io/badge/python-3.5-yellow.svg
  :target: https://www.python.org/

.. |python36| image:: https://img.shields.io/badge/python-3.6-yellow.svg
  :target: https://www.python.org/

:Version: 0.0.1

:Date: 26/10/2018

ShapePipe is a galaxy shape measurement pipeline developed within the
CosmoStat lab at CEA Paris-Saclay.

Contents
========

1. `Dependencies`_

   a. `Pipeline Dependencies`_
   b. `Module Dependencies`_

2. `Installation`_

   a. `Conda`_
   b. `Pip`_

3. `Execution`_

   a. `Running the Pipeline`_
   b. `Configuration`_

4. `Development`_

   a. `Modules`_
   b. `Examples`_

Dependencies
============

Pipeline Dependencies
---------------------
- joblib>=0.13.1
- modopt>=1.2.0
- numpy>=1.16.0

Module Dependencies
-------------------
- astropy>=3.1.1
- GalSim
- matplotlib>=3.0.2
- ngmix
- psfex
- pybind11>=2.2
- sf-tools>=2.0.0

Installation
============

The ShapePipe package can be installed by running:

.. code-block:: bash

  $ python setup.py install

Conda
-----

A ShapePipe Conda environment can be built and activated by running:

.. code-block:: bash

  $ conda env create -f environment.yml
  $ source activate shapepipe

Pip
---

The ShapePipe dependencies can also be installed using ``pip`` as follows:

.. code-block:: bash

  $ pip install -r requirements.txt
  $ pip install -r requirements_git.txt

Execution
=========

Running the Pipeline
--------------------

The pipeline can be run as follows:

.. code-block:: bash

  $ ./shapepipe_run

A list of command line arguments can be displayed using the ``--help``
option:

.. code-block:: bash

  $ ./shapepipe_run --help

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

**File Options**

1. ``LOG_NAME`` : (``str``) Current run log file name. Default value is
   ``shapepipe``.
2. ``RUN_LOG_NAME`` : (``str``) Run history log file name. Default value is
   ``shapepipe_runs``.
3. ``INPUT_DIR`` : (``str`` or ``list``) A valid directory containing input
   files for the first module or a comma separated list of directories. This
   parameter also recognizes the following special strings:

   a. ``last:MODULE`` : This will point to the output directory of the last run
      of the specified module.
   b. ``PATTERN:MODULE`` : This will point to the output directory of a
      specified module from a run matching the specified pattern.

4. ``OUTPUT_DIR`` : (``str``) A valid directory to write the pipeline output
   files.
5. ``FILE_PATTERN`` : (``str`` or ``list``) A list of string patterns to
   identify input files for the first module.
6. ``FILE_EXT`` : (``str`` or ``list``) A list of file extensions to identify
   input files for the first module.

**Job Options**

1. ``BATCH_SIZE`` : (``int``) Number of jobs to run in parallel. Default value
   is ``1``.
2. ``TIMEOUT`` : (``int``) Timeout limit in seconds for a given job.

**Module Options**

Additional module options can be added using the following structure:

.. code-block:: bash

   [MODULE_NAME]
    PARAMETER = PARAMETER VALUE

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

  @module_runner()
  def example_module(*args)

    # DO SOMETHING

    return stdout, stderr

The module runner decorator takes the following keyword arguments:

1. ``input_module`` :  (``str`` or ``list``) The name of a preceding module(s)
   whose output provide(s) the input to this module. Default value is ``None``.
2. ``version`` : (``str``) The module version. Default value is ``'0.0'``.
3. ``file_pattern`` : (``str`` or ``list``) The input file pattern(s) to look
   for. Default value is ``''``.
4. ``file_ext`` : (``str`` or ``list``) The input file extensions(s) to look
   for. Default value is ``''``.
5. ``depends`` : (``str`` or ``list``) The Python package(s) the module depends
   on. Default value is ``[]``.
6. ``executes`` : (``str`` or ``list``) The system executable(s) the module
   implements. Default value is ``[]``.

The arguments passed to the module runner are the following:

1. ``input_file_list`` : The list of input files.
2. ``output_dir`` : The directory for the module output files.
3. ``job_name`` : The name given by the pipeline to the current job.
4. ``config`` : The config parser instance, which provides access to the
   configuration file parameter values. Module specific parameters can be passed
   using the following structure:

   .. code-block:: python

     parameter_value = config.get('MODULE_NAME', 'PARAMETER')

5. ``w_log`` : The worker log instance, which can be used to record additional
   messages in the module output logs using the following structure:

   .. code-block:: python

      self.w_log.info('MESSAGE')

Examples
--------

The following example module runners are provided in ``shapepipe.modules``.

**Python Example**

In this example a Python script using a ``Dummy()`` class is implemented. This
module does not read inputs from any preceding module, but looks for files
in the ``INPUT_DIR`` that match the file patterns ``'numbers'`` and
``'letters'`` with file extension ``'.txt'``. This module depends on
``numpy``.

As this module does not implement any system executable, it is not
necessary to return a ``stderr``. Instead any output content that should be
recorded in the log can be returned, otherwise the module runner should simply
return ``None, None``.

.. code-block:: python

  @module_runner(version='1.0', file_pattern=['numbers', 'letters'],
                file_ext='.txt', depends='numpy')
  def python_example(input_file_list, output_dir, job_name, config, w_log):

      output_file_name = ('{}/{}.cat'.format(output_dir, job_name))
      message = config.get('PYTHON_EXAMPLE', 'MESSAGE')

      inst = Dummy()
      inst.read_files(*input_file_list)
      inst.write_file(output_file_name, message)

      return inst.content, None

**Executable Example**

In this example the module runner call the system executable ``head``. This
module read input files from the ``python_example`` module output that match
the file pattern ``'process'`` with file extension ``'.cat'``.

.. code-block:: python

  @module_runner(input_module='python_example', version='1.0',
                file_pattern='process', file_ext='.cat', executes='head')
  def execute_example(input_file_list, output_dir, job_name, *args):

      command_line = 'head {}'.format(input_file_list[0])
      output_file_name = '{}/{}_head_out.txt'.format(output_dir, job_name)

      stdout, stderr = execute(command_line)

      text_file = open(output_file_name, 'w')
      text_file.write(stdout)

      return stdout, stderr
