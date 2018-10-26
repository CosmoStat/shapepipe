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

   1. `Required Packages`_
   2. `Optional Packages`_

2. `Installation`_
3. `Execution`_

   1. `Running the Pipeline`_
   2. `Configuration`_

4. `Development`_

   1. `Modules`_

Dependencies
============

Required Packages
-----------------
TBD

Optional Packages
-----------------
TBD

Installation
============

The ShapePipe package can be installed by running:

.. code-block:: bash

  $ python setup.py install

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
in order to be run.

An example configuration file is provided in the ``example`` directory.

Development
===========

Modules
-------

New modules can be implemented in the pipeline by simply writing a
*module runner*. Example module runners are provided in ``shapepipe.modules``.

The basic requirement for a new module runner is a single function decorated
with the ``module_runner`` wrapper that outputs the module ``stdout`` and
``stderr``. *e.g.*:

.. code-block:: python

 @module_runner()
  def example_module(*args, **kwargs)

    # DO SOMETHING

    return stdout, stderr
