ShapePipe
=========

|CI| |CD| |python312| |release|

.. |CI| image:: https://github.com/CosmoStat/shapepipe/workflows/CI/badge.svg
  :target: https://github.com/CosmoStat/shapepipe/actions?query=workflow%3ACI

.. |CD| image:: https://github.com/CosmoStat/shapepipe/actions/workflows/pages/pages-build-deployment/badge.svg
  :target: https://github.com/CosmoStat/shapepipe/actions/workflows/pages/pages-build-deployment

.. |python312| image:: https://img.shields.io/badge/python-3.12-green.svg
  :target: https://www.python.org/

.. |release| image:: https://img.shields.io/github/v/release/CosmoStat/shapepipe
  :target: https://github.com/CosmoStat/shapepipe/releases/latest

ShapePipe is a galaxy shape measurement pipeline developed within the
CosmoStat lab at CEA Paris-Saclay. It runs the full chain from raw survey
images to calibrated shear catalogues — object detection, PSF modelling, and
shape measurement — and produced the first UNIONS cosmic-shear release.

Quickstart
----------

ShapePipe ships as a container image, so you can run the bundled example
pipeline — a single CFIS tile through the full chain — without installing
anything:

.. code-block:: bash

    # Apptainer (HPC, no root needed):
    apptainer exec docker://ghcr.io/cosmostat/shapepipe:develop-runtime shapepipe_run_example

    # ...or Docker:
    docker run --rm ghcr.io/cosmostat/shapepipe:develop-runtime shapepipe_run_example

The image is published on every push to the `GitHub Container Registry
<https://github.com/CosmoStat/shapepipe/pkgs/container/shapepipe>`_:
``:develop`` tracks the integration branch, release tags (e.g. ``:v1.1.0``) a
stable cut, and the ``-runtime`` suffix selects the slim batch image over the
full interactive one.

Documentation
-------------

Full documentation lives at https://cosmostat.github.io/shapepipe. Good places
to start:

- `Installation <https://cosmostat.github.io/shapepipe/installation.html>`_ — getting ShapePipe onto your machine or cluster.
- `Basic execution <https://cosmostat.github.io/shapepipe/basic_execution.html>`_ and `configuration <https://cosmostat.github.io/shapepipe/configuration.html>`_ — running ``shapepipe_run`` and writing pipeline configs.
- `Container workflow <https://cosmostat.github.io/shapepipe/container.html>`_ — the two image targets and the ``pyproject.toml`` / ``uv.lock`` / ``Dockerfile`` layers.
- `Running on a cluster <https://cosmostat.github.io/shapepipe/clusters.html>`_ — pulling the image and submitting jobs, with a worked candide (SLURM) example.
- `The Snakemake workflow <https://cosmostat.github.io/shapepipe/workflow.html>`_ — the production orchestration for a whole tile list, driven by ``workflow/bin/sp``.

If you use ShapePipe in academic work, please cite Guinot et al. (2022) and
Farrens et al. (2022).
