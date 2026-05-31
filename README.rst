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
CosmoStat lab at CEA Paris-Saclay.

See the `documentation <https://cosmostat.github.io/shapepipe>`_ for details
on how to install and run ShapePipe.

Quickstart on a cluster (candide)
---------------------------------

ShapePipe ships as a container image — the supported way to run it (see
``docs/source/container.md``). On a SLURM cluster such as candide, pull the slim
``runtime`` image once and submit the bundled example, which runs the pipeline
on a single CFIS tile:

.. code-block:: bash

    # 0. Get a clone (holds the example configs, data, and job scripts).
    git clone https://github.com/CosmoStat/shapepipe.git
    cd shapepipe

    # 1. Keep the SIF and Apptainer's scratch off the quota-limited $HOME.
    #    candide's home quota is tight; a pull there fails with "disk quota
    #    exceeded". Point both at a roomy data partition instead.
    export DATA=/n17data/$USER                 # adjust to your data partition
    export APPTAINER_CACHEDIR=$DATA/.apptainer

    # 2. Pull the runtime image (≈850 MB).
    apptainer pull "$DATA/shapepipe-runtime.sif" \
        docker://ghcr.io/cosmostat/shapepipe:develop-runtime

    # 3. Submit the example pipeline (SMP, single node).
    SP_IMAGE="$DATA/shapepipe-runtime.sif" SPDIR="$PWD" \
        sbatch example/pbs/candide_smp.sh

A clean run logs ``A total of 0 errors were recorded`` and exits ``0``. To span
multiple nodes with hybrid MPI, swap in ``example/pbs/candide_mpi.sh`` (same two
variables) — see the comments in each script for the SLURM directives.

The ``:develop-runtime`` tag tracks the integration branch; for a stable cut use
a release tag (e.g. ``:v1.1.0-runtime``). The interactive ``dev`` image (no
``-runtime`` suffix) carries ``vim``, ``pytest``, and the full toolchain for
working *inside* the container; ``docs/source/container.md`` covers both.
