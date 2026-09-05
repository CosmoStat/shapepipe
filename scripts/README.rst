Scripts directory
=================

Helpers that operate on pipeline inputs and outputs but are **not** pipeline
modules: nothing here is run by ``shapepipe_run``, and nothing here is part of
the Snakemake workflow's rule graph (that is ``workflow/scripts/``).

Everything under ``scripts/*/`` with a ``.py``, ``.sh`` or ``.bash`` extension is
symlinked onto ``$PATH`` inside the container image under its bare name, so
``create_star_cat`` and ``python scripts/python/create_star_cat.py`` are the same
command.

Each script's own ``--help`` and module docstring is the reference; this file is
only a map. Run any of them with ``-h`` first.

``python/`` — pipeline inputs
-----------------------------

- ``cfis_field_select.py`` — select CFIS tiles or exposures by sky coordinates
  or by name, and emit the ID list. This is how a tile list for
  ``workflow/config.yaml`` is made.
- ``create_star_cat.py`` — build reference star catalogues for the mask module.
  Superseded for exposures by the workflow's ``exp_star_cat`` rule; its
  ``-k tile`` mode is the prerequisite for a future tile-mask rule.
- ``check_tile_coverage.py`` — flag input simulation tiles whose weight maps are
  mostly empty, as a YAML exclude list. Called from sp_validation's image-sims
  orchestration, not from here.

``python/`` — pipeline outputs
------------------------------

- ``collate_star_cat.py`` — collate the per-exposure PSF validation catalogues
  into the input ``merge_starcat`` wants for the rho/tau statistics.
- ``create_final_cat.py`` — merge per-tile final catalogues into one HDF5 file
  for sp_validation.
- ``stats_global.py`` — histograms and tables from the SETools per-CCD summary
  statistics of a whole run.

``python/`` — image simulations
-------------------------------

The image-simulation chain is still driven by the pre-Snakemake bash layer,
which sp_validation's ``image_sims.smk`` calls into. These exist for it and are
not used by real-data runs.

- ``init_run_v2.0.py`` — lay out a run directory for that layer.
- ``update_runs_log_file.py`` — rebuild ``log_run_sp.txt`` from the run
  directories on disk, which that layer needs because it deletes run dirs.
- ``test_tile_det.py`` + ``test_tile_det.cfg`` — drive tile detection over it.

``python/`` — the ngmix status dashboard
----------------------------------------

Orthogonal to everything above: ``build_status.py``, ``build_history.py``,
``plot_trends.py``, ``run_breakdown_grid.py``, ``plot_breakdown_grid.py`` and
``plot_s4_ablations.py`` build the shape-measurement status page and its
trend/ablation figures.

``sh/``
-------

- ``shapepipe_run_example.sh`` — run the bundled example pipeline against a
  writable copy of ``/app/example``. This is the container's CI smoke test.
- ``apptainer_noslurm.sh`` — ``apptainer exec`` with the SLURM/PMI/PMIX/OMPI
  environment stripped, for MPI inside a SLURM job.
- ``run_job_sp_canfar_v2.0.bash``, ``job_sp_canfar_v2.0.bash``,
  ``job_list_help.bash``, ``functions.sh`` — the bit-coded bash job layer. Kept
  only for image simulations (above); real-data runs go through ``workflow/``.

``jupyter/``
------------

Scratch notebooks: ``wcs.ipynb``, ``test_centroid_shift.py``.

``validation/``
---------------

``centroid/`` — the metacalibration multiplicative-bias validation of the ngmix
centroid handling: ``run_all.sh`` drives ``centroid_bias_v2.py`` over the three
cases in ``configs/``.
