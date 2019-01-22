CFIS example configuration
==========================


Set up
------

Create a new directory and link to the CFIS example configuration directories.

.. code-block:: bash
        mkdir run_cfis
        cd run_cfis
        ln -s ~/ShapePipe/example/cfis/config_exp
        ln -s ~/ShapePipe/example/cfis/config_exp

(This assumes that `ShapePipe` git repository has been cloned in `~/`).

Select CFIS images
------------------

By default, the pipeline uses all files in the input directory (config file entry `INPUT_DIR`)
that match the input file pattern. For an area selection, use e.g. `cfis_field_select.py `.

Have the CFIS tiles (image and weight files) and exposures (image, weight, and flag files) available
in corresponding directories.

Run
---

1. Mask tiles

   .. code-block:: bash
        mkdir -p output_tiles/mask
        ~/ShapePipe/shapepipe_run.py -c config_tiles/config.mask.ini

1a. Link to flag files.

   This step will be made obsolete. For the moment, when only one input directory is possible,
   the flag files procuced in the last step have to be linked.

2. Decect objects on tiles using SExtractor

   .. code-block:: bash
        ~/ShapePipe/shapepipe_run.py -c config_tiles/config.sex.ini
        mkdir -p output_tiles/SExtractor

3. Identify exposures for selected tiles, and write all HDUs to FITS files.


