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
that match the input file pattern. For an area selection, use e.g. `cfis_field_select.py`, example, to run
on a machine where the CFIS images are stored:

.. code-block:: bash
        cfis_field_select.py -i /home/mkilbing/astro/data/CFIS/tiles -t tile -m a --plot -v --area 210deg_55deg_211deg_56deg -o area_W3_1deg


Retrieve to local machine from cc
---------------------------------

On your local machine, write the selected image file basenames into the text file (`tiles.txt`), and copy those from the cc:

.. code-block:: bash
        scp_CFIS_cc.py -i tiles.txt --from_cc -t tile -v

Have the CFIS tiles (image and weight files) and exposures (image, weight, and flag files) available
in corresponding directories.

Run
---


A. Preprocessing
^^^^^^^^^^^^^^^^

1. Identify exposures for selected tiles, and write all HDUs to FITS files.

.. code-block:: bash
        mkdir -p output_tiles/find_exposures
        ~/ShapePipe/shapepipe_run.py -c config_tiles/config.find_exposures.ini


B. Tiles processing
^^^^^^^^^^^^^^^^^^^

1. Mask

.. code-block:: bash
        mkdir -p output_tiles/mask
        ~/ShapePipe/shapepipe_run.py -c config_tiles/config.mask.ini

On output flag files `flag_*.fits` are created.

2. Detect objects

.. code-block:: bash
        mkdir -p output_tiles/SExtractor
        ~/ShapePipe/shapepipe_run.py -c config_tiles/config.sex.ini

On input, the original images and weights, as well as the flag files from the last step (B.1) are read.

On output, SExtractor files `sexcat_*.fits` are created.


