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
   Module `find_exposures`.

.. code-block:: bash
  mkdir -p output_tiles/find_exposures
  ~/ShapePipe/shapepipe_run.py -c config_tiles/config.find_exposures.ini

On input, original tile images are read (their FITS header), and the images, weight, and flag files of the original exposures.

On output, exposure-single-CCD files (images, weights, and flags) are created.


B. Tiles processing
^^^^^^^^^^^^^^^^^^^

1. Mask images
   Module `mask`.

.. code-block:: bash
  mkdir -p output_tiles/mask
  ~/ShapePipe/shapepipe_run.py -c config_tiles/config.mask.ini

On input, the original images and weights are used.

On output flag files `flag_*.fits` are created.

2. Detect objects
   Module `sextractor`.

.. code-block:: bash
  mkdir -p output_tiles/SExtractor
  ~/ShapePipe/shapepipe_run.py -c config_tiles/config.sex.ini

On input, the original images and weights, as well as the flag files from the last step (B.1) are read.

On output, SExtractor files `sexcat_*.fits` are created.

3. Write detected tiles obects as exposure-single-CCD catalogue files

.. code-block:: bash
  mkdir -p output_tiles/tileobj_as_exp
  ~/ShapePipe/shapepipe_run.py -c config_tiles/config.tileobj_as_exp.ini

On input, the original tile images (to read their FITS header), the SExtractor catalogues (step B.2), and
the exposure-single-CCD images (to use their WCS header information; from A.1) are used.

On output, exposure-single-CCD catalogues `cat.exp*.fits` are created.

C. Exposure-single-CCD images processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Mask images

.. code-block:: bash
  mkdir -p output_exp/mask
  ~/ShapePipe/shapepipe_run.py -c config_exp/config.mask.ini

On input, the exposure-single-CCD images, weights, and flag files (step A.1) are used.

On output, flag files `mask_*.fits` are created. Note that their base names should be different
from the original flag files.

2. Detect objects

.. code-block:: bash
  mkdir -p output_exp/SExtractor
  ~/ShapePipe/shapepipe_run.py -c config_exp/config.sex.ini

On input, the exposure-single-CCD images and  weights (step A.1), and the exposure-single-CCD flags (C.1) are used.

On output, SExtractor catalogue files `sexcat_*.fits` are created.

3. Select stars

.. code-block:: bash
  mkdir -p output_exp/setools
  ~/ShapePipe/shapepipe_run.py -c config_exp/config.setools.ini

On input, the SExtractor catalogue fies from the previous step (C.2) are used.

On output, star candidate catalogues `star_selection_*.fits` are created.

4. Create PSF model

.. code-block:: bash
  mkdir -p output_exp/PSFEx
  ~/ShapePipe/shapepipe_run.py -c config_exp/config.psfex.ini

On input, the star candidate catalogues from the previous step (C.3) are used.

On output, PSF files `*.psf` are created.
