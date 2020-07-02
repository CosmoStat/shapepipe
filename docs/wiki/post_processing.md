# Post-processing

This page specifically deals with the post-processing step of results obtained
on [canfar](./canfar.md).

1. Canfar-access steps.

   1. Check availability of results

      A `canfar` job can submit a large number of tiles, and many but not all might be processed
      simultaneously. To check which tiles are finished, and whose results have been uploaded, use
      ```bash
      $SP_ROOT/scripts/python/canfar_avail_results.py -i <ID_files> -v
      ```
      E.g  with `-i $SP_ROOT/aux/CFIS/tiles_202007/tiles_W3.txt` as input ID file.

   2. Download results

      All results files will be downloaded with
      ```bash
      $SP_ROOT/scripts/sh/canfar_download_results.sh
      ```
      This command can be run in the same directory at subsequent times: Only newer files will be downloaded
      from the vos directory.


2. Local post-processing steps.

   On [candide](./candide.md) it is advisable to perform the following steps not on the loging node, but
   to log interactively onto some other node. Internet access is not required.

   Type
   ```bash
   $SP_ROOT/scripts/sh/canfar_post_proc.sh
   ```
   to automatically perform a number of post-processing steps. In detail, these are (and can also be done individually
   by hand):
   
   a. Un-tar all result `.tgz` files in the current directory:
      ```bash
      $SP_ROOT/scripts/sh/untar_results.sh 
      ``` 
      As a result, the directory `output` is created containing the `ShapePipe` outputs from all canfar runs,
      in uniquely named subdirectories.
      
   b. Analyse psf validation files
      ```bash
      $SP_ROOT/scripts/sh/canfar_psf_residuals.sh
      ```
      This script identifies all psf validation files (from all processed tiles downloaded to `pwd`), creates symbolic links,
      merges the catalogues, and creates plots of PSF ellipticity, size, and residuals over the focal plane.

   c. Prepare output directory with links to all 'final_cat' result files:
      ```bash
      $SP_ROOT/scripts/sh/canfar_prep_tiles.sh
      ```

   d. Merge final output files to single mother catalog
      ```bash
      input_final=output/run_sp_combined/make_catalog_runner/output
      $SP_ROOT/scripts/python/merge_final_cat.py -i $input_final -p $SP_CONFIG/final_cat.param -v
      ```
  
