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

   An optional intermediate step is to create directories for sub-samples, for example one directory
   for each patch on the sky. This will create symbolic links to the results `.tgz` files downloaded in
   the previous step. For example, to create the subdir `tiles_W3` with links to result files to `all` for
   those tiles contained in the list `tiles_W3.txt`, do:
   ```bash
    $SP_ROOT/scripts/python/create_sample_results.py --input_IDs tiles_W3.txt -i . all -o tiles_W3 -v
    ```
    The following steps will then be done in the directort `tiles_W3`.
    
   Type
   ```bash
   $SP_ROOT/scripts/sh/canfar_post_proc.sh
   ```
   to automatically perform a number of post-processing steps. In detail, these are (and can also be done individually
   by hand):
   
   1. Un-tar all result `.tgz` files in the current directory:
      ```bash
      $SP_ROOT/scripts/sh/untar_results.sh 
      ``` 
      As a result, the directory `output` is created containing the `ShapePipe` outputs from all canfar runs,
      in uniquely named subdirectories.
      
   2. Analyse psf validation files
      ```bash
      $SP_ROOT/scripts/sh/canfar_psf_residuals.sh
      ```
      This script identifies all psf validation files (from all processed tiles downloaded to `pwd`), creates symbolic links,
      merges the catalogues, and creates plots of PSF ellipticity, size, and residuals over the focal plane.

   3. Prepare output directory with links to all 'final_cat' result files:
      ```bash
      $SP_ROOT/scripts/sh/canfar_prep_tiles.sh
      ```

   4. Merge final output files to single mother catalog
      ```bash
      input_final=output/run_sp_combined/make_catalog_runner/output
      $SP_ROOT/scripts/python/merge_final_cat.py -i $input_final -p $SP_CONFIG/final_cat.param -v
      ```
  
