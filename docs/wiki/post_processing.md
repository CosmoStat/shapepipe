# Post-processing

This page shows all required steps of post-processing a `ShapePipe` run. Some details pertain specifically to
runs carried out on [canfar](./canfar.md), but most are general.

1. Retrieve files with `ShapePipe` results

   For a local run on the same machine as for post-processing, nothing needs to be done.
   In some cases, the run was carried out on a remote machine or cluster, and the resulting `ShapePipe`
   output files need to be retrieved.
   
   In the specific case of canfar_avail_results.py, this is done as follows.
   
   A. Check availability of results

      A `canfar` job can submit a large number of tiles, whose processing time can vary a lot.
      We assume that the submitted tile ID list is available locally via the ascii file `tile_numbers.txt`. 
      To check which tiles have finished running, and whose results have been uploaded, use
      ```bash
      canfar_avail_results -i tile_numbers.txt -v
      ```
      See `-h` for options.

   B. Download results

      All results files will be downloaded with
      ```bash
      canfar_download_results
      ```
      Use the same options as for same as for `canfar_avail_results`.
      
      This command can be run in the same directory at subsequent times, to complete an ongoing run: Only newer files will be downloaded
      from the `vos` directory.
      
   C. Un-tar results
     ```bash
      untar_results
      ```
      On success, `ShapePipe` output `fits` and `log` files will be now in various subdirs of the `output` directory.

At this step all required `ShapePipe` resulting output files are available in `.`.

2. Optional: Split output in sub-samples

   An optional intermediate step is to create directories for sub-samples, for example one directory
   for each patch on the sky. This will create symbolic links to the results `.tgz` files downloaded in
   the previous step. For example, to create the subdir `tiles_W3` with links to result files to `all` for
   those tiles contained in the list `tiles_W3.txt`, do:
   ```bash
    $SP_ROOT/scripts/python/create_sample_results.py --input_IDs tiles_W3.txt -i . all -o tiles_W3 -v
    ```
    The following steps will then be done in the directory `tiles_W3`.

3. Run PSF diagnostics, create merged catalogue

   Type
   ```bash
   post_proc_sp
   ```
   to automatically perform a number of post-processing steps. Chose the PSF model with the option `-p psfex|mccd`. In detail, these are (and can also be
   done individually
   by hand):
   
   A. Analyse psf validation files
   
      ```bash
      psf_residuals
      ```
      with options as for `post_proc_sp`.
      This script identifies all psf validation files (from all processed tiles downloaded to `pwd`), creates symbolic links,
      merges the catalogues, and creates plots of PSF ellipticity, size, and residuals over the focal plane.

   B. Prepare output directory
   
      Create links to all 'final_cat' result files with 
      ```bash
      prepare_tiles_for_final
      ```
      The corresponding output directory that is created is `output/run_sp_combined/make_catalog_runner/output`.
      On success, it contains links to all `final_cat` output catalogues

   C. Merge final output files
   
      Create a single main catalog:
      ```bash
      merge_final_cat -i <input_dir> -p <param_file> -v
      ```
      Chose as input directory `input_dir` the output of step B. A default parameter file <param_file> is `/path/to/shapepipe/example/cfis/final_cat.param`. 
      On success, the file `./final_cat.npy` is created. Depending on the number of input tiles, this file can be several tens of Gb large. 
