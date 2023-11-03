# Post-processing

This page shows all required steps of post-processing the results from one or
more `ShapePipe` runs. Post-processing combines various individual `ShapePipe`
output files, and creates joint results, for example combining individual tile
catalogues in a large sky area. The output of post-processing is a joint _shape
catalogue_, containing all required information to create a calibrated shear
catalogue via _metacalibration_), a joint star catalogue, and PSF diagnostic plots.

Some of the following steps pertain specifically to runs carried out on [canfar](https://www.canfar.net/en),
but most are general.

1. Retrieve `ShapePipe` result files

   For a local run on the same machine as for post-processing, nothing needs to be done.
   In some cases, the run was carried out on a remote machine or cluster, and the resulting `ShapePipe`
   output files need to be retrieved.
   
   In the specific case of canfar_avail_results.py, this is done as follows.
   
   A. Check availability of results

      A `canfar` job can submit a large number of tiles, whose processing time can vary a lot.
      We assume that the submitted tile ID list is available locally via the ascii file `tile_numbers.txt`. 
      To check which tiles have finished running, and whose results have been uploaded, use
      ```bash
      canfar_avail_results -i tile_numbers.txt -v -p PSF --input_path INPUT_PATH
      ```
      where PSF is one in [`psfex`|`mccd`], and INPUT_PATH the input path on vos, default `vos:cfis/cosmostat/kilbinger/results`.
      See `-h` for all options.

   B. Download results

      All results files will be downloaded with
      ```bash
      canfar_download_results -i tile_numbers.txt -v -p PSF --input_vos INPUT_VOS
      ```
      Use the same options as for same as for `canfar_avail_results`.
      
      This command can be run in the same directory at subsequent times, to complete an ongoing run: Only newer files will be downloaded
      from the `vos` directory. This also assures that partially downloaded or corrupt files will be replaced.

      Checking the `vos` directorty can be slow for large patches.
      To only download files that are not yet present locally (in `.`), first write the missing ones to an ascii file, using again the
      script `canfar_avail_results`, but this time with `.` as input path:
      ```bash
      canfar_avail_results -i tile_numbers.txt --input_path . -p PSF -v -o missing.txt
      '''
      Then, download only the missing files with
      ```bash
      canfar_download_results -i missing.txt --input_vos cosmostat/kilbinger/results_mccd_oc2 -p mccd -v
      ```

   C. Un-tar results
     ```bash
      untar_results -p PSF
      ```
      On success, `ShapePipe` output `fits` and `log` files will be now in various subdirs of the `output` directory.

At this step all required `ShapePipe` resulting output files are available in the current working directory.

2. Optional: Split output in sub-samples

   An optional intermediate step is to create directories for sub-samples, for example one directory
   for each patch on the sky. This will create symbolic links to the results `.tgz` files downloaded in
   the previous step. For example, to create the subdir `tiles_W3` with links to result files to `all` for
   those tiles contained in the list `tiles_W3.txt`, do:
   ```bash
    create_sample_results --input_IDs tiles_W3.txt -i . all -o tiles_W3 -v
    ```
    The following steps will then be done in the directory `tiles_W3`.

3. Run PSF diagnostics, create merged catalogue

   Type
   ```bash
   post_proc_sp -p PSF
   ```
   to automatically perform a number of post-processing steps. Chose the PSF model with the option
   `-p psfex|mccd`. In detail, these are (and can also be done individually
   by hand):
   
   A. Analyse psf validation files
   
      ```bash
      prepare_star_cat -p PSF
      ```
      with options as for `post_proc_sp`.
      This script identifies all psf validation files (from all processed tiles downloaded to `pwd`), creates symbolic links,
      merges the catalogues, and creates plots of PSF ellipticity, size, and residuals over the focal plane.

   B. Create plots of the PSF and their residuals in the focal plane, as a diagnostic of the overall PSF model.
     As a scale-dependend test, which propagates directly to the shear correlation function, the rho statistics are computed,
     see {cite:p}`rowe:10` and {cite:p}`jarvis:16`,
      ```bash
      shapepipe_run -c /path/to/shapepipe/example/cfis/config_MsPl_PSF.ini
      ``` 

   C. Prepare output directory
   
      Create links to all 'final_cat' result files with 
      ```bash
      prepare_tiles_for_final
      ```
      The corresponding output directory that is created is `output/run_sp_combined/make_catalog_runner/output`.
      On success, it contains links to all `final_cat` output catalogues

   D. Merge final output files
   
      Create a single main shape catalog:
      ```bash
      merge_final_cat -i <input_dir> -p <param_file> -v
      ```
      Choose as input directory `input_dir` the output of step C. A default
      parameter file `<param_file>` is `/path/to/shapepipe/example/cfis/final_cat.param`. 
      On success, the file `./final_cat.npy` is created. Depending on the number of
      input tiles, this file can be several tens of Gb large. 
