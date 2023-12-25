# Post-processing

This page shows all required steps of post-processing the results from one or
more `ShapePipe` runs. Post-processing combines various individual `ShapePipe`
output files, and creates joint results, for example combining individual tile
catalogues into a large sky area. The output of post-processing is a joint _shape
catalogue_, containing all required information to create a calibrated shear
catalogue via _metacalibration_), a joint star catalogue, and PSF diagnostic plots.

If main ShapePipe processing happened at the old canfar VM system (e.g. CFIS v0 and v1), go
[here](vos_retrieve.md) for details how to retrieve the ShapePipe output files.

1. Optional: Split output into sub-samples

   An optional intermediate step is to create directories for sub-samples, for example one directory
   for each patch on the sky. This will create symbolic links to the results `.tgz` files downloaded in
   the previous step. For example, to create the subdir `tiles_W3` with links to result files to `all` for
   those tiles contained in the list `tiles_W3.txt`, do:
   ```bash
    create_sample_results --input_IDs tiles_W3.txt -i . all -o tiles_W3 -v
    ```
    The following steps will then be done in the directory `tiles_W3`.

2. Run PSF diagnostics, create merged catalogue

   Type
   ```bash
   post_proc_sp -p PSF
   ```
   to automatically perform a number of post-processing steps. Choose the PSF model with the option
   `-p psfex|mccd`. In detail, these are (and can also be done individually
   by hand):
   
   1. Analyse psf validation files
   
      ```bash
      combine_runs -t psf -p PSF
      ```
      with options as for `post_proc_sp`.
      This script creates a new combined psf run in the ShapePipe `output` directory, by identifying all psf validation files
      and creating symbolic links. The run log file is updated.

   3. Merge individual psf validation files into one catalogue. Create plots of the PSF and their residuals in the focal plane,
      as a diagnostic of the overall PSF model.
      As a scale-dependend test, which propagates directly to the shear correlation function, the rho statistics are computed,
      see {cite:p}`rowe:10` and {cite:p}`jarvis:16`,
      ```bash
      shapepipe_run -c /path/to/shapepipe/example/cfis/config_MsPl_PSF.ini
      ``` 

   4. Prepare output directory
   
      Create links to all 'final_cat' result files with 
      ```bash
      prepare_tiles_for_final
      ```
      The corresponding output directory that is created is `output/run_sp_combined/make_catalog_runner/output`.
      On success, it contains links to all `final_cat` output catalogues

   5. Merge final output files
   
      Create a single main shape catalog:
      ```bash
      merge_final_cat -i <input_dir> -p <param_file> -v
      ```
      Choose as input directory `input_dir` the output of step C. A default
      parameter file `<param_file>` is `/path/to/shapepipe/example/cfis/final_cat.param`. 
      On success, the file `./final_cat.npy` is created. Depending on the number of
      input tiles, this file can be several tens of Gb large. 
