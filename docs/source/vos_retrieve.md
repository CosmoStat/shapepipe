## Retrieve files from VOspace

This page describes how ShapePipe output files can be retrieved via the Virtual Observatory Space
on canfar. This system was used for the CFIS v0 and v1 runs, and is now obsolete.

1. Retrieve ShapePipe result files 

   For a local run on the same machine as for post-processing, nothing needs to be done. In some cases, the run was carried out on a remote machine or cluster, and the resulting ShapePipe output files  
  need to be retrieved.

   In the specific case of canfar_avail_results.py, this is done as follows.

   1. Check availability of results  

      A canfar job can submit a large number of tiles, whose processing time can vary a lot. We assume that the submitted tile ID list is available locally via the ascii file tile_numbers.txt. To check 
      which tiles have finished running, and whose results have been uploaded, use 
      ```bash
      canfar_avail_results -i tile_numbers.txt -v -p PSF --input_path INPUT_PATH
      ```
      where PSF is one in [`psfex`|`mccd`], and INPUT_PATH the input path on vos, default `vos:cfis/cosmostat/kilbinger/results`.
      See `-h` for all options.

   2. Download results

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
      ```
      Then, download only the missing files with
      ```bash
      canfar_download_results -i missing.txt --input_vos cosmostat/kilbinger/results_mccd_oc2 -p mccd -v
      ```

   3. Un-tar results
     ```bash
      untar_results -p PSF
      ```
      On success, `ShapePipe` output `fits` and `log` files will be now in various subdirs of the `output` directory.

At this step all required `ShapePipe` resulting output files are available in the current working directory.
