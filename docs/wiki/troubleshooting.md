[Home](./shapepipe.md)

# Troubleshooting

This page lists some problems found while attempting to run the pipeline.

## Cloning Error

- Error when trying to clone the pipeline :
`fatal: unable to access 'https://drf-gitlab.cea.fr/cosmostat/ShapePipe/': server certificate verification failed. CAfile: /etc/ssl/certs/ca-certificates.crt CRLfile: none`
  - (linux only ?) Run the command : `export GIT_SSL_NO_VERIFY=1` for more details [here](https://stackoverflow.com/questions/21181231/server-certificate-verification-failed-cafile-etc-ssl-certs-ca-certificates-c)

## Run-time errors

### All modules

- **OBSOLETE** Wrong file format:

   `*** ERROR *** Process-4 - An error occurred while generating catalog: <img_000-0>  ('NoneType' object has no attribute 'rfind')`

   Possible reasons:

   - Wrong input file type. E.g. a FITS catalogue instead of an image is required, and a certain keyword or element is not found
   - Input files are not found, e.g. the data directory is incorrect, or the file names do not match with the INPUT_FILENAME_FORMATS key (e.g. flag files are requested that do not exist).
   - File names contain '.', the file extension cannot be extracted correctly, and the file name matching fails.

- **OBSOLETE** Unknown error:

   `*** ERROR *** Process-28 - An error occurred while generating catalog: <img_000-0> (list index out of range)`

   Potential reasons:
      * required executable not found (e.g. SEXtractor, findgsc2.2)
      * input file(s) not found. Verify the BASE_DIR and FILE_PATTERNS keywords in the [PRIMARY_DATASET] section, and the INPUT_FILENAME_FORMATS key in the [CODE] section of the package config file, and check for the the correct file types (images, weights, flags, catalogues, ...)

- **OBSOLETE** Unspecified error:
   `*** ERROR *** Process-1 - An error occurred while generating catalog: <img_381-0>`

   This might be due to a file check, and can be ignored in many cases.

- **OBSOLETE** AttributeError:

   `AttributeError: 'NoneType' object has no attribute 'write'`
   If earlier warning was issued:
   `mask_SMP.py *** Warning ***: [Errno 2] No such file or directory: 'output/mask'`
   the directory `ouput` does not exist and needs to be created.

- **OBSOLETE** Index error:

   `tuple index out of range`
   This was due to a bug in `SETools`:_log_output_success, via a call to 'log_error_p' (wrong bracket). Could still exist in other modules.

- **OBSOLETE** Local variable:

   `*** ERROR ***: local variable 'img_num' referenced before assignment`
   This can occur if a file in the input data directory is not the correct format. Make sure that the input file pattern `[PRIMARY_DATASET]:FILE_PATTERNS` match only the desired input files.

### Mask module

- File error:

   `External flag file not found`.

   Flag file needs to be provided (e.g. `CFIS.flag-*.fits`). If this does not exist, e.g. for tiles, set the keyword in `config.mask` to False,
   ```ini
   EF_MAKE = False
   ```

- Time out/findgsc exits with error code:

   `****aclient: too long to connect (TIME_MAX=60)` or
   `findgsc* [...] returned with error code 1`

   This happens on `candide` when the pipeline is run from qsub, but not on the terminal.

- **Not sure whether this can still occur...** Missing flag catalogue:

   `SCatalog *** ERROR ***: file .temp/halo_spike_flag-002-0.fits no found`

   Another error should have occured before, preventing weightwatcher (ww) to run.

### SExtractor module

- **I think this is not relevant any more, Axel?** Range error:

   `DETECT_MINAREA keyword out of range`

   In FITS header IQFINAL might be 0, use @IQFINAL%0.5 or something similar in the `sextractor` config file.

- **Not sure...** File error:

   `*Error*: cannot open None`
   Some required files cannot be found, e.g. flag file produced my mask package.

- Unknown function:

   E.g.: `Unknown function : max`

   This might occur if a keyword used as function argument is not found in image header.


### PSFEx module

- Segmentation fault. Sometimes `psfex` exits with a seg fault and no further error message, even when adding the verbose and debug flags to the command line. This could be due to a duplication of HDUs in the FITS file of star candidates. E.g. this can happen if `SETools` runs more than once with the same output file, which adds additional HDUs. This can happen when in the `SETools` input directory not only `SExtractor` output object files but also e.g. background files are, that pass the file name matching.


### PSFExInterpolate module

- **OBSOLETE** Range error:

    `An error occurred while generating catalog: <img_381-0> (list index out of range)`

    Reason: Possibly wrong catalogue file format, e.g. on SEXtractor format with only one main HDU.

- **OBSOLETE** Buffer error:

    '*** ERROR *** Process-5 - An error occurred while generating catalog: <img_022-0> (expected string or buffer)`

- **OBSOLETE** Run error, process does not finish.

    Reason: unknown

- **OBSOLETE** File error:

    `*** ERROR *** Process-6 - An error occurred while generating catalog: <img_010-0> (expected string or buffer)`

    Possibly missing input file (e.g. object/galaxy FITS file)

- Other errors:

    Potential reasons:  One or both of the input files might not have the correct number of HDU extensions (2 for the PSF, 3 for the objects).

### **OBSOLETE** Pipeline runner `run_pipeline.py`

  - IO error:

     `OSError: [Errno 39] Directory not empty: '/automnt/n08data/mkilbing/astro/Runs/ShapePipe/CFIS/testrun2_W3_1deg/input_tile/SExtractor1/config'`

     The directory can not be deleted if some process is still active, e.g. an editor with an open
     config file.
