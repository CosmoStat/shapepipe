# Running `shapepipe` on the canfar science portal

## Introduction

## Steps from testing to parallel running

Before starting a batch remote session job on a large number of images (step 6.),
it is recommended to perform some or all of the testing steps (1. - 5.).


1. Run the basic `shapepipe` runner script to test (one or several) modules in question, specified by a given config file, on one image.
   This step has to be run in the image run directory. The command is
   ```bash
   shapepipe_run -c config.ini
   ```

2. Run the job script to test the job management, on one image.
   This step has to be run in the image run directory. The command is
   ```bash
   job_sp_canfar -j JOB [OPTIONS]
   ```

5. In base directory, run single image

  First in dry mode

  init_run_exclusive_canfar.sh -j 64 -e 271.281 -p psfex -k tile -n

  Then for real, remove option -n

4. In base dir, run single image with curl

    First in dry-mode=2 (showing curl command)
    curl_canfar_local.sh -j 64 -e 271.282 -p psfex -k tile -n 2

    First in dry-mode=1 (using curl command)
    curl_canfar_local.sh -j 64 -e 271.282 -p psfex -k tile -n 1

    Then real run
    curl_canfar_local.sh -j 64 -e 271.282 -p psfex -k tile

5. In base dir, run collection of images

    curl_canfar_local.sh -j 64 -f tile_numbers.txt -p psfex -k tile
