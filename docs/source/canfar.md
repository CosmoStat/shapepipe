
Steps from testing to parallel running:

1. In ID directory, run shapepipe with config file

2. In ID directory, run job

   job_sp_canfar -j JOB

3. In base directory, run single image

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
