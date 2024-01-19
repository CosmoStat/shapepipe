# Running `shapepipe` on the canfar science portal

## Introduction

## Steps from testing to parallel running

Before starting a batch remote session job on a large number of images (step 5.),
it is recommended to perform some or all of the testing steps (1. - 4.).


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

3. Run the pipeline script to test the processing step(s), on one image.
   This step has to be run in the patch base directory.

   1. First, run in dry mode:
      ```bash
      init_run_exclusive_canfar.sh -j JOB -e ID -p [psfex|mccd] -k [tile|exp] -n
      ```
   2. Next, perform a real run with
      ```bash
      init_run_exclusive_canfar.sh -j JOB -e ID -p [psfex|mccd] -k [tile|exp] -n
      ```

4. Run remote session script to test job submission using docker images, on one image.
   This step has to be run in the patch base directory.
   1. First, run in dry mode=2, to display curl command, with
      ```bash
      curl_canfar_local.sh -j JOB -e ID -p [psfex|mccd] -k [tile|exp] -n 2
      ```

   2. Next, run in dry mode=1, to use curl command without processing:
      ```bash
      curl_canfar_local.sh -j JOB -e ID -p [psfex|mccd] -k [tile|exp] -n 1
      ```
   3. Then, perform a real run, to use curl with processing:
      ```bash
      curl_canfar_local.sh -j JOB -e ID -p [psfex|mccd] -k [tile|exp]
      ```
   
5. Full run: Call remote session script and docker image with collection of images
      ```bash
      curl_canfar_local.sh -j JOB -f path_IDs -p [psfex|mccd] -k [tile|exp]
      ```
      with `path_IDs` being a text file with one image ID per line.

## Monitoring

Monitoring of the currently active remote session can be performed using the session IDs `session_IDs.txt` written by the
remote session script `curl_canfar_local.sh`. In the patch main directory, run
```bash
curl_canfar_monitor.sh logs
```
to print `stdout` of the remotely run pipeline script, and
```bash
curl_canfar_monitor.sh events
```
to display the remotely started docker image status.

## Post-hoc summary

In the patch main directory, run
```bash
summary_run PATCH
```
to print a summary with missing image IDs per job and module.
