# Running `ShapePipe` processing and post-processing pipelines on CANFAR

Documentation to create ShapePipe output products for catalogues v1.x.

## Initial Setup

### CANFAR Login

Login to the canfar system with

```bash
canfar auth login
```

This can be done from at notebook or terminal within the canfar science portal,
or any remote terminal that has the canfar library installed.

Check authentication status with

```bash
canfar auth list
```

If not on "default", run

```bash
canfar auth switch default
```

### Set variables (optional)

Set the current patch in the shell as

```bash
patch=P[1-9]
```

For convenience, the current PSF model can be set as environment variable, e.g.:

```bash
psf="psfex"
```

Allowed are `psfex` and `mccd`.

Setting the terminal title to display the patch can be useful for long jobs, to keep track of which terminal
runs which patch:

```bash
echo -ne "\033]0;$patch\007"
```

### Prepare run directory

First, go to the dedicated directory with

```bash
cd /path/to/version/$patch
```

Next, set links to the tile number list and configuration directory:

```bash
ln -s ~/shapepipe/auxdir/CFIS/tiles_202106/tiles_$patch.txt tile_numbers.txt
ln -s ~/shapepipe/example/cfis
```

Create output and debug log directories

```bash
mkdir -p output
mkdir -p debug
```

Finally, create and link to central image storage directories for tiles and exposures:

```bash
mkdir -p ~/cosmostat/v2/data_tiles/$patch
ln -s ~/cosmostat/v2/data_tiles/$patch data_tiles
mkdir -p ~/cosmostat/v2/data_exp/$patch
ln -s ~/cosmostat/v2/data_tiles/$patch data_exp
```

## `ShapePipe` processing

Now, everything should be ready to start running `ShapePipe` for the weak lensing processing. The following
details all necessary steps.

### Get Images

We first download images, and in a second run create symbolic links with the proper pipeine naming scheme.

#### Download and move tiles

When running the main `ShapePipe` script `shapepipe_run`, the following env variable needs to point
to the current working directory

```bash
export SP_RUN=`pwd`
```

Now we run the first module (`get_images_runner`) to download the tile images together with the weight files.
This run can get interrupted by VOSpace I/O or connection errors. In that case, 
we move new files to the image storage directory, remove the previous (now void of images) run directory,
and update the run log file. We also check the number of previous and new tiles.

```bash
shapepipe_run -c cfis/config_Git_vos.ini
ls -l data_tiles/ | wc
mv -i output/run_sp_Git_*/get_images_runner/output/CFIS.???.???.*fits* data_tiles
ls -l data_tiles/ | wc
rm -rf output/run_sp_Git_*
update_runs_log_file.py
```

Repeat the above block as needed.

### Find Exposures

With all tile images (= stacks) downloaded, we can inquire their headers to identify the exposures that were used
to create the stacks. This call to the pipeline also creates the symbolic links to the downloaded tile images.

```bash
shapepipe_run -c cfis/config_GitFe_symlink.ini
```

(One could also run `Fe` alone.)

### Download and Move Exposures

The last module create exposure lists on output. These are now used to download all exposures. As for the tile downloads,
we have to account for VOSpace errors.

```bash
shapepipe_run -c cfis/config_Gie_vos.ini
mv -i output/run_sp_Gie_*/get_images_runner/output/*.fits*fz data_exp
rm -rf output/run_sp_Gie_*
update_runs_log_file.py
```

Repeat the above by hand, or peform it in an automatic loop:

export SP_RUN=`pwd`
shapepipe_run -c ~/shapepipe/example/cfis/config_Ms_psfex_conv.ini


# Extra stuff

## Delete jobs
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
for ID in `cat session_IDs.txt`; do echo $ID; curl -X DELETE -E $SSL $SESSION/$ID; done

## Run in terminal in parallel (-e needs to be last arg)
cat all.txt | xargs -P 16 -n 1  init_run_exclusive_canfar.sh -j 64 -p psfex -n -e

## Get missing jobs that are not currently running
stats_jobs_canfar.sh
grep -F -v -f jobs_running.txt summary/missing_job_128_ngmix_runner_3.txt > all3.txt
patch="P7"
psf="psfex"

# Terminal title
echo -ne "\033]0;$patch\007"

# Run directory
dir=~/cosmostat/v2/pre_v2/$psf/$patch
cd $dir

# Get tile number list
ln -s ~/shapepipe/auxdir/CFIS/tiles_202106/tiles_$patch.txt tile_numbers.txt


# Get images

## Download and link separately

### Download
### Create and link to central image storage directory
mkdir -p ~/cosmostat/v2/data_tiles/$patch
ln -s ~/cosmostat/v2/data_tiles/$patch data_tiles
mkdir -p ~/cosmostat/v2/data_exp/$patch
ln -s ~/cosmostat/v2/data_tiles/$patch data_exp

### Download and move tiles 
ln -s ~/shapepipe/example/cfis
mkdir -p output
export SP_RUN=`pwd`

shapepipe_run -c cfis/config_Git_vos.ini
ls -l data_tiles/ | wc; mv -i output/run_sp_Git_*/get_images_runner/output/CFIS.???.???.*fits* data_tiles; ls -l data_tiles/ | wc
rm -rf output/run_sp_Git_*; update_runs_log_file.py
# repeat the above block

### Find exposures; this run can be stopped after Fe
shapepipe_run -c cfis/config_GitFe_symlink.ini
# You can also run Fe alone

### Download and move exposures

shapepipe_run -c cfis/config_Gie_vos.ini
mv -i output/run_sp_Gie_*/get_images_runner/output/*.fits*fz data_exp
rm -rf output/run_sp_Gie_*
update_runs_log_file.py
# repeat the above; or:
while true; do shapepipe_run -c cfis/config_Gie_vos.ini; ls -l data_exp/ | wc; mv -i output/run_sp_Gie_*/get_images_runner/output/*.fits*fz data_exp;  ls -l data_exp/ | wc; rm -rf output/run_sp_Gie_*; update_runs_log_file.py; done
# Make sure that after all images are downloaded there is no Gie run. This would
# mess up later modules since last:get_image_runner could point to this run.

### Create links (and re-run Fe, not necessary)
job_sp_canfar.bash -p $psf `cat tile_numbers.txt` -j 1 -r symlink

# Get single-HDU single-exposure IDs file (from missing 32 job) 
~/shapepipe/scripts/python/summary_run.py P$patch [32]

# Mask tiles

## Run repeatedly if necessary
job_sp_canfar.bash -p $psf -n $OMP_NUM_THREADS -j 4

## Combine all runs
combine_runs.bash -c flag_tile

# Tile detection
curl_canfar_local.sh -j 16 -f tile_numbers.txt -p $psf -N $OMP_NUM_THREADS

# Option 0, global split and exp masks: sp_local=0
# Todo: split Uz and SpMh

# For sp_local=- both mh_local (0, 1) are ok
export mh_local=0
#export mh_local=1

## Uncompress weights,  split exposures into single HDUs
job_sp_canfar.bash -p $psf -n $OMP_NUM_THREADS -j 2

# Mask exposures

## Run repeatedly if necessary
job_sp_canfar.bash -p $psf -n $OMP_NUM_THREADS -j 8

# Combine all runs
combine_runs.bash -c flag_exp

# Option 1: sp_local=1, local split and mask exp
export mh_local=1

# Split exposures
curl_canfar_local.sh -j 2 -f all.txt -p $psf -N $OMP_NUM_THREADS

# Mask exposures
curl_canfar_local.sh -j 8 -f all.txt -p $psf -N $OMP_NUM_THREADS

# Exposure detection

cp summary/missing_job_32_sextractor.txt all.txt
curl_canfar_local.sh -j 32 -m $mh_local -f all.txt -p $psf -N $OMP_NUM_THREADS

# Tile preparation
curl_canfar_local.sh -j 64 -f tile_numbers.txt -p $psf -N $OMP_NUM_THREADS

# Tile shape measurement
curl_canfar_local.sh -j 128 -f tile_numbers.txt -p $psf -N 8

# Merge subcatalogues
curl_canfar_local.sh -j 256 -f tile_numbers.txt -p $psf -N 8

# Create final cat
curl_canfar_local.sh -j 512 -f tile_numbers.txt -p $psf -N $OMP_NUM_THREADS
# Run in parallel
cat mc.txt | xargs -I {} -P 16 bash -c 'init_run_exclusive_canfar.sh -p psfex -j 512 -e {} --n_smp 1'

# Combine all final cats in common output dir as links
combine_runs.bash -c final -p psfex

# Merge all final cats per patch
# (W3: 140GB RAM)
# in /path/to/$psf
patchnum=`tr $patch P ''`
create_final_cat.py -m final_cat_$patch.hdf5 -i . -p $patch/cfis/final_cat.param -P $patchnum -o $patch/n_tiles_final.txt -v

# Star catalogue
combine_runs.bash  -p $psf -c psf
shapepipe_run -c $SP_CONFIG/config_Ms_$psf.ini
shapepipe_run -c $SP_CONFIG/config_Pl_$psf.ini

# Convert star cat to WCS
## Convert all input validation psf files and create directories par patch
## psf_conv_all/P?
cd ../star_cat

# Create files validation_psf_conv-<patchnum>-<idx>.fits
# (for the v1.4 setup only one file)
 convert_psf_pix2world.py -i .. -P $patchnum -v

# Combine previously created files as links within one SP run dir
# (for the v1.4 setup only one link
cd P$patch
combine_runs.bash -p psfex -c psf_conv

# Merge all converted star catalogues and create final-starcat.fits
export SP_RUN=`pwd`
shapepipe_run -c ~/shapepipe/example/cfis/config_Ms_psfex_conv.ini


# Extra stuff

## Delete jobs
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
for ID in `cat session_IDs.txt`; do echo $ID; curl -X DELETE -E $SSL $SESSION/$ID; done

## Run in terminal in parallel (-e needs to be last arg)
cat all.txt | xargs -P 16 -n 1  init_run_exclusive_canfar.sh -j 64 -p psfex -n -e

## Get missing jobs that are not currently running
stats_jobs_canfar.sh
grep -F -v -f jobs_running.txt summary/missing_job_128_ngmix_runner_3.txt > all3.txt

