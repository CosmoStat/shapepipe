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

# Merge all final cats
# (W3: 140GB RAM)
merge_final_cat -i output/run_sp_combined_final/make_catalog_runner/output -p cfis/final_cat.param -v


# Star catalogue
combine_runs.bash  -p $psf -c psf
shapepipe_run -c $SP_CONFIG/config_Ms_$psf.ini
shapepipe_run -c $SP_CONFIG/config_Pl_$psf.ini

# Convert star cat to WCS
## Convert all input validation psf files and create directories par patch
## psf_conv_all/P?
cd ../star_cat
 ~/shapepipe/scripts/python/convert_psf_pix2world.py -i ../ -P $patch -v

# Combine previously created files within one SP run dir
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
