patch="P7"
psf="psfex"
N_SMP=16

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

### Download and move tiles 
ln -s ~/shapepipe/example/cfis
mkdir -p output
export SP_RUN=`pwd`

shapepipe_run -c cfis/config_Git_vos.ini
mv -i output/run_sp_Git_*/get_images_runner/output/CFIS.???.???.*fits* data_tiles
rm -rf output/run_sp_tiles_Git_*
update_run_log_file.py
# repeat the above block

### Find exposures; this run can be stopped after Fe
shapepipe_run -c cfis/config_GitFe_symlink.ini

### Download and move exposures

shapepipe_run -c cfis/config_Gie_vos.ini
mv -i output/run_sp_Gie_*/get_images_runner/output/*.fits*fz data_exp
rm -rf  output/run_sp_Gie_*
update_run_log_file.py
# repeat the above

### Create links (and re-run Fe, not necessary)
job_sp_canfar.bash -p $psf `cat tile_numbers.txt` -j 1 -r symlink

# Uncompress weights,  split exposures into single HDUs
job_sp_canfar.bash -p $psf -n $N_SMP -j 2

# Mask tiles
job_sp_canfar.bash -p $psf -n $N_SMP -j 4

# If not finshed:
combine_runs.bash -p psfex -c flag
mv output/run_sp_combined_flag output/run_sp_exp_Ma

# Mask exposures
job_sp_canfar.bash -p $psf -n $N_SMP -j 8


# Tile detection
curl_canfar_local.sh -j 16 -f tile_numbers.txt -k tile -p $psf -N $N_SMP


# Exposure detection
## Get single-HDU single-exposure IDs
~/shapepipe/scripts/python/summary_run.py

cp summary/missing_job_32_sextractor.txt all.txt
curl_canfar_local.sh -j 32 -f all.txt -k exp -p $psf -N $N_SMP

# Tile preparation
curl_canfar_local.sh -j 64 -f tile_numbers.txt -k tile -p $psf -N $N_SMP

# Tile shape measurement
curl_canfar_local.sh -j 128 -f tile_numbers.txt -k tile -p $psf -N 8

# Merge subcatalogues, and create final cat
job_sp_canfar.bash -p $psf -n 1 -j 256

# Combine all final cats in common output dir as links
combine_runs.bash -c final -p psfex

# Merge all final cats
# (use 192GB RAM)
merge_final_cat -i output/run_sp_combined_final/make_catalog_runner/output -p cfis/final_cat.param -v


# Delete jobs
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
for ID in `cat session_IDs.txt`; do echo $ID; curl -X DELETE -E $SSL $SESSION/$ID; done
