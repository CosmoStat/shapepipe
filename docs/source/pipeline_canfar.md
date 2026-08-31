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

```bash
while true; do
  shapepipe_run -c cfis/config_Gie_vos.ini
  ls -l data_exp/ | wc
  mv -i output/run_sp_Gie_*/get_images_runner/output/*.fits*fz data_exp
  ls -l data_exp/ | wc
  rm -rf output/run_sp_Gie_*
  update_runs_log_file.py
done
```

**Note:** Make sure that after all images are downloaded there is no `Gie` run in the output directory.
This would mess up later modules since `last:get_image_runner` could point to this run.

### Create tile links again (necessary?)

If necessary, e.g. because a previous `Git` run is no longer valid, re-create the symbolic links to the downloaded tiles with

```bash
job_sp_canfar.bash -p $psf `cat tile_numbers.txt` -j 1 -r symlink
```

### Uncompress tile weights

The downloaded tile weights are compressed. The following call uncompresses all.

```bash
shapepipe_run -c cfis/config_tile_Uz.ini
```

### Masks

There is no masking step. `ShapePipe` generates no masks: the sky-fixed
healsparse maps are queried once per object, by `mask_query` on the exposure
catalogues (`FLAG_EXT`, cut by `setools`) and by `make_cat` on the final tile
catalogue (`MASK_<band>` columns). Point the `MASK_PATHS` / `MASK_EXT_PATHS`
config entries at the maps and nothing else is needed — no star-catalogue
download, no rasterization, no `combine_runs.bash -c flag_*`. The only mask that
touches pixels is the instrument flag image shipped with each exposure, which
`split_exp` splits per CCD.

## Tile detection

We can finally run our first module using the canfar submission system. First, determine the number of optimal jobs such that at a given
time the allowed maximum of 512 running jobs is not exceeded.
Set as `N_PAR` (number of parallel jobs) a number between 1 and 8.

```bash
canfar_submit_job -j 16 -f tile_numbers.txt -P N_PAR -v -s
```

Now, run the previous command with that number `JMAX`

```bash
canfar_submit_job -j 16 -f tile_numbers.txt -P N_PAR -v -J JMAX
```

### Exposure Processing

#### Option 0: Global split (deprecated; used for earlier v1.x patch runs)

For this option, set `sp_local=0`.

**TODO: ** Split Uz and SpMh

For `sp_local=-` both `mh_local` (0, 1) are ok:

```bash
export mh_local=0
```

### Option 1: Local split exposures (recommended)

Optional: Enable flags for local split processing and merge header runs as

```bash
export sp_local=1
export mh_local=1
```

These flags are automatically set to 1 in the new job scripts.


Get single-HDU single-exposure IDs file (from missing 32 job):

```bash
summary_run P$patch 32
cp summary/missing_job_32_all.txt exp_shdu.txt
```

### Split exposures

First, determine the number of maximum jobs with the option `-s` (see above). Then, submit with

```bash
canfar_submit_job -j 2 -v -f exp_shdu.txt -v -P N_PAR -J JMAX
```

### Exposure detection

```bash
canfar_submit_job -j 32 -f exp_shdu.txt -v -P N_PAR -J JMAX
```

### Tile preparation

```bash
canfar_submit_job -j 64 -f tile_numbers.txt
```

### Tile shape measurement

```bash
canfar_submit_job -j 128 -f tile_numbers.txt
```

### Merge sub-catalogues

```bash
canfar_submit_job -j 256 -f tile_numbers.txt
```

### Create final catalogues

```bash
canfar_submit_job -j 512 -f tile_numbers.txt
```

This was the last `ShapePipe` module to run for main processing.


### Merge all final catalogues

The last step of `ShapePipe` processing is, per patch, to merget all final catalogues. This is done via a python script, as follows.
First, change to parent directory `/path/to/version` and run the following command for all patches

```bash
patchnum=`tr $patch P ''`
create_final_cat.py -m final_cat_$patch.hdf5 -i . -p $patch/cfis/final_cat.param \
    -P $patchnum -o $patch/n_tiles_final.txt -v
```

## Additional `ShapePipe` processing


### Create star Catalogue

We can additionaly create a combined star catalogue, with star shapes projecte from detector to world coordinates.
This is useful for validation and galaxy-PSF/star correlation diagnostics.

#### Combine all PSF runs

In each patch directory /path/to/version/$patch, run

```bash
combine_runs.bash -p $psf -c psf
```

to create a single output directory of PSF files (symbolic links).

Optionally, to create and plot results for this patch only:

```bash
shapepipe_run -c $SP_CONFIG/config_Ms_$psf.ini
shapepipe_run -c $SP_CONFIG/config_Pl_$psf.ini
```

#### Collate star catalogues

```{warning}
**This step no longer ships.** `collate_star_cat.py` gathered each validation
PSF run's star positions (X/Y/RA/DEC) and MCCD CCD id into the
`validation_psf_conv-*.fits` files the merge below consumes, and it was deleted
with the star-catalogue tooling the removed mask module needed (PR #847). The
merge step that follows therefore has no producer for its input in this repo
until the collation is reimplemented, and the `star_cat/` + `combine_runs.bash
-c psf_conv` staging around it does not apply either.
```

HSM shapes are **not** rotated into world coordinates at this step, and were not
before it was removed. The PSF and star ellipticities and sizes are measured
directly in sky coordinates during PSF interpolation (galsim
`FindAdaptiveMom(use_sky_coords=True)`) — see `tests/module/test_hsm_sky_coords.py`,
which pins that convention. The collation was a pass-through for shapes, which is
both why it carried no science and why it is straightforward to rebuild.

Merge all converted star catalogues and create `final-starcat.fits` (this
reads the `validation_psf_conv` files the collation above used to produce):

```bash
export SP_RUN=`pwd`
shapepipe_run -c ~/shapepipe/example/cfis/config_Ms_psfex_conv.ini
```

Rename to general PSF and star catalogue used for all ("a") sub-versions:


```bash
cp output/run_sp_Ms/merge_starcat_runner/output/full_starcat-0000000.fits \
    unions_shapepipe_psf_2024_v1.6.a.fits 
```

The FITS file `CATTYPE` (newer version) should be `validation_psf_conf`.

## Post-processing

The following post-processing steps are performed with the library `sp_validation`.

### Extract Information

First, we extract all information from the final catalogue, per patch. We copy
the parameter file and set links to the catalogues and `ShapePipe` config directory.

```bash
cd /path/to/version/$patch
cp ~/astro/repositories/github/sp_validation/notebooks/params.py .
ln -s /path/to/final_cat_$patchnum.hdf5  # not relative path ../final_cat_P$patchnum.hdf5 !
ln -s output/run_sp_MsPl/mccd_merge_starcat_runner/output/full_starcat-0000000.fits
ln -s ~/astro/repositories/github/shapepipe/example/cfis
```

Then edit `params.py`: Set patch name; set `wrap_ra` for P2.

Now we can run the script, recommended via job submission on candide. For large patches,
this requies a job with a large memory, e.g. with `mem=380000`
 

```bash
[squeue] python ~/astro/repositories/github/sp_validation/notebooks/extract_info.py
```

This creates a patch-wise comprehensive catalogue.

### Create global comprehensive catalogues

```bash
cd /patch/to/version
[squeue] python ~/astro/repositories/github/sp_validation/scripts/create_joint_comprehensive_cat.py \
    -v v1.6.c -v -p P1+P2+P3+P4+P5+P6+P7+P8+P9
```

This creates the file `unions_shapepipe_comprehensive_2024_v1.6.c.hdf5`.


### Apply structural masks

First, edit the Python script `~/astro/repositories/github/sp_validation/notebooks/demo_apply_hsp_masks.py`
to match catalogue name. Check the coverage mask input file (see below).
Run the script to apply the healsparse structural masks:

```bash
[squeue] python ~/astro/repositories/github/sp_validation/notebooks/demo_apply_hsp_masks.py
```

This creates the file `unions_shapepipe_comprehensive_struct_2024_v1.6.c.hdf5`.


### Define sample, calibrate catalogue

We are close to finally perform the last post-processing step, which is the calibration. First, the final galaxy sample
in question needs to be defined, with masks and cuts to apply from a `yaml` config file. A number of pre-defined files
can be found in `~/astro/repositories/github/sp_validation/calibration`.

For example, to create `v1.6.6`, the steps are:

```bash
cd /path/to/version
mkdir -p v1.6.6
cd v1.6.6
ln -s ~/astro/repositories/github/sp_validation/calibration/mask_v1.X.6.yaml config_mask.yaml
ln -s ..//unions_shapepipe_comprehensive_struct_2024_v1.6.c.hdf5 unions_shapepipe_comprehensive_struct_2024_v1.X.c.hdf5
[squeue] python ~/astro/repositories/github/sp_validation/calibrate_comprehensive_cat.py
```

calibrate_comprehensive


```{note}
The matched-star-catalogue and coverage-mask diagnostics in the next two
sections have largely moved out of ShapePipe into
[`sp_validation`](https://github.com/CosmoStat/sp_validation) / `cosmo_val`.
Several helpers referenced below — `merge_psf_cat.py`, `download_headers`,
`extract_field_corners`, `build_coverage_map`, `plot_coverage_map` — are no
longer shipped in this repository; see `sp_validation` for their current
equivalents. (`get_ccds_with_psf` and `summary_run` are still part of ShapePipe.)
```

### Create matched star catalogue

For diagnostics, a catalogue with multi-epoch shapes measured by ngmix matched with the validation star catalogue is used.
This is created as follows:

```bash
cd /path/to/version
merge_psf_cat.py [-V v1.6|-P P1+P2+...] -v
```

This creates the joint catalogue unions_shapepipe_star_2024_v1.6.a.fits .

### Create coverage mask

First, on canfar, move to the directory that has the patch subdirectories.

```bash
cd /path/to/version
```

#### Get exposure numbers

If the file `$patch/exp_numbers.txt` does not exist for a given patch, create it with the summary program

```bash
summary_run $patch 1
```

Now, create the list of CCDs that have PSF information with

```bash
get_ccds_with_psf -v -V v1.6
```

Next, download exposures headers

```bash
download_headers -i ccds_with_psfs_v1.6.txt -o headers_v1.6 -v
```

From the headers, the CCD corner coordinates are extracted with
```bash
extract_field_corners -i headers_v1.6 -v
```

Then, build the healsparse coverage mask file as
```bash
build_coverage_map
```

Use `plot_coverage_map` to create plots of the coverage mask.

## Extra Utilities

### Run in Terminal in Parallel

```bash
cat IDs.txt | xargs -I {} -P 16 bash -c 'init_run_exclusive_canfar.sh -j 512 -e {}'
```
