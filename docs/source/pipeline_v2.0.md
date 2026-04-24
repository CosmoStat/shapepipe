# Running `ShapePipe` processing and post-processing pipelines on CANFAR

Documentation to create ShapePipe output products for catalogues v2.x.

## Initialise directory structure

```bash
init_run_v2.0.sh
```
sets up the directory structure. This will be

v2.0/
├── tiles/
│   ├── 301/
│   │   ├── 301.278/
│   │   ├── 301.279/
│   │   └── ...
├── exp/
│   ├── 21/
│   │   ├── 21163916
│   │   └── ...
├── cfis -> /arc/home/kilbinger/shapepipe/example/cfis
├── tile_numbers -> /arc/home/kilbinger/shapepipe/auxdir/CFIS/tiles_202604/tiles_r.txt
└── debug/


### Interactive job from the terminal for a single tile

Run bit-coded jobs
```bash
run_job_canfar_v2.0.sh -e ID -j <job>
```
with job processing tiles:
- 1: download tiles
- 2: uncompress tile weights
- 4: find exposures
then exposures:
- 8: download exposures
- 16: split exposures into single-CCD HDUs
- 32: mask exposures
- 64: process stars (selection, PSF movel)
then back to tiles:
- 128: select objects (using external catalogue)
- 256: create object postage stamps

## CANFAR Setup

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


## From previous setup

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

#### Convert star catalogue to wCS

Convert all input validation PSF files and create directories per patch `P?`.
Create files `validation_psf_conv-<patchnum>-<idx>.fits` (for the v1.4 setup only one file):

```bash
cd /path/to/version
mkdir stat_car
cd star_cat
```

For each patch run

```bash
convert_psf_pix2world.py -i .. -P $patchnum -v
```

Combine previously created files as links within one ShapePipe run directory (for the v1.4 setup only one link).
First (and optiohnal), create a subdir for a run and link to the input patches:

```bash
cd /path/to/version/star_cat
mkdir v1.6
ln -s ../P1
ln -s ../P2
...
```

Next, create links to all `validation_conv` runs:

```bash
combine_runs.bash -p psfex -c psf_conv
```

Merge all converted star catalogues and create `final-starcat.fits`:

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
