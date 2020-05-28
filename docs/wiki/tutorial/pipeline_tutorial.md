[Home](./shapepipe.md)

# ShapePipe tutorial

## Index

1. [Introduction](#introduction)
    1. [File types and names](#file-types-and-names)
    1. [CFIS processing](#cfis-processing)
    1. [Running the pipeline](#running-the-pipeline)
1. [Prepare input images](#prepare-input-images)
    1. [Select tiles](#select-tiles)
    1. [Download tiles and modify names](#download-tiles-and-modify-names)
    1. [Uncompress tile weights](#uncompress-tile-weights)
1. [Process single exposure images](#process-single-exposure-images)
    1. [Split images](#split-images)
    1. [Merge WCS headers](#merge-wcs-headers)
    1. [Mask images](#mask-images)
    1. [Extract sources](#extract-sources)
    1. [Select stars](#select-stars)
    1. [Model the PSF](#model-the-psf)
    1. [Validation tests](#validation tests)
1. [Process stacked images](#process-stacked-images)
    1. [Mask stacks](#mask-stacks)
    1. [Extract sources on stacks](#extract-sources-on-stacks)
    1. [Interpolate multi-epoch PSF](#interpolate-multi-epoch-psf)
    1. [Create weight postage stamps](#create-weight-postage-stamps)
    1. [Compute spread model](#compute-spread-model)
    1. [Prepare shape measurement](#Prepare-shape-measurement)
    1. [NGMIX : Shape measurement](#NGMIX-:-Shape measurement)
    1. [Make final catalog](#Make-final-catalog)

## Introduction

The `ShapePipe` pipeline can process single-exposures images, and stacked images. The input images have to be calibrated beforehand for astrometry and photometry.

***WARNING /!\ :*** All file paths for the following examples are relative. When running on a cluster, you need to make sure that these paths are accessible on all computing nodes.
Absolute paths are recommended to avoid problems.

### File types and names

The `ShapePipe` pipeline handles different image and file types, some of which
are created by the pipeline during the analysis. These file types are listed below. All
files follow a (configurable) naming and numbering convention, to facilitate bookkeeping for
tracking relevant image information. We adopt a numbering schemes as follows.

- Single-exposure mosaic image.  
  Multi-HDU FITS file containing a mosaic from multiple CCDs of a single exposure (an exposure is also called epoch).
  Each CCD is stored in a different HDU.
  These files are used on input by `ShapePipe`. The pixel data can contain the observed image, a weight map, or a flag map.
  These images are typically created by a telescope analysis software (e.g.~`pitcairn`).  
  Convention: None. The file names are in general determined by this software, e.g.~they contain the run ID. and do not
  need to be changed to be read by `ShapePipe`.  
  Examples from CFIS: `2228303p.fits`, `2214439p.flag.fits`.

- Single-exposure single-CCD image.  
  FITS file containing a single CCD from an individual exposure. The pixel data can contain the observed image, a weight map, or a flag map.  
  Default convention: **<image_type>-<exposure_name>-<CCD_number>.fits**  
  Examples: `image-2079614-9.fits`, `weight-2079614-3.fits`

- Stacked images  
  FITS file containing a stack by co-adding different single exposures, created by software such as `swarp`.
  A stacked image is also called *tile*.
  The pixel data can contain the observed image, a weight map, or a flag map.  
  Default convention: **<image_type>-<number>.fits**  
  Examples: `CFIS_image-51.fits`, `pipeline_flag-2134.fits`

### CFIS processing

`ShapePipe' splits the processing of CFIS images into three parts: 1.) Preparation of the input images; 2.) Processing of single exposure images;
3.) Processing of stacked images. The single exposures are first split into single-CCD images, which are processed in turn and
independently.

The preparation of input images is done before running the actual pipeline, using auxilliary scripts.

The processing of single exposure images contains the following steps:
  * Split exposure into single-exposure single-CCD images
  * Create masks for bright stars, spikes, borders, Messier objects, ...
  * Detect stars
  * Model the PSF
  * Validate the PSF model (optional)

The processing of stacked images has the following tasks:
  * Create mask for bright stars, spikes, border, Messier objects, ...
  * Detect all sources
  * Interpolate the PSF model at the location of each source for all contributing exposures
  * Create postage stamps necessary for the *spread model*, for galaxy selection
  * Compute the spread model for each source
  * Create postage stamps, for the shape measurement
  * Measure galaxy shapes
  * Merge all results into one parent catalog

The following flowchart visualised these processes:

![ShapePipe_FlowChart](./ShapePipe_v0.0.1.png)

In the following, the individual processing steps are described in detail.

### Running the pipeline

See the main `ShapePipe` readme for more details.

In the following, to have consistent paths, we assume the existence of the following directories or links:
- `CFIS`: Path of downloaded CFIS images.
- `$SP_CONFIG`: Path to configuration files. In our example this is `/path/to/shapepipe/example/GOLD`.

In `~/ShapePipeRun` the following subdirectories need to be created by the user:
- `input_tiles`: Symbolic links to CFIS input tile images and weight files
- `input_exposures`: Symbolic links to CFIS single-exposure images, weights, and flag file
- `output`: General output generated by the pipeline (log files,
  diagnostics, statistics, output images, catalogues)
- `output_headers`: Single-exposure headers with WCS information`
- `output_star_cat`: Star catalogues

In general, a call to the pipeline is done as follows:

```bash
shapepipe_run -c $SP_CONFIG/config_<module[s]>.ini
```

The config file `config_<module[s]>.ini` contains the configuration for one or more modules.


## Prepare input images

### Select tiles

The selection of images on input can be done in the config files of the relevant modules, by specifying input
path(s) and input file name patterns. Thus, a sub-selection of images in a given input directory can be made.
However, one might want to pre-select specific images before the pipeline is run, for example to find all available images (exposures and stacks) in a given sky area. The resulting files can then be copied to a new, dedicated directory (or alternatively linked using symbolic links), or downloaded to a local machine.

Images can also be selected to cover a given sky area, with the script `cfis_field_select`.

First, find the tile(s) covering a given coordinate or area. For example, the tile for a Planck cluster at R.A.=255.66 deg, dec= 34.05 deg can be found with the `--coord` option:
```bash
cfis_field_select -i ~/CFIS/tiles+weights_DR2.txt --coord 255.66deg_34.05deg -v -t tile
```
The input text file (with `-i`) contains a list of CFIS tiles.

We also need to get the weight files for the tile.
```bash
cfis_field_select.py -i ~/CFIS/tiles+weights_DR2.txt --coord 255.66deg_34.05deg -v -t weight
```

### Download tiles and modify names

The tile images and weights selected in the previous section need to be findable by `ShapePipe` in the tiles input directory `input_tiles`. Either download the images and weights to, or, if they are already stored locally on a hard disk, create symbolic links in `input_tiles`. Now is a good time to make a necessary small change to the file names. Any dot (`.`) that does not indicate a file extension needs to be replaced. In addition, file type specifiers need to appear before the tile number. Therefore, images and weights need to be renamed, for example according to the following scheme:
```bash
mv CFIS.424.248.r.fits CFIS_image_424_248.r.fits
mv CFIS.424.248.r.weight.fits.fz CFIS_weight_424_248.r.weight.fits.fz
```
Of course the above renaming can be done at the same time as copying/creating links.

### Uncompress tile weights

The weights of the stacks are compressed .fits.fz files, they need to be uncompressed before the pipeline is run.
This can be done with the script `cfis_weight_format`. For example:

```bash
for fzweight in input_tiles/CFIS_weight-*.fits.fz; do
  weight=`basename $fzweight .fz`
  cfis_weight_format -i $fzweight -o input_tiles/$weight
done
```


### Find exposures

Once the resulting tiles and weight images are downloaded, we need to get the exposure images that where co-added to produce the tiles.
These can be found from the tile header, with the `--tile` option. We need all three single-exposure types, data, weights, and flags:
```bash
~/ShapePipe/scripts/python/cfis_field_select.py -i ~/CFIS --tile -v -t exposure
~/ShapePipe/scripts/python/cfis_field_select.py -i ~/CFIS --tile -v -t exposure_weight.fz
~/ShapePipe/scripts/python/cfis_field_select.py -i ~/CFIS --tile -v -t exposure_flag.fz
```

The resulting files need to be downloaded.



## Process single exposure images

### Split images

**Module:** split_exp_runner   
**Parent:**  none  
**Input:** single-exposure images, weights, flags  
**Output:** single_exposure single-CCD files for input images, weights, flags

The first step of single-exposure processing is to split the single-exposures images into
files that contain one CCD each.

The example config file is `~/ShapePipe/example/GOLD/config_split_exp.ini`.
On input, we need to specify the three input types (exposures, weights, flags),
and their extensions. This happens in the `[FILE]` section:
```ini
[FILE]
FILE_PATTERN = image,weight,flag
FILE_EXT = .fitsfz,.fitsfz,.fitsfz
```
On output, the same three file types are required. The number of MegaCAM CCDs is 40:
```ini
[SPLIT_EXP_RUNNER]
OUTPUT_SUFFIX = image,weight,flag
N_HDU = 40
```

Run
```bash
mkdir -p output
~/ShapePipe/shapepipe_run.py -c ~/ShapePipe/example/GOLD/config_split_exp.ini
```

On success, files accordingt to the three output types are created.

### Merge WCS headers

**Module:** merge_headers  
**Parent:** split_exp_runner  
**Input:** single-exposure single_CCD images, weights, flags  
**Output:** Single SQL file with combined header information  

This pipeline module saves the WCS information (image
transformation and distortions, computed during astrometrical calibration)
for each CCD. At the end, this information has to be merged back into a single file.  
Specify the output path:
```ini
[MERGE_HEADER_RUNNER]
OUTPUT_PATH = $HOME/ShapePipeRun/output_headers
```
Create the output directory, and run the pipeline:

```bash
mkdir -p output_headers
~/ShapePipe/shapepipe_run.py -c ~/ShapePipe/example/GOLD/config_merge_headers.ini
```
Since this produces a single output file
instead of a file per input image, it is more convenient to have this file in
a separated directory for later use.

On success, a single `.sqlite` file is created.

> Note: On cc this module failed to run alone. It should be run together with
the previous one, `split_exp_runner`.

### Mask images

**Module:** mask_runner   
**Parent:** split_exp_runner  
**Input:** single-exposure single-CCD images, weights, flags [, star catalogue]  
**Output**: single-exposure single-CCD flag files

This module creates masks for bright stars, diffraction spikes, Messier objects,
borders, and other artifacts. It joins the newly created mask with the already
existing masks (from the input flag files) of cosmic rays and various artifacts.  

Those mask parameters are read from a second config file, whose path
needs to be specified:
```ini
[MASK_RUNNER]
MASK_CONFIG_PATH = $HOME/ShapePipe/example/GOLD/config.mask
```
In this mask config file the default parameters can be kept in the most part.
These parameters specify the mask properties for border, halos, spikes, Messier
objects, and external flag input (in our case provided from CFIS pre-processing).

It points to various default parameter files for the different mask types,
make sure that that paths are correct, in our case
`$HOME/ShapePipe/example/GOLD/mask_default/` in front of each file name.

To distinguish the newly created output flag files from the input ones,
a suffix is added:
```ini
SUFFIX = pipeline
```
Next, this module requires a star catalogue containing position and magnitude
of bright stars. By default this is automatically created by accessing online
star catalogues. Since in some cases computing nodes on clusters might not have
internet access, such a catalog can also be created for each image, before running
this module as follows:
```bash
mkdir -o output_star_cat
~/ShapePipe/scripts/python/create_star_cat.py input_exposures output_star_cat exp
```
Then, the star catalogue needs to be specified as input in the config file,
and a flag has to be set::
```ini
[FILE]
INPUT_DIR = last:split_exp_runner,${HOME}/ShapePipeRun/output_star_cat
[MASK_RUNNER]
USE_EXT_STAR = True
```

If instad the star catalogues can be accessed during the pipeline running,
the config files looks as follows:
```ini
[FILE]
INPUT_DIR = last:split_exp
[MASK_RUNNER]
USE_EXT_STAR = False
```

Finally, run the module:
```bash
~/ShapePipe/shapepipe_run.py -c ~/ShapePipe/example/GOLD/config_mask.ini
```

On success, pipeline-flag files are created.

**Diagnostics:** Open a single-exposure single-CCD image and the corresponding pipeline flag
in `ds9`, and display both frames next to each other. Example
```bash
ds9 output/shapepipe_run_2020-03-03_15-31-00/split_exp_runner/output/image-2113737-10.fits output/shapepipe_run_2020-03-03_17-29-34/mask_runner/output/pipeline_flag-2113737-10.fits
```
Choose `zoom fit` for both frames, click `scale zscale` for the image, and `color aips0` for the flag, to display something like this:

<img width="250" src="diag_mask.png">

By eye the correspondence between the different flag types and the image can be
seen. Note that the two frames might not match perfectly, since (a) WCS
information is not available in the flag file FITS headers; (b) the image can
have a zero-padded pixel border, which is not accounted for by `ds9`.

### Extract sources

**Module:** sextractor_runner  
**Parent:** split_exp_runner, mask_runner  
**Input:** single-exp_single-CCD image, weights, flags  
**Output:** sextractor catalogue

The purpose of source extraction/source identification on single exposures is
to select stars in the next step. Therefore, a relatively high
detection threshold is chosen to avoid to detect too many low-SNR
artifacts, and to reduce the output catalogue size. The following config entry
is
```ini
DETECT_THRESH    2.             # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
```
in the file `$HOME/ShapePipe/example/GOLD/sextractor_default/default.sex`.

On success, SEXtractor catalogue FITS files are produced.

### Select stars

**Module:** setools  
**Parent:** sextractor_runner  
**Input:** sextractor catalog  
**Output:** masked sextractor catalogue

For the star selection we use the size-magnitude plane. We first find the
stellar locus, by computing the FWHM mode of all objects. Objects are selected
within some range in FWHM around the mode, and within a magnitude range.

The selection criteria are given in a selection configuration file, whose name is specified
in the `setools` section:
```ini
[SETOOLS_RUNNER]
SETOOLS_CONFIG_PATH = $HOME/ShapePipe/example/GOLD/star_selection.setools
```
The selection config file `star_selection.setools` first defined a pre-selectione (or filter, or mask),
such that the subsequent computation of the mode is more stable:
```ini
[MASK:preselect]
MAG_AUTO > 0
MAG_AUTO < 20
FWHM_IMAGE > 0.3 / 0.186
FWHM_IMAGE < 1.5 / 0.186
FLAGS == 0
IMAFLAGS_ISO == 0
NO_SAVE
```
> Note the normalisation by the pixel scale to express **FWHM_IMAGE** in arc seconds.

The star sample is then selected as follows:
```ini
[MASK:star_selection]
MAG_AUTO > 17.
MAG_AUTO < 21.5
# NOTE : unit is pixel
FWHM_IMAGE <= mode(FWHM_IMAGE{preselect}) + 0.2
FWHM_IMAGE >= mode(FWHM_IMAGE{preselect}) - 0.2
FLAGS == 0
IMAFLAGS_ISO == 0
```
In addition, the selection config file can contain instructions to create plots and
statistics of the selected population(s). The former can be scatter plots and histograms,
the former can include mean, mode, extrema, and standard deviation
of any quantity from the SExtractor input catalogue, the number of selected objects, etc..

On success, masked SEXtractor catalogues are created in `mask`, plots are put in `plots`,
and text files with the computed statistics in `stats`.

The following plots show an example, CCD #10 of exposure 2113737.


| Size-magnitude plot | Star magnitude histogram | Stars in CCD (mag) | Stars in CCD (size) |
| --- | --- | --- | --- |
| <img width="250" src="size_mag-2113737-10.png" title="Size-magnitude plot with star selection"> | <img width="250" src="hist_mag_stars-2113737-10.png" title="Magnitude histogram of selected stars"> | <img width="250" src="mag_star_field-2113737-10.png" title="Magnitude distribution in CCD"> | <img width="250" src="fwhm_field-2113737-10.png" title="Size distribution in CCD"> |
| The stellar locus is well-defined | Magnitude distribution looks reasonable | Stars are relatively homogeneously distributed over the CCD | The uniform and small seeing of CFHT is evident |

To contrast the last plot, here is the case of CCD #35 (lower right corner), which shows a known (but yet unexplained) lack of stars
in the lower parts:

<img width="250" src="fwhm_field-2113737-35.png" title="Size distribution in CCD">

The statistics output file, also for CCD #10 is:
```bash
(shapepipe)  dap ~/ShapePipeRun $ cat output/shapepipe_run_2020-03-05_10-00-26/setools_runner/output/stat/star_stat-2113737-10.txt
# Statistics
Nb objects full cat = 1267
Nb stars = 160
stars/deg^2 = 6345.70450519073
Mean star fwhm selected (arcsec) = 0.7441293125152588
Standard deviation fwhm star selected (arcsec) = 0.014217643037438393
Mode fwhm used (arcsec) = 0.7345179691314697
Min fwhm cut (arcesec) = 0.7159179691314698
Max fwhm cut (arcsec) = 0.7531179691314697
```

### Model the PSF

**Module:** psfex  
**Parent:** setools_runner  
**Input:** setools_star_selection  
**Output:** star catalogue, psf file

The PSF modeling is done with `PSFEx`. The psfex module configuration section
has to point to the executable and the default psfex config file, which does not
need to be changed.
```ini
[PSFEX_RUNNER]
EXEC_PATH = psfex
DOT_PSFEX_FILE = ./example/test_psfex/default.psfex
```

On success, FITS files containing the star catalalogue (`psfex_cat-*.cat`) and the PSF at
the stars' positions (`star_selection-2159358-9.psf`) are created.

### Validation tests

#### Interpolate PSF to star positions

**Module:** psfinterp_runner  
**Parent**: psfex, setools  
**Input:** star catalogue, psfex_catalog  
**Output:** star catalogue

The interpolation of the PSF on single exposures alone is not required for shape
measurement. But to carry out validation tests on the model we need the
know the PSF at the position of the stars used. For that we run the module
`psfinterp_runner` in `VALIDATION` mode. The following options are required:
```ini
# Define the way psfexinter will run.
# CLASSIC for classical run.
# MULTI-EPOCH for multi epoch.
# VALIDATION for output allowing validation (only on single epoch !)
MODE = VALIDATION
# Column names of position parameters
POSITION_PARAMS = XWIN_IMAGE,YWIN_IMAGE
# If True, measure and store ellipticity of the PSF
GET_SHAPES = True
# Number minimal of stars require to interpolate the PSF on the CCD
STAR_THRESH = 20
# Maximum value allowed for the global chi2 of the model for the CCD
CHI2_THRESH = 2
```

On success, validation PSF catalogues are created.

#### Merge PSF catalogues

**Module:** merge_star_cat_runner  
**Parent**: psfexinterp  
**Input:** psfex cataloguefex  
**Output:** merged SExtractor catalogue

The star and PSF catalogues from the previous module are now merged
into a single FITS file, encoding the CCD number, with output directory given
in the module section:
```ini
OUTPUT_PATH = $SP_RUN/psf_validation
```
The (single) output file is then `SP_RUN/psf_validation/full_starcat.fits`.


#### Plot focal-plane ellipticity, size and residuals

Outside the pipeline, create plots of PSF, model, and residual
ellipticity and shape:
```bash
~/ShapePipe/scripts/python/MeanShapes.py -o psf_validation -x 10 -y 20 -i psf_validation/full_starcat.fits -v
```


## Process stacked images

### Mask stacks

**Module:** mask_runner   
**Parent:** none  
**Input:** stack image, stack weight [, star catalogue)   
**Output:** stack flag

This is the analogue of the single-exposure mask module(#mask-images), but for stacks.
The mask configuration file needs to be the tile-specific one:
```ini
[MASK_RUNNER]
MASK_CONFIG_PATH = $HOME/ShapePipe/example/GOLD/config.mask
```
This configuration file has a few differences compared to the single-exposure one.
First, the border region can be much smaller. Second, no external flag files exist,
and third the temporary directory should be different:
```ini
[BORDER_PARAMETERS]
BORDER_WIDTH = 5
[EXTERNAL_FLAG]
F_MAKE = False
[OTHER]
TEMP_DIRECTORY = .temp_tiles
```

Run the package:
```bash
~/ShapePipe/shapepipe_run.py -c ~/ShapePipe/example/GOLD/config_mask_tile.ini
```

> Note: On the cc the star catalogues created by `create_star_cat.py` did not
work with the pipeline. Thus for the moment, the mask package needs to be run
on the login node.

### Extract sources on stacks

**Module:** sextractor_runner (in multi-epoch mode)  
**Parent:** mask_runner  
**Input:** tile_image, tile_weight, tile_flag  
**Output:** SExtractor catalogue with multi-epoch information

To detect a maximum of sources on the tiles, we set a low detection threshold.
In addition, a a post-processing step is run to find all epochs that contributed
to a given detected object. The different entries compared to the
single-exposure case (#extract-sources) are thus:
parameter file.
```ini
[SEXTRACTOR_RUNNER]
# Point to the tile-specific SExtractor parameter file
DOT_SEX_FILE = $HOME/ShapePipe/example/GOLD/default_tile.sex

# Necessary for tiles, to enable multi-exposure processing
MAKE_POST_PROCESS = TRUE

# Multi-epoch mode: Path to file with single-exposure WCS header information
LOG_WCS = ${HOME}/ShapePipeRun/output_headers/log_exp_headers.sqlite

# World coordinate keywords, SExtractor output. Format: KEY_X,KEY_Y
WORLD_POSITION = XWIN_WORLD,YWIN_WORLD

# Number of pixels in x,y of a CCD. Format: Nx,Ny
CCD_SIZE = 33,2080,1,4612

```
The different entries in the SExtractor parameter file are
```ini
DETECT_THRESH    1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
DEBLEND_MINCONT  0.0005         # Minimum contrast parameter for deblending
BACK_TYPE        MANUAL         # AUTO or MANUAL
```
On success, output FITS SEXtractor files are created. They have the following format:
```python

HDU  Name        Ver Type         Cards   Dimensions    Format
  0  PRIMARY       1 PrimaryHDU       4   ()      
  1  LDAC_IMHEAD   1 BinTableHDU     12   1R x 1C       [8560A]   
  2  LDAC_OBJECTS  1 BinTableHDU    180   25133R x 45C  [1J, 1I, 1E, 1E, ...]   
  3  EPOCH_0       1 BinTableHDU     16   25133R x 3C   [1J, 7A, 1J]   
  4  EPOCH_1       1 BinTableHDU     16   25133R x 3C   [1J, 7A, 1J]   
 ...
```

The first 3 HDUs correspond to the usual SExtractor output. The following HDUs
`EPOCH_X` contain the single-exposure information contributing to the objects
on the stack, one HDU for each epoch. For those HDUs the columns are:
* **NUMBER** : object id atributed by SExtractor
* **EXP_NAME** : name of the single exposure (same for all objects of one HDU)
* **CCD_N** : CCD number in which the object is following the MegaCam numbering (set to -1 if the object is not on the single exposure)

Those additionnal HDUs contain all the multi-epoch informations we need for the
following processing steps.

### Interpolate multi-epoch PSF

**Module:** psfexinterp_runner   
**Parent:** SExtractor (in multi-epoch postprocessing mode)  
**Input:** SExtractor catalog with multi-epoch information   
**Output:** tile PSF dictionary

This step interpolates the PSF to the position of all detected sources on the
tiles for all epochs where the object appears.
The run mode has to be set to multi-epoch. In addition, the single-exposure
PSF information needs to be read. Unfortunately, at present, this cannot be
provided automatically pointing to the corresponding run. The simplest way
is to set a symbolic link to the corresponding run output directory:
```bash
ln -s output/shapepipe_run_2020-03-19_18-28-03/psfex_runner/output input_psf
```
and to indicate the link name in the config file:
```ini
[PSFEXINTERP_RUNNER]
MODE = MULTI-EPOCH
ME_DOT_PSF_DIR = input_psfex
```

On success, `sqlite` output catalogues are created containing the PSF on
vignets (postage stamps). The structure is similar
to a dictionary with the following format:
```python
{'object_id': {'exp_name-CCD_number' : {'VIGNET': numpy.ndarray, 'SHAPES': {}}}
```
For example:
```python
{'1': {'2104127-35': {'VIGNET': numpy.ndarray, 'SHAPES': {}},
       '2105224-13': {'VIGNET': numpy.ndarray, 'SHAPES': {}},
       '2105382-3':  {'VIGNET': numpy.ndarray, 'SHAPES': {}}},
 '2': {'2104127-33': {'VIGNET': numpy.ndarray, 'SHAPES': {}},
       '2105224-12': {'VIGNET': numpy.ndarray, 'SHAPES': {}},
       '2105382-2':  {'VIGNET': numpy.ndarray, 'SHAPES': {}}},
 ...}
```


### Create weight postage stamps

**Module:** vignetmaker_runner   
**Parent:** sextractor_runner (in multi-epoch mode)  
**Input:** SExtractor catalog with multi-epoch information, tile weight  
**Output:**  vignet FITS table

This modules is a pre-processing step to compute the spread model. This galaxy
classification computation, performed in the nex step, requires
* The vignets of the objects on the stack
* The PSF information at the objects' location
* The vignets of the weight image at the objects' location

The first two have already been obtained, thus we only need to extract the
weights. These are needed to weigh the object images for corresponding
comparison to the weighted single-exposure PSFs for the spread model
classification.

On success, FITS tables with vignets containing the weight for each object.

### Compute spread model

**Module:** spread_model_runner  
**Parent:** psfex_runner (single-exposure), psfexinterp_runner (tile), vignetmaker_runner  
**Input:** psfex catalogue, tile psf dictionary, weight vignet  
**Output:** SExtractor catalogue

The spread model for each object is computed, which serves to classify a
sub-set of detected objects on the tiles as galaxies.

> Note: The effective PSF on the tiles is approximated by a weighted sum
of the single-exposure PSFs. The weight for each epoch is the average over
the postage stamp.

On success, SExtractor catalogues with galaxy number, magnitude, spread model,
and an error estimate is produced.

### Create single-exposure postage stamps

**Module:** vignetmaker_runner2  
**Parent:**  sextractor_runner  
**Input:** sextractor_catalog   
**Output:** single-exposure vignet dictionary

This second iteration of the vignet creation module is the last step in
preparation for galaxy shape measurement. Multi-epoch shape measurement equires
* The SExtractor tile catalog with spread-model information
* The vignets of single-exposure images for tile-detected objects
* The vignets of single-exposure weights at the position of tile-detected objects
* The vignets of single-exposure flags at the position of tile-detected objects
* The vignets of the single-exposure background images at the position of tile-
  detected objects
* The vignets of the single-exposure PSFs

The first and last items of the list were obtained in (#compute-spread-model)
and (#interpolate-multi-epoch-psf), respectively. The missing middle three are
thus to be extracted here. For technical reasons, we have to use for the
moment the module `vignetmaker_runner2`, run in `MULTI-EPOCH` mode. To work
in tile coordinates, we need spherical WORLD coordinates instead of Cartesian IMAGE (pixel)
units:
```ini
[VIGNETMAKER_RUNNER2]

MODE = MULTI-EPOCH

# Coordinate frame type, one in PIX (pixel frame), SPHE (spherical coordinates)
COORD = SPHE

# Coordinate frame type, one in PIX (pixel frame), SPHE (spherical coordinates)
COORD = SPHE
POSITION_PARAMS = XWIN_WORLD,YWIN_WORLD

# Additional parameters for path and file pattern corresponding to single-exposure
# run outputs
ME_IMAGE_DIR = input_split_exp,input_split_exp,input_split_exp,input_sextractor
ME_IMAGE_PATTERN = flag,image,weight,sexcat_background
```
The last entries indicate four input paths and corresponding file patterns, for:
single_exposure single-CCD images, weights, flags, created in (#split-images), and
for single-exposure background vignet file, created in (#extract-sources).

On success, `sqlite` dictionaries with vignets for the image, weight, flag, and
background are created.

### Multi-epoch shape measurement with `ngmix`

**Module:** ngmix_runner   
**Parent:**
**Input:** sextractor_catalog, single_exp_image_vignet, single_exp_bkg_vignet,
tile_psf, single_exp_weight_vignet, single_exp_flag_vignet  
**Output:** SExtractor catalogue

Now we run the shape measurement. At the moment it's done with NGMIX. Most of the features are hard coded (will be more flexible in the future). Here is a commented example config file for the pipeline :

```ini
[NGMIX_RUNNER]

# Create with the split_exp_hdu module
LOG_WCS = /path/to/file/containing/WCS/information/log_exp_headers.npy
```

### Make final catalog

**Module :** make_catalog_runner   
**Module inputs :** sextractor_catalog, spread_model, tile_psf, ngmix_catalog

We finally merge all the results into one catalog per tiles. Here is a commented example config file for the pipeline :

```ini
[MAKE_CATALOG_RUNNER]

# If true will add a column in the final catalog : "SPREAD_CLASS"
SM_DO_CLASSIFICATION = True
# The classification is done by computing :
# classif = SPREAD_MODEL + 5/3 * SPREADERR_MODEL
# The cut for the star : |classif| < SM_STAR_STRESH
SM_STAR_STRESH = 0.003
# The cut for the galaxies : classif > SM_GAL_THRESH
SM_GAL_THRESH = 0.01
```

**Tips :** At the end we have one catalog per tiles. Here are a piece of code to merge all of them (assuming you are in the directory with all the catalogs) :

```python
import os
import numpy as np
from shapepipe.pipeline import file_io as io

all_file = os.listdir()

cat = io.FITSCatalog(all_file[0])
cat.open()

data = np.copy(cat.get_data())

cat.close

for f in all_file[1:]:
  cat = io.FITSCatalog(f)
  cat.open()

  data = np.concatenate((data, np.copy(cat.get_data())))

  cat.close()

final_catalog = io.FITSCatalog('final_cat.fits', open_mode=io.BaseCatalog.OpenMode.ReadWrite)

final_catalog.save_as_fits(data, ext_name='RESULTS')

```
