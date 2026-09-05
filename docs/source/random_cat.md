# Create random catalogues and masks

This section describes how to create tile-based random catalogues and healpix
masks, and combined randoms and masks for a selection of tiles.

The masked regions are obtained on input from ShapePipe pixel mask ("pipeline
flag") files.

```{note}
The `random_cat` module itself is current. The staging and joint-mask steps
below predate the Snakemake workflow and are retained for reference; the
helpers they used to call (`prepare_tiles_for_final`, `merge_final_cat`,
`canfar_avail_results`, `canfar_download_results`) are no longer shipped, and
where a step named one it now says what the step has to achieve instead.
```

## Set up

### ID file and shell variables

First, if it does not exist already, create the file ``tile_numbers.txt``: a
tile list, one ID per line. This is the same format as the input file to
``get_images_runner``.

Next, set the run and config paths,
```bash
export SP_RUN=.
export SP_CONFIG=/path/to/config-files
```

### Get images or image headers

We need the footprint of the image tiles. If they have been downloaded for a
``ShapePipe`` run, check that they are accessible as last run of the
``get_images_runner`` module.

If not, we can just download the headers to gain significant download time.
```bash
shapepipe_run -c $SP_CONFIG/config_get_tiles_vos_headers.ini
```

### Collect the pixel mask files

``random_cat_runner`` reads all pixel mask files from a single input
directory. Make sure every tile in ``tile_numbers.txt`` has its
``pipeline_flag`` file there — by pointing the config at the ``mask_runner``
output of the runs being combined, or by linking the files into one directory.

## Create random catalogue and healpix mask per tile

Run
```bash
shapepipe_run -c $SP_CONFIG/config_Rc.ini
```
The random catalogue and, with the config entry ``SAVE_MASK_AS_HEALPIX = True``
a healpix mask FITS file, will be written to disk.

## Create joint random catalogue

The individual tile-based random catalogues are per-tile FITS files under
``output/run_sp_Rc/random_cat_runner/output``. Merging them into a single
numpy binary was done by ``merge_final_cat``, retired at `2ef07e45`; nothing
in the workflow supersedes it, so the concatenation is currently a hand step.

### Results

We can plot the random objects,
```bash
python ~/astro/repositories/github/sp_validation/scripts/plot_rand.py
```
and also compute the effective survey area,
```bash
~/astro/repositories/github/sp_validation/scripts/compute_area.py
```

## Create joint healpix mask

Read all tile mask and WCS information, and create a joint full-sky healpix
mask with
```bash
/path/to/sp_validation/scripts/scripts/combine_hp_masks.py -p -v
```
With the option ``-p`` the mask is plotted in Mollweid projection.
