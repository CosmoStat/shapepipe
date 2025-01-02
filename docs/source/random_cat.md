# Create random catalogues and masks

This section describes how to create tile-based random catalogues and healpix
masks, and combined randoms and masks for a selection of tiles.

The masked regions are obtained on input from ShapePipe pixel mask ("pipeline flag")
files.

## Set up

### ID file and shell variables

First, if if does not exist already, create the file ``tile_numbers.txt`` containing a list of tile IDs,
one per line. This is the same format as the input file to ``get_images_runner``.
For example, link to a patch ID list,
```bash
ln -s tiles_PX.txt tile_numbers.txt
```
Next, set the run and config paths,
```bash
export SP_RUN=.
export SP_CONFIG=/path/to/config-files
```

### Get images or image headers

We need to footprint of the image tiles. If they have been downloaded for a ``ShapePipe`` run,
check that they are accessible as last run of the ``get_images_runner`` module.

If not, we can just download the headers to gain significant download time.
```bash
shapepipe_run -c $SP_CONFIG/config_get_tiles_vos_headers.ini
```

### Check pixel mask files

Make sure that all pixel mask files are present. If they have been downloaded from ``vos`` as ``.tgz`` files,
type
```bash
canfar_avail_results -i tile_numbers.txt --input_path . -v -m -o missing_mask.txt
```
In case of missing mask files, check whether they are present in the ``vos`` remote directory,
```bash
canfar_avail_results -i tile_numbers.txt --input_path vos:cfis/vos-path/to/results -v -m
```
If missing on ``vos``, process those tiles. If processing only up the the mask is necessary,
the following steps can be carried out,
```bash
job_sp -j 7 TILE_ID
job_sp -j 128 TILE_ID
```
The first command processes the tile up to the mask; the second line uploads the mask files
to ``vos``.

Now, download the missing masks with
```bash
canfar_download_results -i missing_mask.txt --input_vos vos-path/to/results -m -v
```
Untar .tgz files if required,
```bash
while read p; do tar xvf pipeline_flag_$p.tgz; done <missing_mask.txt
```

### Prepare combined mask input directory

To combine all mask files into one input directory, easy to find by the subsequent ``ShapePipe`` module
``random_runner``, type
```bash
prepare_tiles_for_final -c flag
```
To make sure everything went well, check the linked .fits files with
```bash
canfar_avail_results -i tile_numbers.txt --input_path output/run_sp_combined_flag/mask_runner/output -x fits -v -m
```

## Create random catalogue and helapix mask per tile

Run
```bash
shapepipe_run -c $SP_CONFIG/config_Rc.ini
```
The random catalogue and, with the config entry ``SAVE_MASK_AS_HEALPIX = True``
a healpix mask FITS file, will be written to disk.

## Create joint random catalogue

The individual tile-based random catalogues can be merged into a numpy
binary (``.npy``) file with
```bash
merge_final_cat -i output/run_sp_Rc/random_cat_runner/output -n random_cat -v
```

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

First, for convenience all image headers with WCS information are
linked from within one directory, with
```bash
prepare_tiles_for_final -i
```

Next, read all tile mask and WCS information, and create a joint full-sky
healpix mask with
```bash
/path/to/sp_validation/scripts/scripts/combine_hp_masks.py -p -v
```
With the option ``-p`` the mask is plotted in Mollweid projection.


