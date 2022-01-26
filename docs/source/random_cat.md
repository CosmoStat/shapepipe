# Create random catalogue

This section describes how to create a random catalogue from ShapePipe
mask files corresponding to a selection of tiles.

## Set up

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

## Get images or image headers

We need to footprint of the image tiles. If they have been downloaded for a ``ShapePipe`` run,
check that they are accessible as last run of the ``get_images_runner`` module.

If not, we can just download the headers to gain significant download time.
```bash
shapepipe_run -c $SP_CONFIG/config_get_tiles_vos_headers.ini
```

## Check mask files

Make sure that all mask files are present. If they have been downloaded from ``vos`` as ``.tgz`` files,
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
canfar_download_results.bash -i missing_mask.txt --input_vos vos-path/to/results -m -v
```
Untar .tgz files if required,
```bash
while read p; do tar xvf pipeline_flag_$p.tgz; done <missing_mask.txt
```

## Prepare combined mask input directory

To combine all mask files into one input directory, easy to find by the subsequent ``ShapePipe`` module
``random_runner``, type
```bash
prepare_tiles_for_final -c flag
```
To make sure everything went well, check the linked .fits files with
```bash
canfar_avail_results -i tile_numbers.txt --input_path output/run_sp_combined_flag/mask_runner/output -x fits -v -m
```

## Create random catalogue

First, create a random catalogue for each input tile and mask file,
```bash
shapepipe_run -c $SP_CONFIG/config_Rc.ini
```
Next, merge those catalogues into a numpy binary (``.npy``) file,
```bash
merge_final_cat -i output/run_sp_Rc/random_cat_runner/output -n random_cat -v
```

## Results

We can plot the random objects,
```bash
python ~/astro/repositories/github/sp_validation/scripts/plot_rand.py
```
and finally compute the effective survey area,
```bash
~/astro/repositories/github/sp_validation/scripts/compute_area.py
```

