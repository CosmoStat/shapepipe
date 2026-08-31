# Create random catalogues and masks

This section describes how to create tile-based random catalogues and healpix
masks, and combined randoms and masks for a selection of tiles.

The masked regions are obtained on input from per-tile pixel mask images.

```{warning}
**ShapePipe no longer produces those images.** The pipeline generates no masks
at all: sky-fixed masks are healsparse maps, queried once per object into
catalogue columns (see [Masks](pipeline_tutorial.md#masks)), and tiles have no
flag image. `random_cat_runner` therefore needs its mask images supplied from
outside the pipeline — point its second `INPUT_DIR` entry at a directory of tile
mask images matching its `NUMBERING_SCHEME`. The healsparse-native replacement
for this whole procedure (an n_epoch / n_pointings survey-window map built from
the maps directly) is issue #797.
```

```{note}
Parts of this procedure use the legacy canfar-VM / `vos` retrieval workflow (see
[VOSpace retrieval](vos_retrieve.md)) and the obsolete `prepare_tiles_for_final`
helper, which is no longer shipped. The input-staging and joint-mask steps now
overlap with [`sp_validation`](https://github.com/CosmoStat/sp_validation). The
steps are retained for reference.
```

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

### Stage the pixel mask files

Collect the tile mask images into one directory that `random_cat_runner`'s
`INPUT_DIR` points at, one file per tile ID in `tile_numbers.txt`, named to
match the module's `FILE_PATTERN` and `NUMBERING_SCHEME` (`mask-<xxx>-<yyy>.fits`
with the committed `config_Rc.ini`). How you obtain them is outside ShapePipe;
older runs of the pipeline's own (now removed) mask module wrote them as
`pipeline_flag-<xxx>-<yyy>.fits`, and those files still work.

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


