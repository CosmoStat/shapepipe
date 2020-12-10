export SP_RUN=`pwd`
export SP_BASE=$HOME/astro/repositories/github/shapepipe
export SP_CONFIG=$SP_BASE/example/tutorial

mkdir output

# Pre-processing

## Select tile IDs to process
cfis_field_select -i $SP_BASE/aux/CFIS/tiles_202007/tiles_all_order.txt --coord 213.68deg_54.79deg -t tile --input_format ID_only --out_name_only --out_ID_only -s -o tile_numbers -v

## Retrieve tiles
shapepipe_run -c $SP_CONFIG/config_tile_Git.ini

## Uncompress tile weights
shapepipe_run -c $SP_CONFIG/config_tile_Uz.ini

## Find single exposures
shapepipe_run -c $SP_CONFIG/config_tile_Fe.ini

## Retrieve single exposures
shapepipe_run -c $SP_CONFIG/config_tile_Gie.ini

# Processing of single exposures

## Split into single-exposure single-HDU files and write FITS headers
shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini

## Merge FITS headers
shapepipe_run -c $SP_CONFIG/config_exp_Mh.ini

## Mask images
shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini

## Detect objects
shapepipe_run -c $SP_CONFIG/config_exp_Sx.ini

## Select star candidates
shapepipe_run -c $SP_CONFIG/config_exp_Se.ini

## Validation: Create stats plots
stats_global -o stats -v -c $SP_CONFIG/config_stats.ini

## Create PSF model
shapepipe_run -c $SP_CONFIG/config_exp_Psm.ini

## Interpolate PSF model to star positions (for validation)
shapepipe_run -c $SP_CONFIG/config_exp_Psi.ini

## Validation: PSF residuals
shapepipe_run -c $SP_CONFIG/config_exp_Cp.ini

mkdir psf_validation
shapepipe_run -c $SP_CONFIG/config_exp_Mst.ini
