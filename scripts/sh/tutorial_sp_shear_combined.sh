
export SP_RUN=`pwd`
export SP_BASE=$HOME/astro/repositories/github/shapepipe
export SP_CONFIG=$SP_BASE/example/tutorial

mkdir output

# Pre-processing

## Select tile IDs to process
cfis_field_select -i $SP_BASE/aux/CFIS/tiles_202007/tiles_all_order.txt --coord 213.68deg_54.79deg -t tile --input_format ID_only --out_name_only --out_ID_only -s -o tile_numbers -v

## Retrieve tiles, uncompress tile weights, find single exposures, retrieve single exposures
shapepipe_run -c $SP_CONFIG/config_tile_GitUzFeGie.ini

# Processing of single exposures
shapepipe_run -c $SP_CONFIG/config_exp_SpMhMaSx.ini
