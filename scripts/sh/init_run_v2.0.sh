#!/usr/bin/env bash

# Name: init_run_v2.0.sh
# Description: Initialise the directory structure for a v2.0 ShapePipe run.
#              Creates tile and exposure staging directories, config symlinks,
#              and a tile number list from the 202604 tile catalogue.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>

# Version
version="2.0"

# Input type: data or image_sims
type="data"

# Subdir for image_sims
subdir="1z2z_grid_1"

force=0
sample=""

# Default base run directory (permanent storage)
#base_dir="$HOME/cosmostat/v2/v${version}"
base_dir=`pwd`

# ShapePipe repository root (for config symlink and tile list)
sp_root="$HOME/shapepipe"

# Source params.py to copy into the run directory
params_src="$HOME/astro/repositories/github/sp_validation/notebooks/params.py"

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\nOptions:\n
   -h\t\tthis message\n
   -t, --type TYPE input type, allowed are 'data', 'image_sims', default='$type'\n
   -s, --subdir SUBDIR subdir for image simulations, default='$subdir'\n
   -d, --dir DIR\tbase run directory, default='$base_dir'\n
   -P, --params PATH\tsource params.py to copy, default='$params_src'\n
   --force\t\trecreate existing symlinks and parameter files\n
   -S, --sample SAMPLE\tsample version for mask config (e.g. 6)\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -t|--type)
      type="$2"
      shift
      ;;
    -s|--subdir)
      subdir="$2"
      shift
      ;;
    -d|--dir)
      base_dir="$2"
      shift
      ;;
    -P|--params)
      params_src="$2"
      shift
      ;;
    --force)
      force=1
      ;;
    -S|--sample)
      sample="$2"
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done

# Check options
if [ "$type" == "data" ]; then

    # Config file directory
    config_dir="$sp_root/example/cfis"

    # Input tile list
    tiles_src="$sp_root/auxdir/CFIS/tiles_202604/tiles_r.txt"

elif [ "$type" == "image_sims" ]; then

    config_dir="$sp_root/example/cfis_image_sims"
    tiles_src="$sp_root/auxdir/CFIS/im_sims_202606/numbers.txt"

    input_dir_base="/n09data/hervas/skills_out"
    input_dir_tiles="$input_dir_base/$subdir/images/SP_tiles"
    input_dir_exp="$input_dir_base/$subdir/images/SP_exp"

else

    echo "Invalid input type $type"
    exit 3

fi


echo "Initialising ShapePipe v${version} run directory: $base_dir"
echo ""

# --- Base directory ---
mkdir -p "$base_dir/$subdir"
cd "$base_dir/$subdir"

if [ "$type" == "image_sims" ]; then
    for link in input_tiles input_exp; do
        src=$([ "$link" == "input_tiles" ] && echo "$input_dir_tiles" || echo "$input_dir_exp")
        [ $force -eq 1 ] && rm -f "$link"
        [ ! -e "$link" ] && ln -s "$src" "$link"
    done

    if [ -n "$sample" ]; then
        mask_src="$HOME/astro/repositories/github/sp_validation/config/calibration/mask_v1.X.${sample}_im_sim.yaml"
        [ $force -eq 1 ] && rm -f config_mask.yaml
        if [ ! -e config_mask.yaml ]; then
            ln -s "$mask_src" config_mask.yaml
            echo "Created symlink: config_mask.yaml -> $mask_src"
        else
            echo "config_mask.yaml already exists, skipping"
        fi
    fi
fi

echo "Creating tiles/ directory..."
mkdir -p tiles
echo "  exp/ created (subdirs e.g. tiles/301/301.279 added at download time)"

# --- Exposure staging directory ---
echo "Creating exp/ directory..."
mkdir -p exp
echo "  exp/ created (subdirs e.g. exp/28/ added at download time)"

# --- Output and working directories ---
echo "Creating log and debug directories..."
mkdir -p logs
mkdir -p debug

# --- Config symlink ---

# Config directory (will be symlinked)

if [ -d cfis ] && [ ! -L cfis ]; then
    echo "WARNING: cfis/ exists as a directory, not creating symlink"
else
    [ $force -eq 1 ] && rm -f cfis
    if [ ! -e cfis ]; then
        ln -s "$config_dir" cfis
        echo "Created symlink: cfis -> $config_dir"
    else
        echo "cfis symlink already exists, skipping"
    fi
fi

# --- Tile number list ---
echo "Creating tile_numbers.txt symlink..."
if [ -f tile_numbers.txt ] && [ ! -L tile_numbers.txt ]; then
    echo "  WARNING: tile_numbers.txt exists as a regular file, not creating symlink"
else
    [ $force -eq 1 ] && rm -f tile_numbers.txt
    if [ ! -e tile_numbers.txt ]; then
        ln -s "$tiles_src" tile_numbers.txt
        echo "  Created symlink: tile_numbers.txt -> $tiles_src"
    else
        echo "  tile_numbers.txt symlink already exists, skipping"
    fi
fi
n_tiles=$(wc -l < tile_numbers.txt)
echo "  $n_tiles tiles"

# --- params.py ---
echo "Creating params.py..."
if [ -f params.py ] && [ $force -eq 0 ]; then
    echo "  params.py already exists, skipping"
else
    cp "$params_src" params.py
    sed -i "s/[\"']IMAFLAGS_ISO[\"'],\s*//g; s/,\s*[\"']IMAFLAGS_ISO[\"']//g" params.py
    if [ "$type" == "image_sims" ]; then
        sed -i \
            -e "s/^name = .*/name = '$subdir'/" \
            -e "s/^star_cat_path = .*/star_cat_path = None/" \
            -e "s/^output_format = .*/output_format = '.hdf5'/" \
            params.py
        echo "  Copied and adapted: params.py <- $params_src"
        echo "    name = '$subdir', star_cat_path = None, output_format = '.hdf5'"
    else
        echo "  Copied: params.py <- $params_src"
    fi
fi

echo ""
echo "Done. Directory structure:"
echo "  $base_dir/$subdir"
echo "  ├── tiles/"
echo "  ├── exp/"
echo "  ├── logs/"
echo "  ├── cfis  ->  ${config_dir}"
echo "  ├── tile_numbers.txt  ->  ${tiles_src}"
echo "  └── params.py  <-  ${params_src}"
if [ "$type" == "image_sims" ]; then
    echo "  ___ input_dir_tiles -> $input_dir_tiles"
    echo "  ___ input_dir_exp -> $input_dir_exp"
fi
