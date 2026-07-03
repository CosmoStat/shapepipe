#!/usr/bin/env bash

# Name: init_run_v2.0.sh
# Description: Initialise the directory structure for a v2.0 ShapePipe run.
#              Creates tile and exposure staging directories, config symlinks,
#              and a tile number list from the 202604 tile catalogue.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>

# Version
version="2.0"

# Default base run directory (permanent storage)
#base_dir="$HOME/cosmostat/v2/v${version}"
base_dir=`pwd`

# ShapePipe repository root (for config symlink and tile list)
sp_root="$HOME/shapepipe"

# Tile list source (full filenames, will be stripped to NNN.MMM)
tiles_src="$sp_root/auxdir/CFIS/tiles_202604/tiles_r.txt"

# Config directory (will be symlinked as $base_dir/cfis)
config_dir="$sp_root/example/cfis"

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\nOptions:\n
   -h\t\tthis message\n
   -d, --dir DIR\tbase run directory, default='$base_dir'\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -d|--dir)
      base_dir="$2"
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

echo "Initialising ShapePipe v${version} run directory: $base_dir"
echo ""

# --- Base directory ---
mkdir -p "$base_dir"
cd "$base_dir"

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
if [ -L cfis ]; then
    echo "cfis symlink already exists, skipping"
elif [ -d cfis ]; then
    echo "WARNING: cfis/ exists as a directory, not creating symlink"
else
    ln -s "$config_dir" cfis
    echo "Created symlink: cfis -> $config_dir"
fi

# --- Tile number list ---
echo "Creating tile_numbers.txt symlink..."
if [ -L tile_numbers.txt ]; then
    echo "  tile_numbers.txt symlink already exists, skipping"
elif [ -f tile_numbers.txt ]; then
    echo "  WARNING: tile_numbers.txt exists as a regular file, not creating symlink"
else
    ln -s "$tiles_src" tile_numbers.txt
    echo "  Created symlink: tile_numbers.txt -> $tiles_src"
fi
n_tiles=$(wc -l < tile_numbers.txt)
echo "  $n_tiles tiles"

echo ""
echo "Done. Directory structure:"
echo "  $base_dir/"
echo "  ├── tiles/"
echo "  ├── exp/"
echo "  ├── logs/"
echo "  ├── cfis  ->  ${config_dir}"
echo "  └── tile_numbers.txt  ->  ${tiles_src}"
