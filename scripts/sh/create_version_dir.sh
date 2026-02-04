#!/bin/bash

# Create ShapePipe shear catalogue version directory with links and job script
# Usage: create_version_dir.sh <version> [year]
# Example: create_version_dir.sh v1.6.5 2024

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <version> [year]"
    echo "Example: $0 v1.6.5 2024"
    exit 1
fi

VERSION="$1"
YEAR="${2:-2024}"

# Extract major and minor version separately
# (e.g., from v1.6.5: MAJOR=1, MINOR=6, SUB=5)
ERSION=${VERSION#v}   # remove leading 'v' if present
IFS='.' read -r MAJOR MINOR SUB _ <<< "$ERSION"

echo "Creating dictionary for major=${MAJOR} minor=${MINOR} sub=${SUB}"

# Create directory
mkdir -p "$VERSION"
cd "$VERSION"

# Create link to HDF5 file (using X instead of MAJOR in link name)
ln -sf "../unions_shapepipe_comprehensive_struc_${YEAR}_v${MAJOR}.${MINOR}.c.hdf5" \
       "unions_shapepipe_comprehensive_struc_${YEAR}_v${MAJOR}.X.c.hdf5"

# Create link to mask config
ln -sf "/home/mkilbing/astro/repositories/github/sp_validation/config/calibration/mask_v${MAJOR}.X.${SUB}.yaml" \
       "config_mask.yaml"

# Copy and modify job script
JOB_SCRIPT="job_sbatch_calibrate.sh"
cp "/home/mkilbing/unions_wl/v1.4.x/v1.4.6/${JOB_SCRIPT}" .

# Get the full path for the cd command
FULL_PATH="/home/mkilbing/unions_wl/v1.${MINOR}.x/${VERSION}"

# Modify the cd line to point to the new directory (line has leading spaces)
sed -i "s|cd /home/mkilbing/unions_wl/v1.4.x/v1.4.6|cd ${FULL_PATH}|" "$JOB_SCRIPT"

# Modify the job name
sed -i "s|#SBATCH --job-name=.*|#SBATCH --job-name=cal_cc${MAJOR}${MINOR}|" "$JOB_SCRIPT"

echo "Created directory: $VERSION"
echo "Created link: unions_shapepipe_comprehensive_struc_${YEAR}_vX.${MINOR}.c.hdf5 -> ../unions_shapepipe_comprehensive_struc_${YEAR}_v${MAJOR}.${MINOR}.c.hdf5"
echo "Created link: config_mask.yaml"
echo "Copied and modified: $JOB_SCRIPT"
