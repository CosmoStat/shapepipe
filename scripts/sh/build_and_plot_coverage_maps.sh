#!/bin/bash

set -euo pipefail

# Build and plot per-CCD coverage (nexp) maps for a range of catalogue
# versions. Each map counts, per sky pixel, the number of exposures with a
# valid PSF model. The per-version flow is:
#   0. get_ccds_with_psf    -> ccds_with_psf_<version>.txt (valid CCD IDs)
#   1. download_headers     -> per-exposure header text files
#   2. extract_field_corners-> exp_ra_dec_<version>.txt (per-CCD corners)
#   3. build_coverage_map   -> coverage_<version>.x.hsp
#   4. plot_coverage_map    -> SGC / NGC plots
# Steps 0-2 are expected to have been run already on canfar (they need
# VOSpace access); this script rebuilds from the extracted per-CCD corners and
# plots. Uncomment the marked lines to run the full chain.

# Configuration variables

## Catalogue versions
VERSIONS=("v1.3" "v1.4" "v1.5" "v1.6")

## Output directory of plots
OUTPUT_DIR="coverage_map_plots"

# healsparse resolution. nside=131072 gives ~0.1" per pixel, chosen to
# match the UNIONS bit-mask resolution so coverage and mask align pixel-
# wise. CoverageMapBuilder defaults to nside=2048 for lighter-weight
# offline use; override here for the production build.
BUILD_NSIDE=131072
BUILD_CHANNELS=128

# Common parameters
VERBOSE="-v"

## Colorbar
PLOT_COLORBAR="-C"
PLOT_MIN=1
PLOT_MAX=5

# Plot parameters
## SGC region
SGC_RA_MIN=-20
SGC_RA_MAX=45
SGC_DEC_MIN=18
SGC_DEC_MAX=40

## NGC region
NGC_RA_MIN=110
NGC_RA_MAX=270
NGC_DEC_MIN=28
NGC_DEC_MAX=90

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Loop over versions
for VERSION in "${VERSIONS[@]}"; do
    echo "Processing ${VERSION}..."

    # Define file paths
    CCD_LIST="ccds_with_psf_${VERSION}.txt"
    HEADER_DIR="headers_${VERSION}"
    INPUT_CORNERS="exp_ra_dec_${VERSION}.txt"
    COVERAGE_MAP="coverage_${VERSION}.x.hsp"

    # Step 0-2: build the per-CCD corner file (canfar / VOSpace only).
    # Uncomment to run the full chain from scratch.
    # get_ccds_with_psf -V "${VERSION}" -o "${CCD_LIST}" ${VERBOSE}
    # download_headers -i "${CCD_LIST}" -o "${HEADER_DIR}" ${VERBOSE}
    # extract_field_corners -i "${HEADER_DIR}" -l "${CCD_LIST}" \
    #     -o "${INPUT_CORNERS}" ${VERBOSE}

    # Build coverage map
    echo "  Building coverage map from ${INPUT_CORNERS}..."
    CMD="build_coverage_map -i ${INPUT_CORNERS} -o ${COVERAGE_MAP} -c ${BUILD_CHANNELS} -n ${BUILD_NSIDE} ${VERBOSE}"
    echo "$CMD"
    $CMD

    # Plot SGC region
    echo "  Plotting SGC region..."
    CMD="plot_coverage_map -i ${COVERAGE_MAP} -o ${OUTPUT_DIR}/coverage_${VERSION}_SGC.png ${VERBOSE} ${PLOT_COLORBAR} -R ${SGC_RA_MIN} -r ${SGC_RA_MAX} -D ${SGC_DEC_MIN} -d ${SGC_DEC_MAX} -m ${PLOT_MIN} -M ${PLOT_MAX}"
    echo "$CMD"
    $CMD

    # Plot NGC region
    echo "  Plotting NGC region..."
    CMD="plot_coverage_map -i ${COVERAGE_MAP} -o ${OUTPUT_DIR}/coverage_${VERSION}_NGC.png ${VERBOSE} ${PLOT_COLORBAR} -R ${NGC_RA_MIN} -r ${NGC_RA_MAX} -D ${NGC_DEC_MIN} -d ${NGC_DEC_MAX} -m ${PLOT_MIN} -M ${PLOT_MAX}"
    echo "$CMD"
    $CMD

    echo "  Done with ${VERSION}"
    echo
done

echo "All versions processed successfully!"
