#!/bin/bash

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
    INPUT_EXP="exp_ra_dec_${VERSION}.txt"
    COVERAGE_MAP="coverage_${VERSION}.x.hsp"

    # Build coverage map
    echo "  Building coverage map from ${INPUT_EXP}..."
    CMD="build_coverage_map -i ${INPUT_EXP} -o ${COVERAGE_MAP} -c ${BUILD_CHANNELS} -n ${BUILD_NSIDE} ${VERBOSE}"
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
