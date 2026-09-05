"""COVERAGE_RUN

Call coverage processing classes.

Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys

from shapepipe.utilities.header_downloader import HeaderDownloader
from shapepipe.utilities.field_corners_extractor import FieldCornersExtractor
from shapepipe.utilities.coverage_map_builder import CoverageMapBuilder
from shapepipe.utilities.coverage_plotter import CoveragePlotter


def run_download_headers(args=None):
    """Run Download Headers.

    Download FITS headers from VOSpace for exposures in a CCD list.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    obj = HeaderDownloader()
    return obj.run(args=args)


def run_extract_corners(args=None):
    """Run Extract Corners.

    Extract per-CCD sky-footprint corner coordinates from FITS headers.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    obj = FieldCornersExtractor()
    return obj.run(args=args)


def run_build_coverage(args=None):
    """Run Build Coverage.

    Build HealSparse coverage maps from per-CCD corner coordinates.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    obj = CoverageMapBuilder()
    return obj.run(args=args)


def run_plot_coverage(args=None):
    """Run Plot Coverage.

    Plot HealSparse coverage map.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    obj = CoveragePlotter()
    return obj.run(args=args)


def main(argv=None):
    """Main.

    Main program.

    """
    # Scripts to call coverage classes are created by pyproject.toml
    return 0
