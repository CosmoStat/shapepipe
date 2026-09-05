"""COVERAGE_MAP_BUILDER

Build HealSparse coverage maps from per-CCD corner coordinates.

Each input row is one CCD footprint. Since the CCDs of a single exposure do
not overlap, stamping value 1 per CCD polygon and accumulating makes the map
count, per sky pixel, the number of exposures with a valid PSF model covering
that pixel.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
from os.path import exists

import numpy as np
import healsparse as hsp
import hpgeom as hpg

from cs_util import args as cs_args
from cs_util import logging

# Declination beyond which a polygon is considered too close to a pole for
# HealSparse's planar polygon fill to be reliable. CCD footprints are ~10
# arcmin, so this is a guard against pathological headers, not a common path.
_DEC_POLE_LIMIT = 89.0


def unwrap_ra(ra):
    """Unwrap RA.

    Shift RA values onto a common branch when a polygon straddles the RA=0
    seam, so the corners describe the small footprint rather than its ~360°
    complement. Values above 180 deg are shifted by -360 deg when the raw
    spread exceeds 180 deg; otherwise the values are returned unchanged.

    Parameters
    ----------
    ra : array_like
        polygon RA corners, in degrees

    Returns
    -------
    numpy.ndarray
        RA corners on a common branch, in degrees

    """
    ra = np.asarray(ra, dtype=float)
    if ra.max() - ra.min() > 180.0:
        return np.where(ra > 180.0, ra - 360.0, ra)
    return ra


def load_corners(input_file):
    """Read per-CCD corners from an ``exp_ra_dec.txt``.

    The nine-column text format the pre-Snakemake chain appended to: column 0 is
    a ``<expnum>-<ccd_idx>`` string ID, the remaining eight are the 4 RA then the
    4 Dec corners. The Snakemake workflow feeds :func:`build_map` arrays
    directly instead (one JSON record per exposure, written by the
    ``exp_footprint`` rule), so this loader exists only for the CLI.

    Parameters
    ----------
    input_file : str
        path to the per-CCD corner file

    Returns
    -------
    tuple
        ``(ccd_ids, ra, dec)``: a length-N string array and two ``(N, 4)``
        float arrays, in degrees

    """
    ccd_ids = np.atleast_1d(np.loadtxt(input_file, usecols=(0), dtype=str))
    corners = np.atleast_2d(
        np.loadtxt(input_file, usecols=range(1, 9), dtype=float)
    )
    return ccd_ids, corners[:, :4], corners[:, 4:]


def build_map(
    ccd_ids, ra, dec, nside_coverage, nside, verbose=False
):
    """Stamp one polygon per CCD footprint into a HealSparse nexp map.

    The array entry point, and the one both callers go through: the
    ``build_coverage_map`` CLI (which reads its arrays out of a text file with
    :func:`load_corners`) and the Snakemake ``coverage_map`` rule (which reads
    them out of the per-exposure ``exp_footprint.json`` records). Keeping the
    stamping in one place is what keeps the two ways of getting the same map
    from drifting — the RA-seam and pole guards below are survey geometry, not
    input-format handling.

    Since the CCDs of a single exposure do not overlap, accumulating value 1 per
    CCD polygon makes the map count, per sky pixel, the number of exposures with
    a valid PSF model covering it.

    Parameters
    ----------
    ccd_ids : array_like
        length-N CCD IDs, used only in the skipped-polygon warnings
    ra : array_like
        ``(N, 4)`` polygon RA corners, in degrees
    dec : array_like
        ``(N, 4)`` polygon Dec corners, in degrees
    nside_coverage : int
        HealSparse coverage nside
    nside : int
        HealSparse map nside
    verbose : bool, optional
        print progress

    Returns
    -------
    healsparse.HealSparseMap
        the nexp map

    """
    ra = np.atleast_2d(np.asarray(ra, dtype=float))
    dec = np.atleast_2d(np.asarray(dec, dtype=float))

    if verbose:
        print(
            f"Creating HealSparse map (nside_coverage={nside_coverage},"
            f" nside={nside})"
        )
    m = hsp.HealSparseMap.make_empty(nside_coverage, nside, np.uint16)

    if verbose:
        print("Adding polygons to map")

    n_added = 0
    n_skipped = 0
    for i in range(len(ccd_ids)):
        dec_i = dec[i]

        # Pole guard: HealSparse's planar polygon fill degrades near the
        # poles. CCD footprints never reach here, so warn and skip.
        if np.any(np.abs(dec_i) >= _DEC_POLE_LIMIT):
            print(
                f"Warning: skipping CCD {ccd_ids[i]} with |dec| >= "
                f"{_DEC_POLE_LIMIT} (too close to a pole)"
            )
            n_skipped += 1
            continue

        # RA-wrap guard: put corners on a common branch across the seam.
        ra_i = unwrap_ra(ra[i])

        m += hsp.Polygon(ra=list(ra_i), dec=list(dec_i), value=1)
        n_added += 1

        if verbose and i % 1000 == 0:
            print(f"{i:6d} / {len(ccd_ids):6d}")

    print(f"Added {n_added} polygons to map")
    if n_skipped > 0:
        print(f"Skipped {n_skipped} polygons near a pole")

    return m


class CoverageMapBuilder(object):
    """Coverage Map Builder Class.

    Builds HealSparse coverage maps from field corner coordinates.
    """

    def __init__(self):
        """Initialize the builder."""
        self.params_default()

    def params_default(self):
        """Set default parameters and command line options."""

        self._params = {
            "input_file": "exp_ra_dec.txt",
            "output_file": "coverage.hsp",
            "nside_coverage": 32,
            "nside": 2048,
            "apply_median_filter": False,
            "n_median_iterations": 2,
            "create_boolean": False,
            "boolean_threshold": 3,
            "boolean_output": None,
            "create_plot": False,
            "plot_output": None,
            "plot_region": None,
            "verbose": False,
        }

        self._short_options = {
            "input_file": "-i",
            "output_file": "-o",
            "nside_coverage": "-c",
            "nside": "-n",
            "apply_median_filter": "-m",
            "n_median_iterations": "-N",
            "create_boolean": "-b",
            "boolean_threshold": "-t",
            "boolean_output": "-B",
            "create_plot": "-p",
            "plot_output": "-P",
            "plot_region": "-g",
        }

        self._types = {
            "nside_coverage": "int",
            "nside": "int",
            "apply_median_filter": "bool",
            "n_median_iterations": "int",
            "create_boolean": "bool",
            "boolean_threshold": "int",
            "create_plot": "bool",
            "verbose": "bool",
        }

        self._help_strings = {
            "input_file": "input file with field corners; default is {}",
            "output_file": "output HealSparse map file; default is {}",
            "nside_coverage": "HealSparse coverage nside; default is {}",
            "nside": "HealSparse map nside; default is {}",
            "apply_median_filter": "apply median filter to smooth map; default is {}",
            "n_median_iterations": "number of median filter iterations; default is {}",
            "create_boolean": "create boolean coverage map; default is {}",
            "boolean_threshold": "threshold value for boolean map; default is {}",
            "boolean_output": "output file for boolean map; default is <output_file>_bool.hsp",
            "create_plot": "create plot of coverage map; default is {}",
            "plot_output": "output file for plot; default is <output_file>.png",
            "plot_region": "predefined region for plot (NGC, SGC, fullsky); default is {}",
        }

    def set_params_from_command_line(self, args):
        """Set Params From Command line.

        Only use when calling using python from command line.
        Does not work from ipython or jupyter.

        Parameters
        ----------
        args : list
            command line arguments

        """
        # Read command line options
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
            args=args,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def update_params(self):
        """Update parameters.

        Set derived parameters based on input parameters.
        """
        # Set boolean output filename if not specified
        if self._params["boolean_output"] is None:
            output_file = self._params["output_file"]
            if output_file.endswith(".hsp"):
                base = output_file[:-4]
            else:
                base = output_file
            self._params["boolean_output"] = f"{base}_bool.hsp"

        # Set plot output filename if not specified
        if self._params["plot_output"] is None:
            output_file = self._params["output_file"]
            if output_file.endswith(".hsp"):
                base = output_file[:-4]
            else:
                base = output_file
            self._params["plot_output"] = f"{base}.png"

    def check_params(self):
        """Check parameters for validity."""
        if not exists(self._params["input_file"]):
            raise FileNotFoundError(
                f"Input file not found: {self._params['input_file']}"
            )

        # Check nside values are powers of 2
        nside_coverage = self._params["nside_coverage"]
        nside = self._params["nside"]

        if not (nside_coverage & (nside_coverage - 1) == 0):
            raise ValueError(
                f"nside_coverage must be a power of 2, got {nside_coverage}"
            )

        if not (nside & (nside - 1) == 0):
            raise ValueError(f"nside must be a power of 2, got {nside}")

    def median_filter(self, hsp_map):
        """Median Filter.

        Apply median filter to HealSparse map using neighbors.

        Parameters
        ----------
        hsp_map : healsparse.HealSparseMap
            input map

        Returns
        -------
        healsparse.HealSparseMap
            filtered map

        """
        nside = hsp_map.nside_sparse

        new_hsp = hsp_map.copy()
        pixs = hsp_map.valid_pixels
        n = hpg.neighbors(nside, pixs)
        n = np.where(n >= 0, n, 0)
        pix2 = pixs.reshape(len(pixs), 1)
        n = np.hstack([n, pix2])
        mn = hsp_map[n]
        new = np.median(mn, axis=1).astype(np.uint16)
        new_hsp[pixs] = new

        return new_hsp

    def run(self, args=None):
        """Run.

        Main execution method.

        Parameters
        ----------
        args : list, optional
            command line arguments

        Returns
        -------
        int
            exit code (0 for success)

        """
        if args is None:
            args = sys.argv[1:]

        # Set parameters from command line
        self.set_params_from_command_line(args)
        self.update_params()
        self.check_params()

        # Get parameters
        input_file = self._params["input_file"]
        output_file = self._params["output_file"]
        nside_coverage = self._params["nside_coverage"]
        nside = self._params["nside"]
        apply_median_filter = self._params["apply_median_filter"]
        n_median_iterations = self._params["n_median_iterations"]
        create_boolean = self._params["create_boolean"]
        boolean_threshold = self._params["boolean_threshold"]
        boolean_output = self._params["boolean_output"]
        create_plot = self._params["create_plot"]
        plot_output = self._params["plot_output"]
        plot_region = self._params["plot_region"]
        verbose = self._params["verbose"]

        if verbose:
            print(f"Reading per-CCD corners from: {input_file}")

        ccd_ids, ra_all, dec_all = load_corners(input_file)

        print(f"Loaded {len(ccd_ids)} CCD footprints")

        m = build_map(
            ccd_ids, ra_all, dec_all, nside_coverage, nside, verbose=verbose
        )

        # Apply median filter if requested
        if apply_median_filter:
            if verbose:
                print(
                    f"Applying median filter ({n_median_iterations} iterations)"
                )

            for i in range(n_median_iterations):
                m = self.median_filter(m)
                if verbose:
                    print(f"  Iteration {i+1}/{n_median_iterations} complete")

        # Write main coverage map
        if verbose:
            print(f"Writing coverage map to: {output_file}")

        m.write(output_file, clobber=True)
        print(f"Coverage map saved to {output_file}")

        # Create and write boolean map if requested
        if create_boolean:
            if verbose:
                print(
                    f"Creating boolean map (threshold={boolean_threshold})"
                )

            c = hsp.HealSparseMap.make_empty(nside_coverage, nside, bool)
            c[m.valid_pixels] = m[m.valid_pixels] >= boolean_threshold

            if verbose:
                print(f"Writing boolean map to: {boolean_output}")

            c.write(boolean_output, clobber=True)
            print(f"Boolean map saved to {boolean_output}")

        # Create plot if requested
        if create_plot:
            if verbose:
                print(f"Creating plot of coverage map")

            try:
                from shapepipe.utilities.coverage_plotter import CoveragePlotter
                plotter = CoveragePlotter()
                plotter.plot_coverage_map(
                    m,
                    plot_output,
                    region=plot_region,
                    vmax=boolean_threshold if create_boolean else 3,
                    colorbar=True,
                    colorbar_label="Coverage depth",
                )
                print(f"Plot saved to {plot_output}")
            except ImportError as e:
                print(f"Warning: Could not create plot: {e}")
                print("Install the cs_util package for plotting support.")

        return 0
