"""COVERAGE_MAP_BUILDER

Build HealSparse coverage maps from field corner coordinates.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
from os.path import exists

import numpy as np
import healsparse as hsp
import hpgeom as hpg

from cs_util import args as cs_args
from cs_util import logging



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
            print(f"Reading field corners from: {input_file}")

        # Load field corner data
        dtype = np.dtype(
            [
                ("expnum", int),
                ("tlra", float),
                ("trra", float),
                ("brra", float),
                ("blra", float),
                ("tldec", float),
                ("trdec", float),
                ("brdec", float),
                ("bldec", float),
            ]
        )
        rows = np.loadtxt(input_file, dtype=dtype)

        print(f"Loaded {len(rows)} exposure field corners")

        if verbose:
            print(
                f"Creating HealSparse map (nside_coverage={nside_coverage}, nside={nside})"
            )

        # Create empty map
        m = hsp.HealSparseMap.make_empty(nside_coverage, nside, np.uint16)

        # Add polygons for each exposure
        if verbose:
            print("Adding polygons to map")

        for i, r in enumerate(rows):
            ra = [r["tlra"], r["trra"], r["brra"], r["blra"]]
            dec = [r["tldec"], r["trdec"], r["brdec"], r["bldec"]]
            poly = hsp.Polygon(ra=ra, dec=dec, value=1)
            m += poly

            if verbose and i % 1000 == 0:
                print(f"{i:6d} / {len(rows):6d}")

        print(f"Added {len(rows)} polygons to map")

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
