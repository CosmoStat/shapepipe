"""COVERAGE_PLOTTER

Plot HealSparse coverage maps.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
from os.path import exists

import numpy as np
import healsparse as hsp
import matplotlib
import matplotlib.pyplot as plt

try:
    from cs_util.plots import FootprintPlotter
    HAS_FOOTPRINT_PLOTTER = True
except ImportError:
    HAS_FOOTPRINT_PLOTTER = False

from cs_util import args as cs_args
from cs_util import logging


class CoveragePlotter(object):
    """Coverage Plotter Class.

    Plots HealSparse coverage maps.
    """

    def __init__(self):
        """Initialize the plotter."""
        self.params_default()

    def params_default(self):
        """Set default parameters and command line options."""

        self._params = {
            "input_file": None,
            "output_file": "coverage_plot.png",
            "region": None,
            "ralo": 100.0,
            "rahi": 270.0,
            "declo": 28.0,
            "dechi": 85.0,
            "vmin": 0,
            "vmax": 3,
            "cmap": "rainbow",
            "draw_milky_way": True,
            "colorbar": True,
            "colorbar_label": "Coverage depth",
            "figsize_x": 12,
            "figsize_y": 8,
            "dpi": 200,
        }

        self._short_options = {
            "input_file": "-i",
            "output_file": "-o",
            "region": "-g",
            "ralo": "-R",
            "rahi": "-r",
            "declo": "-D",
            "dechi": "-d",
            "vmin": "-m",
            "vmax": "-M",
            "cmap": "-c",
            "draw_milky_way": "-w",
            "colorbar": "-C",
            "colorbar_label": "-L",
            "figsize_x": "-x",
            "figsize_y": "-y",
            "dpi": "-p",
        }

        self._types = {
            "ralo": "float",
            "rahi": "float",
            "declo": "float",
            "dechi": "float",
            "vmin": "int",
            "vmax": "int",
            "draw_milky_way": "bool",
            "colorbar": "bool",
            "figsize_x": "float",
            "figsize_y": "float",
            "dpi": "int",
        }

        self._help_strings = {
            "input_file": "input HealSparse map file; required",
            "output_file": "output plot file (png, pdf, etc.); default is {}",
            "region": "predefined region (NGC, SGC, fullsky); overrides RA/Dec limits if specified; default is {}",
            "ralo": "minimum RA for plot extent; default is {}",
            "rahi": "maximum RA for plot extent; default is {}",
            "declo": "minimum Dec for plot extent; default is {}",
            "dechi": "maximum Dec for plot extent; default is {}",
            "vmin": "minimum value for colormap; default is {}",
            "vmax": "maximum value for colormap; default is {}",
            "cmap": "matplotlib colormap name; default is {}",
            "draw_milky_way": "draw Milky Way outline; default is {}",
            "colorbar": "add colorbar to plot; default is {}",
            "colorbar_label": "colorbar label; default is {}",
            "figsize_x": "figure width in inches; default is {}",
            "figsize_y": "figure height in inches; default is {}",
            "dpi": "figure DPI; default is {}",
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
        pass

    def check_params(self):
        """Check parameters for validity."""
        if self._params["input_file"] is None:
            raise ValueError("Input file is required (use -i or --input_file)")

        if not exists(self._params["input_file"]):
            raise FileNotFoundError(
                f"Input file not found: {self._params['input_file']}"
            )

        if not HAS_FOOTPRINT_PLOTTER:
            raise ImportError(
                "cs_util package with FootprintPlotter is required for plotting. "
                "Install with: pip install -e /path/to/cs_util"
            )

        # Check region validity
        region = self._params["region"]
        if region is not None:
            valid_regions = list(FootprintPlotter._regions.keys())
            if region not in valid_regions:
                raise ValueError(
                    f"Invalid region '{region}'. Valid regions are: {', '.join(valid_regions)}"
                )

    def plot_coverage_map(
        self,
        hsp_map,
        output_file,
        region=None,
        ralo=100.0,
        rahi=270.0,
        declo=28.0,
        dechi=85.0,
        vmin=0,
        vmax=3,
        title=None,
        colorbar=True,
        colorbar_label="Coverage depth",
        figsize=(10, 10),
        dpi=200,
    ):
        """Plot Coverage Map.

        Create a plot of the HealSparse coverage map using FootprintPlotter.

        Parameters
        ----------
        hsp_map : healsparse.HealSparseMap
            coverage map to plot
        output_file : str
            output file path
        region : str, optional
            predefined region name (NGC, SGC, fullsky)
        ralo : float
            minimum RA for plot extent (used if region is None)
        rahi : float
            maximum RA for plot extent (used if region is None)
        declo : float
            minimum Dec for plot extent (used if region is None)
        dechi : float
            maximum Dec for plot extent (used if region is None)
        vmin : int
            minimum value for colormap
        vmax : int
            maximum value for colormap
        title : str, optional
            plot title
        colorbar : bool
            add colorbar to plot
        colorbar_label : str
            colorbar label
        figsize : tuple
            figure size (width, height) in inches
        dpi : int
            figure DPI

        """
        # Set up matplotlib parameters
        font = {"size": 16}
        matplotlib.rc("font", **font)
        matplotlib.rcParams["savefig.dpi"] = dpi
        matplotlib.rcParams["figure.dpi"] = dpi
        matplotlib.rcParams["figure.figsize"] = figsize

        # Create FootprintPlotter instance
        plotter = FootprintPlotter(
            nside_coverage=hsp_map.nside_coverage,
            nside_map=hsp_map.nside_sparse
        )

        # Plot based on region or manual extent
        if region is not None:
            # Use predefined region
            region_params = FootprintPlotter._regions[region]
            plotter.plot_region(
                hsp_map,
                region_params,
                outpath=output_file,
                title=title,
                colorbar=colorbar,
                colorbar_label=colorbar_label
            )
        else:
            # Use manual RA/Dec limits
            # Calculate ra_0 for projection
            if ralo < rahi:
                ra_0 = (ralo + rahi) / 2.0
            else:
                ra_0 = 130.0

            extend = [ralo, rahi, declo, dechi]
            plotter.plot_area(
                hsp_map,
                ra_0=ra_0,
                extend=extend,
                vmax=vmax,
                outpath=output_file,
                title=title,
                colorbar=colorbar,
                colorbar_label=colorbar_label
            )

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
            args = sys.argv

        # Set parameters from command line
        self.set_params_from_command_line(args)
        self.update_params()
        self.check_params()

        # Get parameters
        input_file = self._params["input_file"]
        output_file = self._params["output_file"]
        region = self._params["region"]
        ralo = self._params["ralo"]
        rahi = self._params["rahi"]
        declo = self._params["declo"]
        dechi = self._params["dechi"]
        vmin = self._params["vmin"]
        vmax = self._params["vmax"]
        colorbar = self._params["colorbar"]
        colorbar_label = self._params["colorbar_label"]
        figsize_x = self._params["figsize_x"]
        figsize_y = self._params["figsize_y"]
        dpi = self._params["dpi"]
        verbose = self._params["verbose"]

        if verbose:
            print(f"Reading coverage map from: {input_file}")

        # Load HealSparse map
        hsp_map = hsp.HealSparseMap.read(input_file)

        print(f"Loaded map with nside={hsp_map.nside_sparse}")
        print(
            f"Valid pixels: {len(hsp_map.valid_pixels)} / {hsp_map.nside_sparse**2 * 12}"
        )

        if verbose:
            if region:
                print(f"Creating plot for region: {region}")
            else:
                print(f"Creating plot with extent: RA=[{ralo}, {rahi}], Dec=[{declo}, {dechi}]")

        # Create plot
        self.plot_coverage_map(
            hsp_map,
            output_file,
            region=region,
            ralo=ralo,
            rahi=rahi,
            declo=declo,
            dechi=dechi,
            vmin=vmin,
            vmax=vmax,
            colorbar=colorbar,
            colorbar_label=colorbar_label,
            figsize=(figsize_x, figsize_y),
            dpi=dpi,
        )

        print(f"Plot saved to {output_file}")

        return 0
