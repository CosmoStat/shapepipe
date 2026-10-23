#! /usr/bin/env python3

"""COLLATE STAR CATALOGUES.

Collate the per-exposure PSFEx validation catalogues into star catalogues:
gather positions (X/Y/RA/DEC) and the HSM PSF and star shapes, tag each row
with its CCD number, and merge one ``validation_psf_conv-<idx>.fits`` per
single-exposure run.

HSM ellipticities and sizes are not rotated here. ``psfex_interp`` measures
adaptive moments directly in world coordinates (galsim
``FindAdaptiveMom(use_sky_coords=True)``), so the WCS-Jacobian shape rotation
this script used to perform is redundant and has been removed; only positions
are collated.

Caveat -- no workflow rule produces these catalogues yet. The input layout
below is the pre-workflow one: a link farm of ``psfex_interp`` run directories
under a single input root. The Snakemake workflow instead persists
``validation_psf-<exp>-<ccd>.fits`` into the per-exposure ``exp_persist``
tars on the products root, and ``merge_starcat_runner`` still wants the
collated ``validation_psf_conv-*`` files for the rho/tau statistics.
Repointing this script at the persist tars (or replacing it with a rule) is a
follow-up, and a bigger change than this strip.
"""

import sys
import os
import re
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import gc

import numpy as np
from astropy.io import fits

from cs_util import args as cs_args
from cs_util import logging


# Input layout: single-exposure run directories matching ``EXP_RUN_GLOB`` under
# ``<input_base_dir>/output``, each holding the PSF interpolation output in
# ``PSF_INTERP_SUBDIR``.
EXP_RUN_GLOB = "run_sp_combined_psf*"
PSF_INTERP_SUBDIR = "psfex_interp_runner/output"


def output_filename(file_pattern, idx):
    """Output Filename.

    Build the collated catalogue filename.

    Parameters
    ----------
    file_pattern : str
        input file pattern (e.g. ``validation_psf``)
    idx : int
        exposure run index

    Returns
    -------
    str
        output catalogue file name

    """
    return f"{file_pattern}_conv-{idx}.fits"


class Convert(object):

    def __init__(self):

        self.params_default()

    def set_params_from_command_line(self, args):
        """Set Params From Command line.

        Only use when calling using python from command line.
        Does not work from ipython or jupyter.

        """
        # Read command line options
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def params_default(self):

        self._params = {
            "input_base_dir": ".",
            "output_base_dir": ".",
            "mode": "merge",
            "file_pattern_psfint": "validation_psf",
        }

        self._short_options = {
            "input_base_dir": "-i",
            "mode": "-m",
        }

        self._types = {}

        self._help_strings = {
            "input_base_dir": (
                "input base dir; single-exposure runs are expected in"
                + " <input_base_dir>/output; default is {}"
            ),
            "mode": (
                "run mode, allowed are 'merge', 'test'; default is" + " '{}'"
            ),
        }

        # Output column names with types
        self._dt = [
            ("X", float),
            ("Y", float),
            ("RA", float),
            ("DEC", float),
            ("E1_PSF_HSM", float),
            ("E2_PSF_HSM", float),
            ("SIGMA_PSF_HSM", float),
            ("FLAG_PSF_HSM", float),
            ("E1_STAR_HSM", float),
            ("E2_STAR_HSM", float),
            ("SIGMA_STAR_HSM", float),
            ("FLAG_STAR_HSM", float),
            ("CCD_NB", int),
        ]

    def run(self):
        """Run.

        Main processing function.

        """
        do_parallel = True

        output_dir = self._params["output_base_dir"]
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        subdirs = f"{self._params['input_base_dir']}/output/{EXP_RUN_GLOB}"
        exp_run_dirs = glob.glob(subdirs)
        n_exp_runs = len(exp_run_dirs)
        print(f"Found {n_exp_runs} input single-exposure run(s) ({subdirs})")

        if self._params["mode"] == "test":
            exp_run_dirs = exp_run_dirs[:2]
            n_exp_runs = len(exp_run_dirs)
            print(
                f"test mode: only using {n_exp_runs} input single-exposure"
                + f" runs"
            )

        # Loop over exposure runs
        if not do_parallel:
            for idx_exp, exp_run_dir in tqdm(
                enumerate(exp_run_dirs),
                total=n_exp_runs,
                disable=self._params["verbose"],
            ):
                self.transform_exposures(output_dir, idx_exp, exp_run_dir)
        else:
            res = Parallel(n_jobs=-1, backend="loky")(
                delayed(self.transform_exposures)(
                    output_dir, idx_exp, exp_run_dir
                )
                for idx_exp, exp_run_dir in tqdm(
                    enumerate(exp_run_dirs),
                    total=n_exp_runs,
                    disable=self._params["verbose"],
                )
            )

    def transform_exposures(self, output_dir, idx, exp_run_dir):
        """Transform Exposures.

        Collate the PSF validation catalogues of one single-exposure run
        (input exp run dir) into one output catalogue.

        """
        output_path = (
            f"{output_dir}/"
            + output_filename(self._params["file_pattern_psfint"], idx)
        )
        if os.path.exists(output_path):
            print(f"Skipping transform_exposures, file {output_path} exists")
            return

        psf_dir = f"{exp_run_dir}/{PSF_INTERP_SUBDIR}"
        try:
            all_files = os.listdir(psf_dir)
            if self._params["verbose"]:
                print(f"Found {len(all_files)} file(s) in {psf_dir}")
        except Exception:
            if self._params["verbose"]:
                print(f"Found zero PSFEx files in {psf_dir}, skipping")
            return

        cat_list = []
        for file_name in all_files:
            if self._params["file_pattern_psfint"] not in file_name:
                continue

            tmp = re.findall(r"\d+", file_name)

            exp_name, ccd_id = int(tmp[0]), int(tmp[1])

            if self._params["verbose"]:
                print("Match found ", exp_name, ccd_id)

            psf_file_path = f"{psf_dir}/{file_name}"

            try:
                psf_file_hdus = fits.open(psf_file_path, memmap=False)
                psf_file = psf_file_hdus[2].data
                psf_file_hdus.close()
            except Exception:
                continue

            # HSM ellipticities and sizes are measured directly in world
            # coordinates upstream (FindAdaptiveMom use_sky_coords=True), so
            # they are passed straight through; only positions are collated.
            exp_cat = np.array(
                list(
                    map(
                        tuple,
                        np.array(
                            [
                                psf_file["X"],
                                psf_file["Y"],
                                psf_file["RA"],
                                psf_file["DEC"],
                                psf_file["E1_PSF_HSM"],
                                psf_file["E2_PSF_HSM"],
                                psf_file["SIGMA_PSF_HSM"],
                                psf_file["FLAG_PSF_HSM"],
                                psf_file["E1_STAR_HSM"],
                                psf_file["E2_STAR_HSM"],
                                psf_file["SIGMA_STAR_HSM"],
                                psf_file["FLAG_STAR_HSM"],
                                np.ones_like(psf_file["RA"], dtype=int)
                                * ccd_id,
                            ]
                        ).T.tolist(),
                    )
                ),
                dtype=self._dt,
            )
            cat_list.append(exp_cat)

            del psf_file

        if len(cat_list) == 0:
            return

        # Finalize catalogue
        star_cat = np.concatenate(cat_list)
        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(fits.BinTableHDU(star_cat))

        # Write catalogue
        hdul.writeto(
            output_path,
            overwrite=True,
        )

        del cat_list
        del hdul
        gc.collect()


def run_convert(*args):

    # Create instance
    obj = Convert()

    obj.set_params_from_command_line(args)

    obj.run()


def main(argv=None):
    """Main

    Main program
    """
    if argv is None:
        argv = sys.argv[1:]
    run_convert(*argv)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
