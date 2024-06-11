#! /usr/bin/env python3

import sys
import os
import re
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import gc

import numpy as np
from astropy.io import fits
import galsim

from cs_util import args as cs_args
from cs_util import logging


def transform_shape(mom_list, jac):
    scale, shear, theta, flip = jac.getDecomposition()

    sig_tmp = mom_list[2] * scale
    shape = galsim.Shear(g1=mom_list[0], g2=mom_list[1])
    if flip:
        print("FLIP!")
        shape = galsim.Shear(g1=-shape.g1, g2=shape.g2)
    shape = galsim.Shear(g=shape.g, beta=shape.beta + theta)
    shape = shear + shape

    return shape.g1, shape.g2, sig_tmp


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
            "sub_dir_pattern" : "run_sp_exp_202",
            "file_pattern_psfint": "validation_psf",
            "sub_dir_psfint": "psfex_interp_runner/output/",
        }

        self._short_options = {
            "input_base_dir": "-i",
            "mode": "-m",
            "sub_dir_pattern" : "-s",
        }

        self._types = {}

        self._help_strings = {
            "input_base_dir" : (
                "input base dir, runs are expected in"
                + " <input_base_dir>/P<patch?>/output;"
                " default is {}"
            ),
            "sub_dir_pattern" : (
                "First subdir name of PSF interpolation, default is {}"
            ),
            "mode" : (
                "run mode, allowed are 'create', 'merge', 'test'; default is"
                + " '{}':"
            )
        }

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

        if self._params["mode"] == "test":
            patch_nums = ["3", "4"]
        else:
            n_patch = 7                                                             
            #patch_nums = [idx for idx in np.arange(n_patch) + 1]
            patch_nums = [1, 3, 4]

        do_parallel = True

        # Loop over patches
        for patch in patch_nums:

            output_dir = f"{self._params['output_base_dir']}/P{patch}"
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir, exist_ok=False)

            print("Running patch:", patch)

            patch_dir = f"{self._params['input_base_dir']}/P{patch}/output/"
            subdirs = f"{patch_dir}/{self._params['sub_dir_pattern']}*"
            exp_run_dirs = glob.glob(subdirs)
            n_exp_runs = len(exp_run_dirs)
            print(
                f"Found {n_exp_runs} input single-exposure run(s) for patch"
                + f" {patch_dir} ({subdirs})"
            )

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
                        enumerate(exp_run_dirs), total=n_exp_runs,
                        disable=self._params["verbose"],
                ):
                    self.transform_exposure(
                        output_dir, patch, idx_exp, exp_run_dir
                    )
            else:
                res = Parallel(n_jobs=-1, backend="loky")(
                    delayed(self.transform_exposure)(
                        output_dir, patch, idx_exp, exp_run_dir
                    )
                    for idx_exp, exp_run_dir in tqdm(
                            enumerate(exp_run_dirs),
                            total=n_exp_runs,
                            disable=self._params["verbose"],
                    )
                )

    def transform_exposure(self, output_dir, patch, idx, exp_run_dir):
        """Transform exposures.

        Transform shapes for exposures for a given run (input exp run dir).

        """
        output_path = (
            f"{output_dir}/{self._params['file_pattern_psfint']}_conv"
            + f"-{patch}-{idx}.fits"
        )
        if os.path.exists(output_path):
            print(f"Skipping transform_exposure, file {output_path} exists")
            return

        mccd_dir = f"{exp_run_dir}/{self._params['sub_dir_psfint']}"
        try:
            all_files = os.listdir(mccd_dir)
            if self._params["verbose"]:
                print(f"Found {len(all_files)} file(s) in {mccd_dir}")
        except Exception:
            if self._params["verbose"]:
                print(f"Found zero PSFEx files in {mccd_dir}, skipping")
            return

        cat_list = []
        for file_name in all_files:
            if self._params["file_pattern_psfint"] not in file_name:
                continue

            tmp = re.findall(r"\d+", file_name)
            exp_name, ccd_id = int(tmp[0]), int(tmp[1])
            if self._params["verbose"]:
                print("Match found ", exp_name, ccd_id)

            psfex_file_path = f"{mccd_dir}/{file_name}"

            try:
                psf_file = fits.open(psfex_file_path, memmap=False)
                psfex_file = psf_file[2].data
                header_file = psf_file[1].data
            except Exception:
                continue

            new_e1_psf = np.zeros_like(psfex_file["RA"])
            new_e2_psf = np.zeros_like(psfex_file["RA"])
            new_sig_psf = np.zeros_like(psfex_file["RA"])
            new_e1_star = np.zeros_like(psfex_file["RA"])
            new_e2_star = np.zeros_like(psfex_file["RA"])
            new_sig_star = np.zeros_like(psfex_file["RA"])

            header = fits.Header.fromstring(
                "\n".join(header_file[0][0]), sep="\n"
            )
            wcs = galsim.AstropyWCS(header=header)

            k = 0
            for ind, obj in enumerate(psfex_file):
                try:
                    # jac = wcs.jacobian(world_pos=galsim.CelestialCoord(
                    #     ra=obj["RA"]*galsim.degrees,
                    #     dec=obj["DEC"]*galsim.degrees
                    # ))
                    jac = wcs.jacobian(
                        image_pos=galsim.PositionD(
                            obj["X"],
                            obj["Y"],
                        )
                    )
                except Exception:
                    continue
                g1_psf_tmp, g2_psf_tmp, sig_psf_tmp = transform_shape(
                    [
                        obj["E1_PSF_HSM"],
                        obj["E2_PSF_HSM"],
                        obj["SIGMA_PSF_HSM"],
                    ],
                    jac,
                )

                new_e1_psf[ind] = g1_psf_tmp
                new_e2_psf[ind] = g2_psf_tmp
                new_sig_psf[ind] = sig_psf_tmp

                g1_star_tmp, g2_star_tmp, sig_star_tmp = transform_shape(
                    [
                        obj["E1_STAR_HSM"],
                        obj["E2_STAR_HSM"],
                        obj["SIGMA_STAR_HSM"],
                    ],
                    jac,
                )
                new_e1_star[ind] = g1_star_tmp
                new_e2_star[ind] = g2_star_tmp
                new_sig_star[ind] = sig_star_tmp
                k += 1

            exp_cat = np.array(
                list(
                    map(
                        tuple,
                        np.array(
                            [
                                psfex_file["X"],
                                psfex_file["Y"],
                                psfex_file["RA"],
                                psfex_file["DEC"],
                                new_e1_psf,
                                new_e2_psf,
                                new_sig_psf,
                                psfex_file["FLAG_PSF_HSM"],
                                new_e1_star,
                                new_e2_star,
                                new_sig_star,
                                psfex_file["FLAG_STAR_HSM"],
                                np.ones_like(psfex_file["RA"], dtype=int)
                                * ccd_id,
                            ]
                        ).T.tolist(),
                    )
                ),
                dtype=self._dt,
            )
            cat_list.append(exp_cat)

            psf_file.close()
            del psf_file
            gc.collect()

        if len(cat_list) == 0:
            return

        # Finalize catalogue
        patch_cat = np.concatenate(cat_list)
        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(fits.BinTableHDU(patch_cat))

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
