"""FAKE PSF.

Create fake PSF postage-stamp dictionaries for image simulations.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import pickle

import numpy as np
from astropy.io import fits
from sqlitedict import SqliteDict


class FakePsf:
    """Fake PSF.

    Parameters
    ----------
    sexcat_path : str
        Path to the tile SExtractor catalogue (multi-epoch FITS format).
    psf_dict_path : str
        Path to the pickled PSF dictionary.
    output_path : str
        Path for the output SqliteDict file.
    w_log : logging.Logger
        Pipeline logger.
    """

    def __init__(self, sexcat_path, psf_dict_path, output_path, w_log):
        self._sexcat_path = sexcat_path
        self._psf_dict_path = psf_dict_path
        self._output_path = output_path
        self._w_log = w_log

    def process(self):
        """Run fake PSF creation."""
        self._w_log.info(f"Reading sexcat: {self._sexcat_path}")
        try:
            sex = fits.open(self._sexcat_path, ignore_missing_simple=True)
        except Exception as e:
            self._w_log.error(f"Error opening FITS file {self._sexcat_path}: {e}")
            raise

        self._w_log.info(f"Catalogue has {len(sex)} HDUs")

        if len(sex) < 4:
            self._w_log.error(
                f"Catalogue file has insufficient HDUs: expected ≥4, got {len(sex)} "
                f"(file: {self._sexcat_path})"
            )
            raise ValueError(f"Invalid catalogue: only {len(sex)} HDUs (need ≥4)")

        self._w_log.info(f"Loading PSF dictionary: {self._psf_dict_path}")
        try:
            with open(self._psf_dict_path, "rb") as f:
                psf_dict = pickle.load(f)
        except Exception as e:
            self._w_log.error(f"Error loading PSF dictionary {self._psf_dict_path}: {e}")
            raise

        try:
            n_gal = len(sex[3].data.field("NUMBER"))
        except Exception as e:
            self._w_log.error(
                f"Error reading catalogue data from HDU 3 in {self._sexcat_path}: {e}"
            )
            raise
        n_exp = len(sex) - 3

        self._w_log.info(f"Processing {n_gal} galaxies over {n_exp} epochs")

        # Build (n_gal, n_exp) array of "expname-ccdnum" strings
        string_array = np.empty((n_gal, n_exp), dtype="object")
        for exp_ind in range(3, len(sex)):
            ccd_n = sex[exp_ind].data.field("CCD_N")
            exp_name = sex[exp_ind].data.field("EXP_NAME")
            string_array[:, exp_ind - 3] = [
                f"{int(exp)}-{int(ccd)}"
                for exp, ccd in zip(exp_name, ccd_n)
            ]

        # Mask invalid entries (exp or ccd == -99 produces strings with "--")
        mask = np.array(
            [["--" in item for item in row] for row in string_array]
        )
        masked = np.ma.masked_array(string_array, mask)

        # Build per-galaxy PSF dictionaries
        output_file = SqliteDict(self._output_path)
        missing = 0
        for idx, gal_row in enumerate(masked):
            galaxy_number = idx + 1  # 1-based, matches NUMBER field
            gal_dict = {}
            for exp_ccd in gal_row.compressed():
                if exp_ccd not in psf_dict:
                    missing += 1
                    self._w_log.warning(
                        f"Galaxy {galaxy_number}: key '{exp_ccd}' not in PSF dict"
                    )
                    continue
                gal_dict[exp_ccd] = psf_dict[exp_ccd]
            output_file[str(galaxy_number)] = gal_dict

        output_file.commit()
        output_file.close()
        sex.close()

        if missing:
            self._w_log.warning(f"{missing} missing PSF entries across all galaxies")
        self._w_log.info(f"Written: {self._output_path}")
