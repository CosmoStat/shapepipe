# -*- coding: utf-8 -*-

"""TILEOBJ_AS_EXP_SCRIPT

This script contains a class to handle processing for the combine_mexp_script
module: Combine information on objects and PSF from multiple
exposure-single-CCD tile-object catalogues.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: February 2019

:Version: 1.0

"""


import os
import re

import numpy as np

from astropy.io import fits
from astropy import wcs

from sf_tools.image.stamp import FetchStamps

import shapepipe.pipeline.file_io as io


class CombineMExpError(Exception):
   """ CombineMExpError

   Generic error that is raised by this script.

   """

   pass


class combine_mexp():

    """Class for combining multiple exposure-single-CCD catalogue and PSF information.

    Parameters
    ----------
    input_file_name: list of 2 string
        path to input files, [image, catalogue]
    config: config class
        config file content
    output_dir: string
        output directory
    image_number: string
        image number for input tile according to numbering scheme from config file
    w_log: log file class
        log file

    Returns
    -------
    None
    """

    def __init__(self, input_file_names, output_dir, image_number, config, w_log):

        self._img_tile_path = input_file_names[0]
        self._output_dir = output_dir
        self._image_number = image_number
        self._config = config
        self._w_log = w_log

        self._vignet = False
        self._do_check_consistency = True
        self._method = 'PSF_size'
        self._params = None

        self._write_all_mobj = False

    def process(self):

        # List of exposures contributing to this tile
        exp_list_uniq = self.get_exposure_list()

        # List of exposure-single-CCD file numbers, default format is
        # 000-0-00 (tile-exposure-hdu)
        sexp_number_list = self.get_sexp_number_list(len(exp_list_uniq))

        # Get detected objects from all input catalogues, and combine
        # into one array
        combined_data, n_exp_ok = self.combine_information(sexp_number_list)
        if combined_data:

            if self._write_all_mobj:
                # Write all selected objects to file (for testing)
                outcat_pattern = 'mobj_all'
                self.write_combined_data(combined_data, outcat_pattern)

            # Perform galaxy selection
            galaxies = self.select_galaxies(combined_data)

            # Write all selected galaxies
            outcat_pattern = 'mobj_gal'
            self.write_combined_data(galaxies, outcat_pattern)

            # For selected galaxies retrieve postage stamps and
            # corresponding IDs from exposure-single-CCD images
            vignets, IDs = self.get_vignets(galaxies, sexp_number_list)

            self.write_vignets(vignets, IDs)

            n_gal = len(galaxies['nexp'])
            n_tot = len(combined_data['ID'])

        else:

            n_gal = n_tot = 0

        self._w_log.info('Processed {} pairs of exposure-single-CCD and '
                         'galaxy PSF files'.format(n_exp_ok))
        self._w_log.info('{}/{} objects selected as galaxies'.format(n_gal, n_tot))

        return None, None

    def get_vignets(self, galaxies, sexp_number_list):
        """Return list of vignets corresponding to positions in list galaxies.

           See vignet_runner.py.

        Parameters
        ----------
        galaxies: FITS_rec
            multi-exposure galaxy catalogue data
        sexp_number_list: list of strings
            exposure-single-CCD file numbers

        Returns
        -------
        vignets: array of class FetchStamps
            list of vignets
        IDs: array of int
            list of object IDs corresponding to the postage stamps,
            including multiples.
        """

        img_basename = 'exp.img'
        input_dir_img = self._config.get('COMBINE_MEXP_RUNNER',
                                         'INPUT_DIR_IMG')
        img_ext = self._config.get('COMBINE_MEXP_RUNNER',
                                   'INPUT_IMG_EXT')

        stamp_size = config.getint("COMBINE_MEXP_RUNNER", "STAMP_SIZE") - 1
        if stamp_size % 2 != 0:
            raise ValueError("COMBINE_MEXP_RUNNER:STAMP_SIZE has to be an odd number.")
        half_size = int(stamp_size/2)

        # Loop over input file numbers, get all vignets
        vignets = []
        IDs = []
        for sexp_number in sexp_number_list:

            img_base_name = self.get_file_name(img_basename, sexp_number)
            img = self.open(input_dir_img, img_base_name, img_ext, hdu_no=0)
            header = img.get_header(hdu_no=0)
            wcs.WCS(header)

            # Get positions of galaxies from current exposure-single-CCD
            indices = np.where(galaxies['sepx_number'] == sexp_number)

            # see vignet_runner.py:get_pos
            ra = galaxies['X_WORLD'][indices]
            dec = galaxies['Y_WORLD'][indices]
            pox_xy = wcs.all_world2pix(ra, dec, 1)

            fs = FetchStamps(img, int(half_size))
            fs.get_pixels(np.round(pos).astype(int))
            vign = fs.scan()

            # Get galaxy IDs along with vignets
            vignets.add(vign)
            IDs.append(galaxies['ID'][indices])

        return vignets, IDs

    def combine_information(self, sexp_number_list):
        """Combine galaxy and PSF information from all input files.

        Parameters
        ----------
        sexp_number_list: list of strings
            exposure-single-CCD file numbers

        Returns
        -------
        combined_data: FITS_rec
            combined exposure-single-CCD and PSF data covering the tile
        n_exp_ok: int
            number of exposure-single-CCDs for which galaxy and PSF
            information contributed to combined_data
        """

        exp_basename = 'cat.exp.img'
        input_dir_exp_cat = self._config.get('COMBINE_MEXP_RUNNER',
                                             'INPUT_DIR_EXP_CAT')
        exp_cat_ext = self._config.get('COMBINE_MEXP_RUNNER',
                                       'INPUT_EXP_CAT_EXT')

        psf_basename = 'galaxy_psf'
        input_dir_psf = self._config.get('COMBINE_MEXP_RUNNER',
                                         'INPUT_DIR_PSF')
        psf_ext = self._config.get('COMBINE_MEXP_RUNNER', 'INPUT_PSF_EXT')

        combined_data = None

        n_exp_ok = 0
        for sexp_number in sexp_number_list:

            exp_base_name = self._get_file_name(exp_basename, sexp_number)
            exp_cat = self.open_or_continue(input_dir_exp_cat,
                                            exp_base_name,
                                            exp_cat_ext, hdu_no=2)
            if not exp_cat:
                continue

            psf_base_name = self._get_file_name(exp_basename, sexp_number)
            psf = self.open_or_continue(input_dir_psf, psf_base_name,
                                        psf_ext, hdu_no=2)
            if not psf:
                continue

            print('{}: both exposure-single CCD and PSF found'.format(
                  exp_base_name))

            exp_cat_data = exp_cat.get_data()
            exp_cat_cols = exp_cat.get_col_names()
            n_data = len(exp_cat_data)

            psf_data = psf.get_data()
            psf_cols = psf.get_col_names()

            if n_data != len(psf_data):
                raise CombineMExpError('Length of objects ({}) and PSF ({}) '
                                       'have to be the same'.
                                       format(n_data, len(psf_data)))

            exp_cat_data_plus = self.combine_exp_psf(exp_cat_data,
                                                     exp_cat_cols,
                                                     psf_data, psf_cols,
                                                     sexp_number)

            exp_cat.close()
            psf.close()

            if combined_data is None:
                combined_data = exp_cat_data_plus
            else:
                for c in exp_cat_data_plus:
                    combined_data[c] = np.concatenate((exp_cat_data_plus[c], combined_data[c]))

            n_exp_ok = n_exp_ok + 1

        return combined_data, n_exp_ok

    def write_vignets(vignets, IDs):
        """Save vignets and IDs to multi-HDU files.

        Parameters
        ----------
        vignets: array of class FetchStamps
            list of vignets
        IDs: array of int
            list of galaxy IDs corresponding to the postage stamps

        Returns
        -------
        None
        """

        IDs_unique, nexp = np.unique(IDs, return_counts=True)
        nexp_max = max(nexp_max)

        # Create nexp_max empty data arrays
        dat = []
        for m in range(nexp_max):
            dat_m = {}
            dat_m['ID'] = []
            dat_m['vignet'] = []
            dat.append(dat_m)

        # Loop over all unique IDs, write vignets of the nexp
        # contributing exposures to first nexp HDUs
        for u, n in zip(IDs_unique, nexp):
            indices = np.where(IDs == u)[0]

            # Append m-th ID and vignet to data array m
            for m in range(n):
                dat['ID'][m] = IDs[indices[m]]
                dat['vignet'][m] = vignets[incides[m]]

    def write_combined_data(self, galaxies, outcat_pattern):
        """Save data to file.

        Parameters
        ----------
        galaxies: FITS_rec
            catalogue data
        outcat_pattern: string
            output catalogue file name pattern

        Returns
        -------
        None
        """

        out_path = '{}/{}{}-0.fits'.format(self._output_dir, outcat_pattern, self._image_number)
        print('Writing tile data to file \'{}\''.format(out_path))
        output = io.FITSCatalog(out_path, open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        output.save_as_fits(galaxies)

    def select_galaxies(self, combined_data):
        """Perform galaxy selection.

        Parameters
        ----------
        combined_data: FITS_rec
            combined exposure-single-CCD and PSF data covering the tile

        Returns
        -------
        dat_gal: FITS_rec
            catalogue with unique galaxies (one entry per multi-exposure
            object)
        """

        # Create empty data array for selected galaxies
        dat_gal = {}
        for c in combined_data:
            dat_gal[c] = []

        # Add new column for number of exposures on which object is observed
        dat_gal['nexp'] = []

        # Find all objects with same ID and their frequency
        IDs = combined_data['ID']
        IDs_unique, nexp = np.unique(IDs, return_counts=True)

        n_is_gal = 0
        for u, n in zip(IDs_unique, nexp):
            indices = np.where(IDs == u)[0]

            self.check_consistency(indices, n)

            is_gal, params_out = self.select_one(combined_data, indices)
            if is_gal:
                n_is_gal += 1

                # Append information from first exposure only.
                # Galaxy info from all exposures is identical, since it came
                # from the same tile-selected object.
                # However, this is not the case for the PSF Info.
                # MKDEBUG TODO: Include PSF info from all exposures, here
                # or in a different module.
                for c in combined_data:
                    dat_gal[c].append(combined_data[c][indices[0]])
                dat_gal['nexp'].append(n)

            # TODO: param_out

        print(n_is_gal)

        return dat_gal

    def check_consistency(self, indices, nexp):
        """Check consistency of multi-exposure galaxy data

        Parameters
        ----------
        indices: numpy array of int
            list of indices
        nexp: int
            number of exposures

        Returns
        -------
        result: bool
            True if consistent
        """

        if not self._do_check_consistency:
            return True

        if len(indices) != nexp:
            raise CombineMExpError('Length of galaxy data {} != number of '
                                   'exposures {}'.format(len(indices), nexp))
        else:
            return True

    def select_one(self, data, indices):
        """Select galaxy according to multi-exposure information.

        Parameters
        ----------
        data: numpy array
            galaxy and PSF data
        indices: numpy array of int
            indices of corresponding to the same object from
            different exposures

        Returns
        -------
        is_gal, params_out: bool
            True if object is selected as galaxy
        """

        is_gal = False
        params_out = {}

        nexp = len(indices)

        if self._method == 'PSF_size':

            # fwhm from first exposure, identical to other exposures (from same tile object)
            fwhm = data['FWHM_IMAGE'][indices[0]]
            n_larger = 0
            for idx in indices:
                if data['HSM_FLAG'][idx] == 0 and fwhm > 2.355 * data['SIGMA_PSF_HSM'][idx]:
                    n_larger += 1
            params_out['n_larger'] = n_larger
            if n_larger == nexp:
                is_gal = True

        else:

            raise CombineMExpError('Invalid selection method \'{}\''
                                   .format(self._method))

        return is_gal, params_out

    def combine_exp_psf(self, exp_cat_data, exp_cat_cols, psf_data,
                        psf_cols, sexp_number):
        """Return data array with combined information from exposure-
           single-CCD and PSF.

        Parameters
        ----------
        exp_cat_data: FITS_rec
            exposure-single-CCD data
        exp_cols: array of string
            exposure-single-CCD column names
        psf_data: FITS_rec
            PSF data
        psf_cols: array of string
            PSF column names
        sexp_number: string
            exposure-single-CCD file number

        Returns
        -------
        exp_cat_data_plut: FITS_rec
            exposure-single-CCD and PSF joined data
        """

        exp_cat_data_plus = {}

        # Copy object data
        for c in exp_cat_cols:
            exp_cat_data_plus[c] = exp_cat_data[c]

        # Add PSF information to object catalogue
        for c in psf_cols:
            if self._vignet or c != 'VIGNET':
                # Add vignet only if argument 'vignet' is True
                exp_cat_data_plus[c] = psf_data[c]

        # Add exposure-single-CCD file number
        n = len(exp_cat_data_plus[exp_cat_cols[0]])
        exp_cat_data_plus['sexp_number'] = np.full(n, sexp_number)

        return exp_cat_data_plus

    def open_or_continue(self, input_dir, base_name, ext, hdu_no=0):
        """Return file content, or continue if file does not exist.

        Parameters
        ----------
        input_dir: string
            input directory
        base_name: string
            file base name
        ext: string
            file extension
        hdu_no: int, optional, default=0
            HDU number

        Returns
        -------
        cat: io.FITSCatalogue
            file content
        """

        try:
            cat = self.open(self, input_dir, base_name, ext, no_hdu=no_hdu)
        except:
            msg = 'File \'{}\' not found, continuing...'.format(cat_name)
            self._w_log.info(msg)
            cat = None

        return cat

    def open(self, input_dir, base_name, ext, hdu_no):
        """Return file content"

        Parameters
        ----------
        input_dir: string
            input directory
        base_name: string
            file base name
        ext: string
            file extension
        hdu_no: int, optional, default=0
            HDU number

        Returns
        -------
        cat: io.FITSCatalogue
            file content
        """

        cat_name = '{}/{}{}'.format(input_dir, base_name, ext)

        cat = io.FITSCatalog(cat_name, hdu_no=hdu_no)
        cat.open()
        return cat

    def get_exposure_list(self):
        """Return list of exposure file used for the tile in process, from tiles FITS header
           MKDEBUG: Note: this function is identical to the find_exposures class routine

        Parameters
        ----------

        Returns
        -------
        exp_list_uniq: list of string
            list of exposure basenames
        """

        try:
            hdu = fits.open(self._img_tile_path)
            hist = hdu[0].header['HISTORY']

        except:
            self._w_log.info('Error while reading tile catalogue FITS file {}, continuing...'.format(self._img_tile_path))

        tile_num = self._image_number

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise CombineMExpError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile #{}'.format(len(exp_list_uniq), tile_num))
        self._w_log.info('{} duplicates were removed'.format(len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq

    def get_file_list(self, n_exp):
        """Return lists of exposure-single-CCD file numbers.

        Parameters
        ----------
        n_exp: int
            number of exposures

        Returns
        -------
        sexp_number_list: list of strings
            exposure-single-CCD file numbers
        """

        n_hdu = int(self._config.get('COMBINE_MEXP_RUNNER', 'N_HDU'))
        sexp_number_list = []
        for i in range(n_exp):
            for k_img in range(n_hdu):
                sexp_number = '{}{}_{:02d}'.format(self._image_number, i, k_img)
                sexp_number_list.append(sexp_number)

        return sexp_number

    def get_file_name(basename, sexp_number):
        """ Return file name.

        Parameters
        ----------
        basename: string
            file base name
        sexp_number: string
            exposure-single-CCD file number

        Returns
        -------
        filename: string
            file name
        """

        return '{}{}'.format(basename, sexp_number)
