"""PSFEX INTERPOLATION SCRIPT.

This module computes the PSFs from a PSFEx model at several galaxy positions.

:Authors: Morgan Schmitz and Axel Guinot

"""

import os
import re

import numpy as np
from astropy.io import fits
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io

try:
    import galsim.hsm as hsm
    from galsim import Image
except ImportError:
    import_fail = True
else:
    import_fail = False


NOT_ENOUGH_STARS = 'Fail_stars'
BAD_CHI2 = 'Fail_chi2'
FILE_NOT_FOUND = 'File_not_found'


def interpsfex(dotpsfpath, pos, thresh_star, thresh_chi2):
    """Interpolate PSFEx.

    Use PSFEx generated model to perform spatial PSF interpolation.

    Parameters
    ----------
    dotpsfpath : str
        Path to ``.psf`` file (PSFEx output)
    pos : numpy.ndarray
        Positions where the PSF model should be evaluated
    thresh_star : int
        Threshold of stars under which the PSF is not interpolated
    thresh_chi2 : int
        Threshold for chi squared

    Returns
    -------
    numpy.ndarray
        Array of PSFs, each row is the PSF image at the corresponding position
        requested

    """
    if not os.path.exists(dotpsfpath):
        return FILE_NOT_FOUND

    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]

    # Check number of stars used to compute the PSF
    if PSF_model.header['ACCEPTED'] < thresh_star:
        return NOT_ENOUGH_STARS
    if PSF_model.header['CHI2'] > thresh_chi2:
        return BAD_CHI2

    PSF_basis = np.array(PSF_model.data)[0][0]
    try:
        deg = PSF_model.header['POLDEG1']
    except KeyError:
        # constant PSF model
        return PSF_basis[0, :, :]

    # scale coordinates
    x_interp, x_scale = (
        PSF_model.header['POLZERO1'],
        PSF_model.header['POLSCAL1']
    )
    y_interp, y_scale = (
        PSF_model.header['POLZERO2'],
        PSF_model.header['POLSCAL2']

    )
    xs, ys = (pos[:, 0] - x_interp) / x_scale, (pos[:, 1] - y_interp) / y_scale

    # compute polynomial coefficients
    coeffs = np.array([[x ** idx for idx in range(deg + 1)] for x in xs])
    cross_coeffs = np.array([
        np.concatenate([
            [(x ** idx_j) * (y ** idx_i) for idx_j in range(deg - idx_i + 1)]
            for idx_i in range(1, deg + 1)
        ])
        for x, y in zip(xs, ys)
    ])
    coeffs = np.hstack((coeffs, cross_coeffs))

    # compute interpolated PSF
    PSFs = np.array([
        np.sum(
            [coeff * atom for coeff, atom in zip(coeffs_posi, PSF_basis)],
            axis=0,
        )
        for coeffs_posi in coeffs
    ])

    return PSFs


class PSFExInterpolator(object):
    """The PSFEx Interpolator Class.

    This class uses a PSFEx output file to compute the PSF at desired
    positions.

    Parameters
    ----------
    dotpsf_path : str
        Path to PSFEx output file
    galcat_path : str
        Path to SExtractor-like galaxy catalogue
    output_path : str
        Path to folder where output PSFs should be written
    img_number : str
        File number string
    w_log : logging.Logger
        Logging instance
    pos_params : list, optional
        Desired position parameters, ff provided, there should be exactly two,
        and they must also be present in the galaxy catalogue; otherwise,
        they are read directly from the ``.psf`` file.
    get_shapes : bool
        If ``True`` will compute shapes for the PSF model
    star_thresh : int
        Threshold of stars under which the PSF is not interpolated
    thresh_chi2 : int
        Threshold for chi squared

    """

    def __init__(
        self,
        dotpsf_path,
        galcat_path,
        output_path,
        img_number,
        w_log,
        pos_params=None,
        get_shapes=True,
        star_thresh=20,
        chi2_thresh=2,
    ):

        # Path to PSFEx output file
        if (
            isinstance(dotpsf_path, type(None))
            or os.path.isfile(dotpsf_path)
        ):
            self._dotpsf_path = dotpsf_path
        else:
            raise ValueError(f'Cound not find file {dotpsf_path}.')
        # Path to catalogue containing galaxy positions
        if os.path.isfile(galcat_path):
            self._galcat_path = galcat_path
        else:
            raise ValueError(f'Cound not find file {galcat_path}.')
        # Path to output file to be written
        self._output_path = output_path + '/galaxy_psf'
        # Path to output file to be written for validation
        self._output_path_validation = output_path + '/validation_psf'
        # if required, compute and save shapes
        self._compute_shape = get_shapes
        # Number of stars under which we don't interpolate the PSF
        self._star_thresh = star_thresh

        self._chi2_thresh = chi2_thresh

        # Logging
        self._w_log = w_log

        # handle provided, but empty pos_params (for use within
        # CosmoStat's ShapePipe)
        if pos_params:
            if not len(pos_params) == 2:
                raise ValueError(
                    f'{len(pos_params)} position parameters were passed on; '
                    + 'there should be exactly two.'
                )
            self._pos_params = pos_params
        else:
            self._pos_params = None
        self.gal_pos = None
        self.interp_PSFs = None

        self._img_number = img_number

    def process(self):
        """Process.

        Process the PSF interpolation single-epoch run.

        """
        if self.gal_pos is None:
            self._get_galaxy_positions()

        if self.interp_PSFs is None:
            self._interpolate()

        if (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == NOT_ENOUGH_STARS
        ):
            self._w_log.info(
                'Not enough stars to interpolate the psf in the file '
                + f'{self._dotpsf_path}.'
            )
        elif (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == BAD_CHI2
        ):
            self._w_log.info(
                f'Bad chi2 for the psf model in the file {self._dotpsf_path}.'
            )
        elif (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == FILE_NOT_FOUND
        ):
            self._w_log.info(f'Psf model file {self._dotpsf_path} not found.')
        else:
            if self._compute_shape:
                self._get_psfshapes()

            self._write_output()

    def _get_position_parameters(self):
        """Get Position Parameters.

        Read position parameters from ``.psf`` file.

        """
        dotpsf = file_io.FITSCatalogue(self._dotpsf_path)
        dotpsf.open()
        self._pos_params = [
            dotpsf.get_header()['POLNAME1'],
            dotpsf.get_header()['POLNAME2']
        ]
        dotpsf.close()

    def _get_galaxy_positions(self):
        """Get Galaxy Positions.

        Extract galaxy positions from galaxy catalogue.

        """
        if self._pos_params is None:
            self._get_position_parameters()

        galcat = file_io.FITSCatalogue(self._galcat_path, SEx_catalogue=True)
        galcat.open()

        try:
            self.gal_pos = np.array([
                [x, y] for x, y in zip(
                    galcat.get_data()[self._pos_params[0]],
                    galcat.get_data()[self._pos_params[1]]
                )
            ])
            self._w_log.info(
                f'Read {self.gal_pos.shape[0]} positions from galaxy catalog'
            )

        except KeyError as detail:
            # extract erroneous position parameter from original exception
            err_pos_param = detail.args[0][4:-15]
            pos_param_err = (
                f'Required position parameter {err_pos_param}'
                + 'was not found in galaxy catalog. Leave '
                + 'pos_params (or EXTRA_CODE_OPTION) blank to '
                + 'read them from .psf file.'
            )
            raise KeyError(pos_param_err)
        galcat.close()

    def _interpolate(self):
        """Interpolate.

        Run Sheldon and Rykoff's PSFEx interpolator method at the desired
        positions.

        """
        self.interp_PSFs = interpsfex(
            self._dotpsf_path,
            self.gal_pos,
            self._star_thresh,
            self._chi2_thresh,
        )

    def _get_psfshapes(self):
        """Get PSF Shapes.

        Compute shapes of PSF at galaxy positions using HSM.

        """
        if import_fail:
            raise ImportError('Galsim is required to get shapes information')

        psf_moms = [
            hsm.FindAdaptiveMom(Image(psf), strict=False)
            for psf in self.interp_PSFs
        ]

        self.psf_shapes = np.array([
            [
                moms.observed_shape.g1,
                moms.observed_shape.g2,
                moms.moments_sigma,
                int(bool(moms.error_message))
            ]
            for moms in psf_moms
        ])

    def _write_output(self):
        """Write Output.

        Save computed PSFs to a FITS file.

        """
        output = file_io.FITSCatalogue(
            self._output_path + self._img_number + '.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True,
        )

        if self._compute_shape:
            data = {
                'VIGNET': self.interp_PSFs,
                'E1_PSF_HSM': self.psf_shapes[:, 0],
                'E2_PSF_HSM': self.psf_shapes[:, 1],
                'SIGMA_PSF_HSM': self.psf_shapes[:, 2],
                'FLAG_PSF_HSM': self.psf_shapes[:, 3].astype(int)
            }
        else:
            data = {'VIGNET': self.interp_PSFs}
        output.save_as_fits(data, sex_cat_path=self._galcat_path)

    def process_validation(self, psfex_cat_path):
        """Process Validation.

        Process validation steps.

        Parameters
        ----------
        str
            Path to PSFEx catalogue

        """
        if not os.path.isfile(psfex_cat_path):
            raise ValueError(f'Cound not find file {psfex_cat_path}.')

        if self.gal_pos is None:
            self._get_galaxy_positions()

        if self.interp_PSFs is None:
            self._interpolate()

        if (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == NOT_ENOUGH_STARS
        ):
            self._w_log.info(
                'Not enough stars to interpolate the psf in the file '
                + f'{self._dotpsf_path}.'
            )
        elif (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == BAD_CHI2
        ):
            self._w_log.info(
                f'Bad chi2 for the psf model in the file {self._dotpsf_path}.'
            )
        elif (
            isinstance(self.interp_PSFs, str)
            and self.interp_PSFs == FILE_NOT_FOUND
        ):
            self._w_log.info(f'Psf model file {self._dotpsf_path} not found.')
        else:
            star_cat = file_io.FITSCatalogue(
                self._galcat_path,
                SEx_catalogue=True,
            )
            star_cat.open()
            star_dict = {}
            star_vign = np.copy(star_cat.get_data()['VIGNET'])
            star_dict['NUMBER'] = np.copy(star_cat.get_data()['NUMBER'])
            star_dict['X'] = np.copy(star_cat.get_data()['XWIN_IMAGE'])
            star_dict['Y'] = np.copy(star_cat.get_data()['YWIN_IMAGE'])
            star_dict['RA'] = np.copy(star_cat.get_data()['XWIN_WORLD'])
            star_dict['DEC'] = np.copy(star_cat.get_data()['YWIN_WORLD'])
            star_dict['MAG'] = np.copy(star_cat.get_data()['MAG_AUTO'])
            star_dict['SNR'] = np.copy(star_cat.get_data()['SNR_WIN'])
            star_cat.close()

            self._get_psfshapes()
            self._get_starshapes(star_vign)
            psfex_cat_dict = self._get_psfexcatdict(psfex_cat_path)

            self._write_output_validation(star_dict, psfex_cat_dict)

    def _get_starshapes(self, star_vign):
        """Get Star Shapes.

        Compute shapes of stars at stars positions using HSM.

        Parameters
        ----------
        numpy.ndarray
            Array containing the star's vignets.

        """
        if import_fail:
            raise ImportError('Galsim is required to get shapes information')

        masks = np.zeros_like(star_vign)
        masks[np.where(star_vign == -1e30)] = 1

        star_moms = [
            hsm.FindAdaptiveMom(Image(star), badpix=Image(mask), strict=False)
            for star, mask in zip(star_vign, masks)
        ]

        self.star_shapes = np.array([
            [
                moms.observed_shape.g1,
                moms.observed_shape.g2,
                moms.moments_sigma,
                int(bool(moms.error_message))
            ]
            for moms in star_moms
        ])

    def _get_psfexcatdict(self, psfex_cat_path):
        """Get PSFEx Catalogue Dictionary.

        Get data from PSFEx ``.cat`` file.

        Parameters
        ----------
        psfex_cat_path : str
            Path to the ``.cat`` file from PSFEx

        Returns
        -------
        dict
            Dictionary containing information from the PFSEx ``.cat`` file

        """
        psfex_cat = file_io.FITSCatalogue(psfex_cat_path, SEx_catalogue=True)
        psfex_cat.open()

        psfex_cat_dict = {}
        psfex_cat_dict['SOURCE_NUMBER'] = np.copy(
            psfex_cat.get_data()['SOURCE_NUMBER']
        )
        psfex_cat_dict['DELTAX_IMAGE'] = np.copy(
            psfex_cat.get_data()['DELTAX_IMAGE']
        )
        psfex_cat_dict['DELTAY_IMAGE'] = np.copy(
            psfex_cat.get_data()['DELTAY_IMAGE']
        )
        psfex_cat_dict['CHI2_PSF'] = np.copy(
            psfex_cat.get_data()['CHI2_PSF']
        )

        return psfex_cat_dict

    def _write_output_validation(self, star_dict, psfex_cat_dict):
        """Write Output Validation.

        Save computed PSFs and stars to fits file.

        Parameters
        ----------
        star_dict : dict
            Dictionary containing star information
        psfex_cat_dict : dict
            Dictionary containing information from the PFSEx ``.cat`` file

        """
        output = file_io.FITSCatalogue(
            self._output_path_validation + self._img_number + '.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True,
        )

        data = {
            'E1_PSF_HSM': self.psf_shapes[:, 0],
            'E2_PSF_HSM': self.psf_shapes[:, 1],
            'SIGMA_PSF_HSM': self.psf_shapes[:, 2],
            'FLAG_PSF_HSM': self.psf_shapes[:, 3].astype(int),
            'E1_STAR_HSM': self.star_shapes[:, 0],
            'E2_STAR_HSM': self.star_shapes[:, 1],
            'SIGMA_STAR_HSM': self.star_shapes[:, 2],
            'FLAG_STAR_HSM': self.star_shapes[:, 3].astype(int)
        }
        data = {**data, **star_dict}

        data['ACCEPTED'] = np.ones_like(data['NUMBER'], dtype='int16')
        star_used = psfex_cat_dict.pop('SOURCE_NUMBER')

        for idx in range(len(data['NUMBER'])):
            if idx + 1 not in star_used:
                data['ACCEPTED'][idx] = 0

        output.save_as_fits(data, sex_cat_path=self._galcat_path)

    def process_me(self, dot_psf_dir, dot_psf_pattern, f_wcs_path):
        """Process Multi-Epoch.

        Process the multi-epoch.

        Parameters
        ----------
        dot_psf_dir : str
            Path to the directory containing the ``.psf`` files
        dot_psf_pattern : str
            Common pattern of the ``.psf`` files
        f_wcs_path : str
            Path to the log file containing the WCS for each CCDs

        """
        if os.path.exists(dot_psf_dir):
            self._dot_psf_dir = dot_psf_dir
        else:
            raise ValueError(f'Cound not find directory {dot_psf_dir}.')

        self._dot_psf_pattern = dot_psf_pattern

        if os.path.isfile(f_wcs_path):
            self._f_wcs_file = SqliteDict(f_wcs_path)
        else:
            raise ValueError(f'Cound not find file {f_wcs_path}.')

        if self.gal_pos is None:
            self._get_galaxy_positions()

        output_dict = self._interpolate_me()

        self._write_output_me(output_dict)

    def _interpolate_me(self):
        """Interpolate Multi-Epoch.

        Interpolate PSFs for multi-epoch run.

        Raises
        ------
        KeyError
            If 'N_EPOCH' key not in input catalogue

        Returns
        -------
        dict
            Dictionnary containing object Ids, the interpolated PSFs and shapes
            (optionally)

        """
        cat = file_io.FITSCatalogue(self._galcat_path, SEx_catalogue=True)
        cat.open()

        all_id = np.copy(cat.get_data()['NUMBER'])
        key_ne = 'N_EPOCH'
        if key_ne not in cat.get_data():
            raise KeyError(
                f'Key {key_ne} not found in input galaxy catalogue, needed for'
                + ' PSF interpolation to multi-epoch data; run previous module'
                + ' (SExtractor) in multi-epoch mode'
            )
        n_epoch = np.copy(cat.get_data()[key_me])

        list_ext_name = cat.get_ext_name()
        hdu_ind = [
            idx for idx in range(len(list_ext_name))
            if 'EPOCH' in list_ext_name[idx]
        ]

        final_list = []
        for hdu_index in hdu_ind:
            exp_name = cat.get_data(hdu_index)['EXP_NAME'][0]
            ccd_list = list(set(cat.get_data(hdu_index)['CCD_N']))
            array_psf = None
            array_id = None
            array_shape = None
            array_exp_name = None
            for ccd in ccd_list:
                if ccd == -1:
                    continue
                dot_psf_path = (
                    f'{self._dot_psf_dir}/{self._dot_psf_pattern}-{exp_name}'
                    + f'-{ccd}.psf'
                )
                ind_obj = np.where(cat.get_data(hdu_index)['CCD_N'] == ccd)[0]
                obj_id = all_id[ind_obj]
                gal_pos = np.array(
                    self._f_wcs_file[exp_name][ccd]['WCS'].all_world2pix(
                        self.gal_pos[:, 0][ind_obj],
                        self.gal_pos[:, 1][ind_obj],
                        0,
                    )
                ).T

                self.interp_PSFs = interpsfex(
                    dot_psf_path,
                    gal_pos,
                    self._star_thresh,
                    self._chi2_thresh,
                )

                if (
                    isinstance(self.interp_PSFs, str)
                    and self.interp_PSFs == NOT_ENOUGH_STARS
                ):
                    self._w_log.info(
                        f'Not enough stars find in the ccd {ccd} of the '
                        + f'exposure {exp_name}. Object inside this ccd will '
                        + 'lose an epoch.'
                    )
                    continue
                if (
                    isinstance(self.interp_PSFs, str)
                    and self.interp_PSFs == BAD_CHI2
                ):
                    self._w_log.info(
                        f'Bad chi2 for the psf model in the ccd {ccd} of the '
                        + f'exposure {exp_name}. Object inside this ccd will '
                        + 'lose an epoch.'
                    )
                    continue
                if (
                    isinstance(self.interp_PSFs, str)
                    and self.interp_PSFs == FILE_NOT_FOUND
                ):
                    self._w_log.info(
                        f'Psf model file {self._dotpsf_path} not found. '
                        + 'Object inside this ccd will lose an epoch.'
                    )
                    continue

                if array_psf is None:
                    array_psf = np.copy(self.interp_PSFs)
                else:
                    array_psf = np.concatenate(
                        (array_psf, np.copy(self.interp_PSFs))
                    )

                if array_id is None:
                    array_id = np.copy(obj_id)
                else:
                    array_id = np.concatenate((array_id, np.copy(obj_id)))

                if self._compute_shape:
                    self._get_psfshapes()
                    if array_shape is None:
                        array_shape = np.copy(self.psf_shapes)
                    else:
                        array_shape = np.concatenate((
                            array_shape,
                            np.copy(self.psf_shapes),
                        ))
                else:
                    array_shape = None

                exp_name_tmp = np.array([
                    exp_name + '-' + str(ccd) for _ in range(len(obj_id))
                ])
                if array_exp_name is None:
                    array_exp_name = exp_name_tmp
                else:
                    array_exp_name = np.concatenate(
                        (array_exp_name, exp_name_tmp)
                    )

            final_list.append([
                array_id,
                array_psf,
                array_shape,
                array_exp_name
            ])

        self._f_wcs_file.close()
        cat.close()

        output_dict = {}
        n_empty = 0
        for id_tmp in all_id:
            output_dict[id_tmp] = {}
            counter = 0
            for j in range(len(final_list)):
                where_res = np.where(final_list[j][0] == id_tmp)[0]
                if (len(where_res) != 0):
                    output_dict[id_tmp][final_list[j][3][where_res[0]]] = {}
                    output_dict[id_tmp][
                        final_list[j][3][where_res[0]]
                    ]['VIGNET'] = final_list[j][1][where_res[0]]
                    if self._compute_shape:
                        shape_dict = {}
                        shape_dict['E1_PSF_HSM'] = (
                            final_list[j][2][where_res[0]][0]
                        )
                        shape_dict['E2_PSF_HSM'] = (
                            final_list[j][2][where_res[0]][1]
                        )
                        shape_dict['SIGMA_PSF_HSM'] = (
                            final_list[j][2][where_res[0]][2]
                        )
                        shape_dict['FLAG_PSF_HSM'] = (
                            final_list[j][2][where_res[0]][3]
                        )
                        output_dict[id_tmp][
                            final_list[j][3][where_res[0]]
                        ]['SHAPES'] = shape_dict
                    counter += 1
            if counter == 0:
                output_dict[id_tmp] = 'empty'
                n_empty += 1

        self._w_log.info(f'{n_empty}/{len(all_id)} PSFs are empty')

        return output_dict

    def _write_output_me(self, output_dict):
        """Write Output Multi-Epoch.

        Save computed PSFs to numpy object file for multi-epoch run.

        Parameters
        ----------
        output_dict : dict
            Dictionnary of outputs to save

        """
        output_file = SqliteDict(
            self._output_path + self._img_number + '.sqlite'
        )
        for idx in output_dict.keys():
            output_file[str(idx)] = output_dict[idx]
        output_file.commit()
        output_file.close()
