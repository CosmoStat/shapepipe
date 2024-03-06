"""MCCD INTERPOLATION SCRIPT.

This module computes the PSFs from a MCCD model at several galaxy positions.

:Author: Tobias Liaudat

"""

import os

import mccd
import numpy as np
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


def interp_MCCD(mccd_model_path, positions, ccd):
    r"""Interpolate MCCD.

    Interpolate MCCD model on requested positions.

    Parameters
    ----------
    mccd_model_path: str
        Path to the correct MCCD model
    positions: numpy.ndarray
        Positions for the PSF recovery in local coordinates with respect to
        the ``ccd``; array shape: ``(n_pos, 2)``
    ccd: int
        CCD from where the positions were taken

    Returns
    -------
    PSFs: numpy.ndarray or str
        Recovered PSFs at the required positions; the returned array has the
        shape: ``(n_pos, n_pix, n_pix)``; a ``str`` is returned if some error
        occurs

    """
    if not os.path.exists(mccd_model_path):
        return FILE_NOT_FOUND

    # Import the model
    mccd_instance = mccd.mccd_quickload(mccd_model_path)

    # Create ccd number array
    ccd_list = (np.ones((positions.shape[0])) * ccd).astype(int)

    # Positions to global coordinates
    loc2glob = mccd.mccd_utils.Loc2Glob()
    glob_pos = np.array([
        loc2glob.loc2glob_img_coord(_ccd, _pos[0], _pos[1])
        for _ccd, _pos in zip(ccd_list, positions)
    ])

    # Interpolate the model
    PSFs = mccd_instance.interpolate_psf_pipeline(glob_pos, ccd)

    del mccd_instance

    # MCCD model returns None if there were not enough stars to train the model
    # in the requested ccd.
    if PSFs is None:
        return NOT_ENOUGH_STARS

    return PSFs


class MCCDinterpolator(object):
    """The MCCD Interpolator Class.

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
    img_number: str
        Corresponds to shapepipe's module input ``file_number_string`` which
        is the input file ID
    w_log : logging.Logger
        Logging instance
    pos_params : list, optional
        Desired position parameters, if provided, there should be exactly two,
        and they must also be present in the galaxy catalogue; otherwise,
        they are read directly from the ``.psf`` file.
    get_shapes : bool
        If ``True`` will compute shapes for the PSF model

    """

    def __init__(
        self,
        dotpsf_path,
        galcat_path,
        output_path,
        img_number,
        w_log,
        pos_params=None,
        get_shapes=True
    ):
        # Path to PSFEx output file
        self._dotpsf_path = dotpsf_path
        # Path to catalogue containing galaxy positions
        self._galcat_path = galcat_path
        # Path to output file to be written
        self._output_path = output_path + '/galaxy_psf'
        # Path to output file to be written for validation
        self._output_path_validation = output_path + '/validation_psf'
        # if required, compute and save shapes
        self._compute_shape = get_shapes

        # Paths to the PSF model
        self._dot_psf_dir = None
        self._dot_psf_pattern = None
        self._f_wcs_file = None

        # Logging
        self._w_log = w_log

        # handle provided, but empty pos_params (for use within
        # CosmoStat's ShapePipe)
        if pos_params:
            if not len(pos_params) == 2:
                raise ValueError(
                    f'{len(pos_params)} position parameters were passed on;'
                    + f' there should be exactly two.'
                )
            self._pos_params = pos_params
        else:
            self._pos_params = None
        self.gal_pos = None
        self.interp_PSFs = None

        self._img_number = img_number

    def _get_position_parameters(self):
        """Get Position Parameters.

        Read position parameters from a ``.psf`` file.

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

        Extract galaxy positions from a galaxy catalogue.

        """
        if self._pos_params is None:
            self._get_position_parameters()

        galcat = file_io.FITSCatalogue(self._galcat_path, SEx_catalogue=True)
        galcat.open()
        try:
            self.gal_pos = np.array(
                [[x, y] for x, y in zip(
                    galcat.get_data()[self._pos_params[0]],
                    galcat.get_data()[self._pos_params[1]]
                )]
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

        Save computed PSFs to fits file.

        """
        output = file_io.FITSCatalogue(
            self._output_path + self._img_number + '.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True
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

    def _get_starshapes(self, star_vign):
        """Get Star Shapes.

        Compute shapes of stars at stars positions using HSM.

        Parameters
        ----------
        star_vign : numpy.ndarray
            Array containing the star's vignets

        """
        if import_fail:
            raise ImportError('Galsim is required to get shapes information')

        masks = np.zeros_like(star_vign)
        masks[np.where(star_vign == -1e30)] = 1

        star_moms = [hsm.FindAdaptiveMom(
            Image(star),
            badpix=Image(mask),
            strict=False)
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

        Get data from a PSFEx ``.cat`` file.

        Parameters
        ----------
        psfex_cat_path : str
            Path to the ``.cat`` file from PSFEx

        Returns
        -------
        psfex_cat_dict : dict
            Dictionary containing information from PFSEx ``.cat`` file

        """
        psfex_cat = file_io.FITSCatalogue(psfex_cat_path, SEx_catalogue=True)
        psfex_cat.open()

        psfex_cat_dict = {}
        psfex_cat_dict['SOURCE_NUMBER'] = np.copy(
            psfex_cat.get_data()['SOURCE_NUMBER'])
        psfex_cat_dict['DELTAX_IMAGE'] = np.copy(
            psfex_cat.get_data()['DELTAX_IMAGE'])
        psfex_cat_dict['DELTAY_IMAGE'] = np.copy(
            psfex_cat.get_data()['DELTAY_IMAGE'])
        psfex_cat_dict['CHI2_PSF'] = np.copy(psfex_cat.get_data()['CHI2_PSF'])

        return psfex_cat_dict

    def _write_output_validation(self, star_dict, psfex_cat_dict):
        """Write Output Validation.

        Save computed PSFs and stars to a FITS file.

        Parameters
        ----------
        star_dict : dict
            Dictionary containing star informations
        psfex_cat_dict : dict
            Dictionary containing information from PFSEx ``.cat`` file

        """
        output = file_io.FITSCatalogue(
            self._output_path_validation + self._img_number + '.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True)

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

        for i in range(len(data['NUMBER'])):
            if i + 1 not in star_used:
                data['ACCEPTED'][i] = 0

        output.save_as_fits(data, sex_cat_path=self._galcat_path)

    def process_me(self, dot_psf_dir, dot_psf_pattern, f_wcs_path):
        """Process Multi-Epoch.

        Parameters
        ----------
        dot_psf_dir : str
            Path to the directory containing the PSF model files.
        dot_psf_pattern : str
            Common pattern of the PSF model files; e.g. ``fitted_model``
        f_wcs_path : str
            Path to the log file containing the WCS for each CCD

        """
        self._dot_psf_dir = dot_psf_dir
        self._dot_psf_pattern = dot_psf_pattern
        self._f_wcs_file = SqliteDict(f_wcs_path)

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
        output_dict: dict
            Dictionnary containing object IDs, the interpolated PSFs and
            shapes (optionally)

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
        n_epoch = np.copy(cat.get_data()[key_ne])

        list_ext_name = cat.get_ext_name()
        hdu_ind = [i for i in range(len(list_ext_name)) if
                   'EPOCH' in list_ext_name[i]]

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
                # dot_psf_path = self._dot_psf_dir + '/' +\
                # self._dot_psf_pattern + '-' + exp_name + '-' + str(ccd) +\
                # '.psf'
                mccd_model_path = self._dot_psf_dir + '/' +\
                    self._dot_psf_pattern + '-' + exp_name + '.npy'

                ind_obj = np.where(cat.get_data(hdu_index)['CCD_N'] == ccd)[0]
                obj_id = all_id[ind_obj]
                gal_pos = np.array(
                    self._f_wcs_file[exp_name][ccd]['WCS'].all_world2pix(
                        self.gal_pos[:, 0][ind_obj],
                        self.gal_pos[:, 1][ind_obj],
                        0
                    )
                ).T

                self.interp_PSFs = interp_MCCD(mccd_model_path, gal_pos, ccd)

                if (
                    isinstance(self.interp_PSFs, str)
                    and (self.interp_PSFs == NOT_ENOUGH_STARS)
                ):
                    self._w_log.info(
                        f'Not enough stars find in the ccd {ccd} of the'
                        + f' exposure {exp_name}. Object inside this '
                        + f'ccd will lose an epoch.'
                    )
                    continue

                if (
                    isinstance(self.interp_PSFs, str)
                    and (self.interp_PSFs == BAD_CHI2)
                ):
                    self._w_log.info(
                        f'Bad chi2 for the psf model in the ccd {ccd} of the'
                        + f' exposure {exp_name}. Object inside this ccd'
                        + f' will lose an epoch.'
                    )

                    continue

                if (
                    isinstance(self.interp_PSFs, str)
                    and (self.interp_PSFs == FILE_NOT_FOUND)
                ):
                    self._w_log.info(
                        f'Psf model file {self._dotpsf_path} not found.'
                        + f' Object inside this ccd will lose an epoch.'
                    )
                    continue

                if array_psf is None:
                    array_psf = np.copy(self.interp_PSFs)
                else:
                    array_psf = np.concatenate(
                        (array_psf, np.copy(self.interp_PSFs)))

                if array_id is None:
                    array_id = np.copy(obj_id)
                else:
                    array_id = np.concatenate((array_id, np.copy(obj_id)))

                if self._compute_shape:
                    self._get_psfshapes()
                    if array_shape is None:
                        array_shape = np.copy(self.psf_shapes)
                    else:
                        array_shape = np.concatenate(
                            (array_shape, np.copy(self.psf_shapes)))
                else:
                    array_shape = None

                exp_name_tmp = np.array(
                    [exp_name + '-' + str(ccd) for i in range(len(obj_id))])
                if array_exp_name is None:
                    array_exp_name = exp_name_tmp
                else:
                    array_exp_name = np.concatenate(
                        (array_exp_name, exp_name_tmp))

            final_list.append(
                [array_id, array_psf, array_shape, array_exp_name])

        self._f_wcs_file.close()
        cat.close()

        output_dict = {}
        for id_tmp in all_id:
            output_dict[id_tmp] = {}
            counter = 0
            for j in range(len(final_list)):
                where_res = np.where(final_list[j][0] == id_tmp)[0]
                if len(where_res) != 0:
                    output_dict[id_tmp][final_list[j][3][where_res[0]]] = {}
                    output_dict[id_tmp][final_list[j][3][where_res[0]]][
                        'VIGNET'] = final_list[j][1][where_res[0]]
                    if self._compute_shape:
                        shape_dict = {}
                        shape_dict['E1_PSF_HSM'] = \
                            final_list[j][2][where_res[0]][0]
                        shape_dict['E2_PSF_HSM'] = \
                            final_list[j][2][where_res[0]][1]
                        shape_dict['SIGMA_PSF_HSM'] = \
                            final_list[j][2][where_res[0]][2]
                        shape_dict['FLAG_PSF_HSM'] = \
                            final_list[j][2][where_res[0]][3]
                        output_dict[id_tmp][final_list[j][3][where_res[0]]][
                            'SHAPES'] = shape_dict
                    counter += 1
            if counter == 0:
                output_dict[id_tmp] = 'empty'

        return output_dict

    def _write_output_me(self, output_dict):
        """Write Output Multi-Epoch.

        Save computed PSFs to numpy object file for multi-epoch run.

        Parameters
        ----------
        output_dict : dict
            Dictionnary of outputs to save

        """
        # np.save(self._output_path+self._img_number, output_dict)

        output_file = SqliteDict(
            self._output_path + self._img_number + '.sqlite')
        for i in output_dict.keys():
            output_file[str(i)] = output_dict[i]
        output_file.commit()
        output_file.close()
