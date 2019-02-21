# -*- coding: utf-8 -*-

"""INTERPOLATION SCRIPT

This module computes the PSFs from a PSFEx model at several galaxy positions.

Was using Erin Sheldon & Eli Rykoff's psfex module, available on GitHub at:
https://github.com/esheldon/psfex

:Author: Morgan Schmitz

:Version: 1.2.0

:Date: 06/02/2018

"""

import numpy as np
# import psfex
import re
from astropy.io import fits

try:
    import galsim.hsm as hsm
    from galsim import Image
except ImportError:
    import_fail = True
else:
    import_fail = False

from shapepipe.pipeline import file_io as sc


NOT_ENOUGHT_STARS = 'Fail'


def interpsfex(dotpsfpath, pos, thresh):
    """Use PSFEx generated model to perform spatial PSF interpolation.

        Parameters
        ----------
        dotpsfpath : string
            Path to .psf file (PSFEx output).

        pos : np.ndaray
            Positions where the PSF model should be evaluated.

        thresh : int
            Threshold of stars under which the PSF is not interpolated

        Returns
        -------
        PSFs : np.ndarray
            Each row is the PSF imagette at the corresponding asked position.

    """

    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]

    if PSF_model.header['ACCEPTED'] < thresh:
        return NOT_ENOUGHT_STARS

    PSF_basis = np.array(PSF_model.data)[0][0]
    try:
        deg = PSF_model.header['POLDEG1']
    except KeyError:
        # constant PSF model
        return PSF_basis[0, :, :]

    # scale coordinates
    x_interp, x_scale = (PSF_model.header['POLZERO1'],
                         PSF_model.header['POLSCAL1'])
    y_interp, y_scale = (PSF_model.header['POLZERO2'],
                         PSF_model.header['POLSCAL2'])
    xs, ys = (pos[:, 0] - x_interp) / x_scale, (pos[:, 1] - y_interp) / y_scale

    # compute polynomial coefficients
    coeffs = np.array([[x**i for i in range(deg+1)] for x in xs])
    cross_coeffs = np.array([np.concatenate([[(x ** j) * (y ** i)
                                              for j in range(deg - i + 1)]
                                             for i in range(1, deg + 1)])
                             for x, y in zip(xs, ys)])
    coeffs = np.hstack((coeffs, cross_coeffs))

    # compute interpolated PSF
    PSFs = np.array([np.sum([coeff * atom for coeff, atom in
                     zip(coeffs_posi, PSF_basis)], axis=0)
                     for coeffs_posi in coeffs])

    return PSFs


class PSFExInterpolator(object):
    """Interpolator class.

    This class uses a PSFEx output file to compute the PSF at desired
    positions.

    Parameters
    ----------
    dotpsf_path : str
        Path to PSFEx output file.
    galcat_path : str
        Path to SExtractor-like galaxy catalog.
    output_path : str
        Path to folder where output PSFs should be written.
    pos_params : list, optional
        Desired position parameters. If provided, there should be exactly two,
        and they must also be present in the galaxy catalog. Otherwise,
        they are read directly from the .psf file.
    get_shapes : boolean
        If True will compute shapes for the PSF model.
    star_thresh : int
        Threshold of stars under which the PSF is not interpolated.

    """

    def __init__(self, dotpsf_path, galcat_path, output_path, img_number,
                 pos_params=None, get_shapes=True, star_thresh=20):

        # Path to PSFEx output file
        self._dotpsf_path = dotpsf_path
        # Path to catalog containing galaxy positions
        self._galcat_path = galcat_path
        # Path to output file to be written
        self._output_path = output_path+'/galaxy_psf'
        # if required, compute and save shapes
        self._make_shape = get_shapes
        # Number of stars under which we don't interpolate the PSF
        self._star_thresh = star_thresh

        # handle provided, but empty pos_params (for use within
        # CosmoStat's ShapePipe)
        if pos_params:
            if not len(pos_params) == 2:
                raise ValueError('{} position parameters were passed on; '
                                 'there should be exactly two.'
                                 ''.format(len(pos_params)))
            self._pos_params = pos_params
        else:
            self._pos_params = None
        self.gal_pos = None
        self.interp_PSFs = None

        self._img_number = img_number

        # # if required, compute and save shapes
        # if get_shapes:
        #     self._get_psfshapes()
        #     self._has_shapes = True
        # else:
        #     self._has_shapes = False

    def process(self):
        """ Process the PSF interpolation single-epoch run.

        """
        if self.gal_pos is None:
            self._get_galaxy_positions()

        if self.interp_PSFs is None:
            self._interpolate()

        if self.interp_PSFs != NOT_ENOUGHT_STARS:
            if self._make_shape:
                self._get_psfshapes()

            self._write_output()

    def _get_position_parameters(self):
        """ Read position parameters from .psf file.

        """

        dotpsf = sc.FITSCatalog(self._dotpsf_path)
        dotpsf.open()
        self._pos_params = [dotpsf.get_header()['POLNAME1'],
                            dotpsf.get_header()['POLNAME2']]
        dotpsf.close()

    def _get_galaxy_positions(self):
        """ Extract galaxy positions from galaxy catalog.

        """
        if self._pos_params is None:
            self._get_position_parameters()

        galcat = sc.FITSCatalog(self._galcat_path, SEx_catalog=True)
        galcat.open()
        try:
            self.gal_pos = np.array([[x, y] for x, y in
                                    zip(galcat.get_data()[self._pos_params[0]],
                                    galcat.get_data()[self._pos_params[1]])])
        except KeyError as detail:
            # extract erroneous position parameter from original exception
            err_pos_param = detail.args[0][4:-15]
            pos_param_err = ('Required position parameter ' + err_pos_param +
                             'was not found in galaxy catalog. Leave '
                             'pos_params (or EXTRA_CODE_OPTION) blank to '
                             'read them from .psf file.')
            raise KeyError(pos_param_err)
        galcat.close()

    def _interpolate(self):
        """

        (Run Sheldon & Rykoff's PSFEx interpolator method at desired
        positions.)

        """
        # if self.gal_pos is None:
        #     self._get_galaxy_positions()

        # pex = psfex.PSFEx(self._dotpsf_path)
        # self.interp_PSFs = np.array([pex.get_rec(x,y) for x,y in
        # zip(self.gal_pos[:,0],
        #                              self.gal_pos[:,1])])

        # self.interp_PSFs = interpsfex(self._dotpsf_path, self.gal_pos)
        self.interp_PSFs = interpsfex(self._dotpsf_path, self.gal_pos, self._star_thresh)

    def _get_psfshapes(self):
        """ Compute shapes of PSF at galaxy positions using HSM.

        """
        if import_fail:
            raise ImportError('Galsim is required to get shapes information')

        # if self.interp_PSFs is None:
        #     self._interpolate()
        psf_moms = [hsm.FindAdaptiveMom(Image(psf), strict=False)
                    for psf in self.interp_PSFs]

        self.psf_shapes = np.array([[moms.observed_shape.g1,
                                     moms.observed_shape.g2,
                                     moms.moments_sigma,
                                     int(bool(moms.error_message))] for moms in psf_moms])
        # self.hsm_flags = np.array([bool(mom.error_message)
        #                            for mom in psf_moms]).astype(int)

    # def write_output(self):
    def _write_output(self):
        """ Save computed PSFs to fits file.

        """
        # if self.interp_PSFs is None:
        #     self._interpolate()
        output = sc.FITSCatalog(self._output_path+self._img_number+'.fits',
                                open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                SEx_catalog=True)
        # if self._has_shapes:
        if self._make_shape:
            data = {'VIGNET': self.interp_PSFs,
                    'E1_PSF_HSM': self.psf_shapes[:, 0],
                    'E2_PSF_HSM': self.psf_shapes[:, 1],
                    'SIGMA_PSF_HSM': self.psf_shapes[:, 2],
                    # 'HSM_FLAG': self.hsm_flags}
                    'HSM_FLAG': self.psf_shapes[:, 3].astype(int)}
        else:
            data = {'VIGNET': self.interp_PSFs}
        output.save_as_fits(data, sex_cat_path=self._galcat_path)

    def process_me(self, dot_psf_dir, dot_psf_pattern, f_wcs_path):
        """ Process the multi-epoch.

        Parameters
        ----------
        dot_psf_dir : str
            Path to the directory containing the ".psf" files.
        dot_psf_pattern : str
            Common pattern of the ".psf" files.
        f_wcs_path : str
            Path to the log file containing the WCS for each CCDs.

        """

        self._dot_psf_dir = dot_psf_dir
        self._dot_psf_pattern = dot_psf_pattern
        self._f_wcs_file = np.load(f_wcs_path).item()

        if self.gal_pos is None:
            self._get_galaxy_positions()

        output_list = self._interpolate_me()

        self._write_output_me(output_list)

    def _interpolate_me(self):
        """ Interpolate PSFs for multi-epoch run.

        Returns
        -------
        list
            List contianing object Ids, the interpolated PSFs and shapes (optionally)

        """

        cat = sc.FITSCatalog(self._galcat_path, SEx_catalog=True)
        cat.open()

        all_id = np.copy(cat.get_data()['NUMBER'])
        n_epoch = np.copy(cat.get_data()['N_EPOCH'])

        l = cat.get_ext_name()
        hdu_ind = [i for i in range(len(l)) if 'EPOCH' in cat.get_ext_name(i)]

        final_list = []
        for i in hdu_ind:
            exp_name = cat.get_data(i)['EXP_NAME'][0]
            ccd_list = list(set(cat.get_data(i)['CCD_N']))
            array_psf = None
            array_id = None
            array_shape = None
            for j in ccd_list:
                if j == -1:
                    continue
                dot_psf_path = self._dot_psf_dir + '/' + self._dot_psf_pattern + '-' + exp_name + '-' + str(j) + '.psf'
                ind_obj = np.where(cat.get_data(i)['CCD_N']==j)[0]
                obj_id = all_id[ind_obj]
                gal_pos = np.array(self._f_wcs_file[exp_name][j].all_pix2world(self.gal_pos[:,0][ind_obj], self.gal_pos[:,1][ind_obj], 0)).T

                self.interp_PSFs = interpsfex(dot_psf_path, gal_pos, self._star_thresh)

                if self.interp_PSFs == NOT_ENOUGHT_STARS:
                    continue

                if array_psf is None:
                    array_psf = np.copy(self.interp_PSFs)
                else:
                    array_psf = np.concatenate((array_psf, np.copy(self.interp_PSFs)))

                if array_id is None:
                    array_id = np.copy(obj_id)
                else:
                    array_id = np.concatenate((array_id, np.copy(obj_id)))

                if self._make_shape:
                    self._get_psfshapes()
                    if array_shape is None:
                        array_shape = np.copy(self.psf_shapes)
                    else:
                        array_shape = np.concatenate((array_shape, np.copy(self.psf_shapes)))
                else:
                    array_shape = None

            final_list.append([array_id, array_psf, array_shape])

        cat.close()

        output_list_id = [[] for i in range(max(n_epoch))]
        output_list_vign = [[] for i in range(max(n_epoch))]
        output_list_shape = [[] for i in range(max(n_epoch))]
        for i in range(len(all_id)):
            k = 0
            for j in range(len(final_list)):
                where_res = np.where(final_list[j][0] == all_id[i])[0]

                if (len(where_res) != 0):
                    output_list_id[k].append(final_list[j][0][where_res])
                    output_list_vign[k].append(final_list[j][1][where_res])
                    if self._make_shape:
                        output_list_shape[k].append(final_list[j][2][where_res])
                    k += 1

        return [output_list_id, output_list_vign, output_list_shape]


    def _write_output_me(self, output_list):
        """ Save computed PSFs to fits file for multi-epoch run.

        Parameters
        ----------
        output_list : list
            List of outputs to save
            
        """

        output_file = sc.FITSCatalog(self._output_path+self._img_number+'.fits',
                                     open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                     SEx_catalog=True)

        for i in range(len(output_list[0])):
            out_dict = {}
            out_dict['NUMBER'] = np.array(output_list[0][i]).squeeze()
            out_dict['VIGNET'] = np.array(output_list[1][i]).squeeze()
            if self._make_shape:
                out_dict['E1_PSF_HSM'] = np.array(output_list[2][i]).squeeze()[:, 0]
                out_dict['E2_PSF_HSM'] = np.array(output_list[2][i]).squeeze()[:, 1]
                out_dict['SIGMA_PSF_HSM'] = np.array(output_list[2][i]).squeeze()[:, 2]
                out_dict['HSM_FLAG'] = np.array(output_list[2][i]).squeeze()[:, 3].astype(int)

            output_file.save_as_fits(out_dict,
                                     ext_name='N_EPOCH_{}'.format(i+1),
                                     sex_cat_path=self._galcat_path)
