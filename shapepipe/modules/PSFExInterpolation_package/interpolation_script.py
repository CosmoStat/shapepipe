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


def interpsfex(dotpsfpath, pos):
    """Use PSFEx generated model to perform spatial PSF interpolation.

        Parameters
        ----------
        dotpsfpath : string
            Path to .psf file (PSFEx output).

        pos : np.ndaray
            Positions where the PSF model should be evaluated.

        Returns
        -------
        PSFs : np.ndarray
            Each row is the PSF imagette at the corresponding asked position.
    """

    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]
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

    """

    def __init__(self, dotpsf_path, galcat_path, output_path, img_number,
                 pos_params=None, get_shapes=True):

        # Path to PSFEx output file
        self._dotpsf_path = dotpsf_path
        # Path to catalog containing galaxy positions
        self._galcat_path = galcat_path
        # Path to output file to be written
        self._output_path = output_path+'/galaxy_psf'

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

        # get number naming convention for this particular run
        #s = re.split(r"\-([0-9]*)\-([0-9]+)\.", self._galcat_path)
        #self._img_number = '-{0}-{1}'.format(s[1], s[2])
        self._img_number = img_number

        # if required, compute and save shapes
        if get_shapes:
            self._get_psfshapes()
            self._has_shapes = True
        else:
            self._has_shapes = False

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
        if self.gal_pos is None:
            self._get_galaxy_positions()

        # pex = psfex.PSFEx(self._dotpsf_path)
        # self.interp_PSFs = np.array([pex.get_rec(x,y) for x,y in
        # zip(self.gal_pos[:,0],
        #                              self.gal_pos[:,1])])

        self.interp_PSFs = interpsfex(self._dotpsf_path, self.gal_pos)

    def _get_psfshapes(self):
        """ Compute shapes of PSF at galaxy positions using HSM.

        """
        if import_fail:
            raise ImportError('Galsim is required to get shapes information')

        if self.interp_PSFs is None:
            self._interpolate()
        psf_moms = [hsm.FindAdaptiveMom(Image(psf), strict=False)
                    for psf in self.interp_PSFs]

        self.psf_shapes = np.array([[moms.observed_shape.g1,
                                     moms.observed_shape.g2,
                                     moms.moments_sigma] for moms in psf_moms])
        self.hsm_flags = np.array([bool(mom.error_message)
                                   for mom in psf_moms]).astype(int)

    def write_output(self):
        """ Save computed PSFs to fits file.

        """
        if self.interp_PSFs is None:
            self._interpolate()
        output = sc.FITSCatalog(self._output_path+self._img_number+'.fits',
                                open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                SEx_catalog=True)
        if self._has_shapes:
            data = {'VIGNET': self.interp_PSFs,
                    'E1_PSF_HSM': self.psf_shapes[:, 0],
                    'E2_PSF_HSM': self.psf_shapes[:, 1],
                    'SIGMA_PSF_HSM': self.psf_shapes[:, 2],
                    'HSM_FLAG': self.hsm_flags}
        else:
            data = {'VIGNET': self.interp_PSFs}
        output.save_as_fits(data, sex_cat_path=self._galcat_path)
