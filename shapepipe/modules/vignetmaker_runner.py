# -*- coding: utf-8 -*-

"""VIGNET MAKER RUNNER

This file contains methods to create postage stamps from images.

:Author: Axel Guinot

"""

import numpy as np
from astropy.wcs import WCS
import re
from sf_tools.image.stamp import FetchStamps
import shapepipe.pipeline.file_io as io
from shapepipe.modules.module_decorator import module_runner


class vignetmaker(object):
    """ Vignetmaker

    This class handle the creation of vignets.

    Parameters
    ----------
    galcat_path : str
        Path to catalog containing the positions.
    pos_type : str
        Must be in ['PIX', 'SPHE']
        PIX : position in pixel coordinates.
        SPHE : position in world coordinates.
    pos_params : list
        Parameters to use as positions ['xpos', 'ypos'].
    output_dir : str
        Path to the output directory.
    image_num : str
        Image numbering.

    """

    def __init__(self, galcat_path, pos_type, pos_params,
                 output_dir, image_num):

        self._galcat_path = galcat_path
        self._output_dir = output_dir
        self._image_num = image_num

        self._pos = self.get_pos(pos_params)
        self._pos_type = pos_type

    def process(self, image_path_list, rad, suffix):
        """ Process

        Main function to create the stamps

        Parameters
        ----------
        image_path_list : list
            List of path for the input images.
        rad : int
            Radius to use for the stamps.
        suffix : str
            Suffix of the output file.

        """

        for _suffix, img in zip(suffix, image_path_list):
            image_path = img

            if self._pos_type == 'PIX':
                pos = self._pos
            elif self._pos_type == 'SPHE':
                pos = convert_pos(image_path)
            else:
                raise ValueError('Coordinates type must be in : PIX (pixel), '
                                 'SPHE (spherical).')

            vign = self._get_stamp(image_path, pos-1, rad)

            save_vignet(vign, self._galcat_path, self._output_dir, _suffix,
                        self._image_num)

    def get_pos(self, pos_params):
        """Get positions

        Get the positions of the given parameters from SExtractor catalog.

        Parameters
        ----------
        pos_params : list
            List of string containing the SExtractor's parameters to use as
            positions

        Returns
        -------
        numpy.ndarray
            Array of the positions

        """

        f = io.FITSCatalog(self._galcat_path, SEx_catalog=True)
        f.open()

        pos = np.array([f.get_data()[pos_params[1]],
                        f.get_data()[pos_params[0]]]).T

        f.close()

        return pos

    def convert_pos(self, image_path):
        """Convert position

        Convert positions from world coordinates to pixel coordinates.

        Parameters
        ----------
        image_path : str
            Path to the image from where the stamp are created.

        Return
        ------
        numpy.ndarray
            New positions in pixel coordinates.

        """

        f = io.FITSCatalog(image_path)
        f.open()
        h = f.get_header(0)
        f.close()

        w = WCS(h)

        pos_tmp = np.copy(self._pos)
        pos_tmp[:, [0, 1]] = pos_tmp[:, [1, 0]]

        new_pos = w.all_world2pix(pos_tmp, 1)

        new_pos[:, [0, 1]] = new_pos[:, [1, 0]]

        return new_pos

    def _get_stamp(self, img_path, pos, rad):
        """Get stamp

        Extract stamp at given positions on the image

        Parameters
        ----------
        img_path : str
            Path to the image
        pos : numpy.ndarray
            Array of positions for the vignet's centers.
        rad : int
            Radius of the stamp, must be odd.

        Returns
        -------
        numpy.array
            Array containing the vignets

        """

        img_file = io.FITSCatalog(img_path)
        img_file.open()
        img = img_file.get_data(0)
        img_file.close()

        fs = FetchStamps(img, int(rad))
        fs.get_pixels(np.round(pos).astype(int))

        vign = fs.scan()

        return vign

    def _get_stamp_me(self, image_dir, image_pattern):
        """ Get stamp ME

        Get stamps for multi-epoch data.

        Parameters
        ----------
        image_dir : str
            Path to the directory where the image are.
        image_pattern : str
            Common part of the files.

        Returns
        -------
        output_dict : dict
            Directory containing object id and vignets for each epoch.

        """

        cat = io.FITSCatalog(self._galcat_path, SEx_catalog=True)
        cat.open()

        all_id = np.copy(cat.get_data()['NUMBER'])
        n_epoch = np.copy(cat.get_data()['N_EPOCH'])

        list_ext_name = cat.get_ext_name()
        hdu_ind = [i for i in range(len(list_ext_name))
                   if 'EPOCH' in list_ext_name[i]]

        final_list = []
        for hdu_index in hdu_ind:
            exp_name = cat.get_data(hdu_index)['EXP_NAME'][0]
            ccd_list = list(set(cat.get_data(hdu_index)['CCD_N']))
            array_vign = None
            array_id = None
            array_exp_name = None
            for ccd in ccd_list:
                if ccd == -1:
                    continue
                img_path = (image_dir + '/' + image_pattern + '-' +
                            exp_name + '-' + str(ccd) + '.fits')
                ind_obj = np.where(cat.get_data(hdu_index)['CCD_N'] == ccd)[0]
                obj_id = all_id[ind_obj]
                pos = np.array(self._f_wcs_file[exp_name][ccd]['WCS'].all_world2pix(self._pos[:, 1][ind_obj], self._pos[:, 0][ind_obj], 1)).T
                pos[:, [0, 1]] = pos[:, [1, 0]]

                tmp_vign = self._get_stamp(img_path, pos-1, self._rad)

                if array_vign is None:
                    array_vign = np.copy(tmp_vign)
                else:
                    array_vign = np.concatenate((array_vign,
                                                 np.copy(tmp_vign)))

                if array_id is None:
                    array_id = np.copy(obj_id)
                else:
                    array_id = np.concatenate((array_id, np.copy(obj_id)))

                exp_name_tmp = np.array([exp_name + '-' + str(ccd)
                                         for i in range(len(obj_id))])
                if array_exp_name is None:
                    array_exp_name = exp_name_tmp
                else:
                    array_exp_name = np.concatenate((array_exp_name,
                                                     exp_name_tmp))

            final_list.append([array_id, array_vign, array_exp_name])

        cat.close()

        output_dict = {}
        for id_tmp in all_id:
            output_dict[id_tmp] = {}
            counter = 0
            for j in range(len(final_list)):
                where_res = np.where(final_list[j][0] == id_tmp)[0]
                if (len(where_res) != 0):
                    output_dict[id_tmp][final_list[j][2][where_res[0]]] = {}
                    output_dict[id_tmp][final_list[j][2][where_res[0]]]['VIGNET'] = final_list[j][1][where_res[0]]
                    counter += 1
            if counter == 0:
                output_dict[id_tmp] = 'empty'

        return output_dict

    def process_me(self, image_dir, image_pattern, f_wcs_path, rad):
        """ Process ME

        Main function to create the stamps in the multi-epoch case.

        Parameters
        ----------
        image_dir : list
            List of directories where the image are.
            If len(image_dir) == 1 -> all images in the same directory.
            Else len(image_dir) must match len(image_pattern).
        image_pattern : list
            Common part of each kind of files.
        f_wcs_path : str
            Path to the log file containing the WCS for each CCDs.
        rad : int
            Radius of the stamp, must be odd.

        """

        self._f_wcs_file = np.load(f_wcs_path, allow_pickle=True).item()
        self._rad = rad

        for i in range(len(image_pattern)):

            if len(image_dir) != len(image_pattern):
                output_dict = self._get_stamp_me(image_dir[0],
                                                 image_pattern[i])
            else:
                output_dict = self._get_stamp_me(image_dir[i],
                                                 image_pattern[i])

            self._save_vignet_me(output_dict, image_pattern[i])

    def _save_vignet_me(self, output_dict, suffix):
        """ Save vignet ME

        Save vignets for the multi-epoch case.

        Parameters
        ----------
        output_dict : dict
            Dictionary containing object id and vignets for each epoch.
        suffix : str
            Suffix to use for the output file name.

        """
        output_name = (self._output_dir + '/' + suffix +
                       '_vignet{}'.format(self._image_num))
        np.save(output_name, output_dict)


def get_original_vignet(galcat_path):
    """Get original vignet

    Get the vignets from the SExtractor catalog

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalog

    Returns
    -------
    numpy.ndarray
        Array containing the vignets

    """

    f = io.FITSCatalog(galcat_path, SEx_catalog=True)
    f.open()

    vign = f.get_data()['VIGNET']

    f.close()

    return vign


def make_mask(galcat_path, mask_value):
    """Make mask

    Change the value of the SExtractor mask on vignet

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalog
    mask_value : float
        New value of the mask

    Returns
    -------
    numpy.array
        Array of the vignets with the new mask value

    """

    vign = get_original_vignet(galcat_path)

    vign[np.where(vign < -1e29)] = mask_value

    return vign


def save_vignet(vign, sexcat_path, output_dir, suffix, image_num):
    """Save vignet

    Save the vignet into a SExtractor format catalog

    Parameters
    ----------
    sexcat_path : str
        Path to the original SExtractor catalog
    output_dir : str
        Path to the output directory
    suffix : str
        Suffix to use for the output file name.
    image_num : str
        Image numbering.

    """

    output_name = (output_dir + '/' + suffix +
                   '_vignet{}.fits'.format(image_num))
    f = io.FITSCatalog(output_name, SEx_catalog=True,
                       open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    f.save_as_fits(vign, names=['VIGNET'], sex_cat_path=sexcat_path)


@module_runner(input_module='sextractor_runner',
               file_pattern=['galaxy_selection', 'image'],
               file_ext=['.fits', '.fits'],
               depends=['numpy', 'astropy', 'sf_tools'])
def vignetmaker_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):

    galcat_path = input_file_list[0]

    do_masking = config.getboolean("VIGNETMAKER_RUNNER", "MASKING")
    if do_masking:
        mask_value = config.getfloat("VIGNETMAKER_RUNNER", "MASK_VALUE")
        vign = make_mask(galcat_path, mask_value)
        save_vignet(vign, galcat_path, run_dirs['output'], 'cat',
                    file_number_string)

    else:
        stamp_size = config.getint("VIGNETMAKER_RUNNER", "STAMP_SIZE") - 1
        if stamp_size % 2 != 0:
            raise ValueError("The STAMP_SIZE must be odd")
        rad = int(stamp_size/2)

        pos_type = config.get("VIGNETMAKER_RUNNER", "COORD")
        pos_params = config.getlist("VIGNETMAKER_RUNNER", "POSITION_PARAMS")

        mode = config.get("VIGNETMAKER_RUNNER", "MODE")

        if mode == 'CLASSIC':
            suffix = config.getlist("VIGNETMAKER_RUNNER", "SUFFIX")
            if len(suffix) != len(input_file_list[1:]):
                raise ValueError('Number of suffixes ({}) has to be equal to '
                                 'the number of input file type ({})'
                                 ''.format(len(suffix),
                                           len(input_file_list[1:])))

            inst = vignetmaker(galcat_path, pos_type, pos_params,
                               run_dirs['output'], file_number_string)
            inst.process(input_file_list[1:], rad, suffix)
        elif mode == 'MULTI-EPOCH':
            image_dir = config.getlist("VIGNETMAKER_RUNNER", "ME_IMAGE_DIR")
            image_pattern = config.getlist("VIGNETMAKER_RUNNER",
                                           "ME_IMAGE_PATTERN")
            f_wcs_path = config.getexpanded("VIGNETMAKER_RUNNER", "ME_LOG_WCS")

            inst = vignetmaker(galcat_path, pos_type, pos_params,
                               run_dirs['output'], file_number_string)
            inst.process_me(image_dir, image_pattern, f_wcs_path, rad)
        else:
            raise ValueError('Invalid MODE=\'{}\''.format(mode))

    return None, None
