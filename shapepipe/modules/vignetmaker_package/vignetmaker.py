"""VIGNET MAKER.

This module contains a class to create postage stamps from images.

:Author: Axel Guinot

"""

import re

import numpy as np
from astropy.wcs import WCS
from sf_tools.image.stamp import FetchStamps
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io


class VignetMaker(object):
    """Vignet Maker.

    This class handles the creation of vignets.

    Parameters
    ----------
    galcat_path : str
        Path to catalogue containing the positions.
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

    def __init__(
        self,
        galcat_path,
        pos_type,
        pos_params,
        output_dir,
        image_num,
    ):

        self._galcat_path = galcat_path
        self._output_dir = output_dir
        self._image_num = image_num
        self._pos = self.get_pos(pos_params)
        self._pos_type = pos_type

    def process(self, image_path_list, rad, suffix):
        """Process.

        Main function to create the stamps.

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
                raise ValueError(
                    'Coordinates type must be in : PIX (pixel), '
                    + 'SPHE (spherical).'
                )

            vign = self._get_stamp(image_path, pos - 1, rad)

            save_vignet(
                vign,
                self._galcat_path,
                self._output_dir,
                _suffix,
                self._image_num,
            )

    def get_pos(self, pos_params):
        """Get Positions.

        Get the positions of the given parameters from SExtractor catalogue.

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
        file = file_io.FITSCatalogue(self._galcat_path, SEx_catalogue=True)
        file.open()

        pos = np.array(
            [file.get_data()[pos_params[1]], file.get_data()[pos_params[0]]]
        ).T

        file.close()

        return pos

    def convert_pos(self, image_path):
        """Convert Position.

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
        file = file_io.FITSCatalogue(image_path)
        file.open()
        head = file.get_header(0)
        file.close()

        wcs = WCS(head)

        pos_tmp = np.copy(self._pos)
        pos_tmp[:, [0, 1]] = pos_tmp[:, [1, 0]]

        new_pos = wcs.all_world2pix(pos_tmp, 1)

        new_pos[:, [0, 1]] = new_pos[:, [1, 0]]

        return new_pos

    def _get_stamp(self, img_path, pos, rad):
        """Get Stamp.

        Extract stamps at given positions in the image.

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
        img_file = file_io.FITSCatalogue(img_path)
        img_file.open()
        img = img_file.get_data(0)
        img_file.close()

        fs = FetchStamps(img, int(rad))
        fs.get_pixels(np.round(pos).astype(int))

        vign = fs.scan()

        return vign

    def _get_stamp_me(self, image_dir, image_pattern):
        """Get Stamp Multi-Epoch.

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
        cat = file_io.FITSCatalogue(self._galcat_path, SEx_catalogue=True)
        cat.open()

        all_id = np.copy(cat.get_data()['NUMBER'])
        n_epoch = np.copy(cat.get_data()['N_EPOCH'])

        list_ext_name = cat.get_ext_name()
        hdu_ind = [
            i for i in range(len(list_ext_name))
            if 'EPOCH' in list_ext_name[i]
        ]

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

                img_path = (
                    image_dir + '/' + image_pattern + '-'
                    + exp_name + '-' + str(ccd) + '.fits'
                )
                ind_obj = np.where(cat.get_data(hdu_index)['CCD_N'] == ccd)[0]
                obj_id = all_id[ind_obj]

                wcs_file = self._f_wcs_file[exp_name][ccd]['WCS']
                pos = np.array(wcs_file.all_world2pix(
                    self._pos[:, 1][ind_obj],
                    self._pos[:, 0][ind_obj],
                    1,
                )).T
                pos[:, [0, 1]] = pos[:, [1, 0]]

                tmp_vign = self._get_stamp(img_path, pos - 1, self._rad)

                if array_vign is None:
                    array_vign = np.copy(tmp_vign)
                else:
                    array_vign = np.concatenate(
                        (array_vign, np.copy(tmp_vign))
                    )

                if array_id is None:
                    array_id = np.copy(obj_id)
                else:
                    array_id = np.concatenate((array_id, np.copy(obj_id)))

                exp_name_tmp = np.array(
                    [exp_name + '-' + str(ccd) for i in range(len(obj_id))]
                )
                if array_exp_name is None:
                    array_exp_name = exp_name_tmp
                else:
                    array_exp_name = np.concatenate(
                        (array_exp_name, exp_name_tmp)
                    )

            final_list.append([array_id, array_vign, array_exp_name])

        cat.close()

        output_dict = {}
        for id_tmp in all_id:

            output_dict[id_tmp] = {}
            counter = 0

            for j in range(len(final_list)):

                where_res = np.where(final_list[j][0] == id_tmp)[0]

                if (len(where_res) != 0):
                    index = final_list[j][2][where_res[0]]
                    output_dict[id_tmp][index] = {}
                    output_dict[id_tmp][index]['VIGNET'] = (
                        final_list[j][1][where_res[0]]
                    )
                    counter += 1

            if counter == 0:
                output_dict[id_tmp] = 'empty'

        return output_dict

    def process_me(self, image_dir, image_pattern, f_wcs_path, rad):
        """Process Multi-Epoch.

        Main function to create the stamps in the multi-epoch case.

        Parameters
        ----------
        image_dir : list
            List of directories where the image are.9
            If len(image_dir) == 1 -> all images in the same directory.
            Else len(image_dir) must match len(image_pattern).
        image_pattern : list
            Common part of each kind of files.
        f_wcs_path : str
            Path to the log file containing the WCS for each CCDs.
        rad : int
            Radius of the stamp, must be odd.

        """
        self._f_wcs_file = SqliteDict(f_wcs_path)
        self._rad = rad

        for i in range(len(image_pattern)):

            if len(image_dir) != len(image_pattern):
                output_dict = self._get_stamp_me(
                    image_dir[0],
                    image_pattern[i],
                )

            else:
                output_dict = self._get_stamp_me(
                    image_dir[i],
                    image_pattern[i],
                )

            self._save_vignet_me(output_dict, image_pattern[i])

        self._f_wcs_file.close()

    def _save_vignet_me(self, output_dict, suffix):
        """Save vignet Multi-Epoch.

        Save vignets for the multi-epoch case.

        Parameters
        ----------
        output_dict : dict
            Dictionary containing object id and vignets for each epoch.
        suffix : str
            Suffix to use for the output file name.

        """
        output_name = f'{self._output_dir}/{suffix}_vignet{self._image_num}'
        output_file = SqliteDict(output_name + '.sqlite')

        for _index in output_dict.keys():
            output_file[str(_index)] = output_dict[_index]

        output_file.commit()
        output_file.close()


def get_original_vignet(galcat_path):
    """Get Original Vignet.

    Get the vignets from the SExtractor catalogue.

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalogue

    Returns
    -------
    numpy.ndarray
        Array containing the vignets

    """
    file = file_io.FITSCatalogue(galcat_path, SEx_catalogue=True)
    file.open()

    vignet = file.get_data()['VIGNET']

    file.close()

    return vignet


def make_mask(galcat_path, mask_value):
    """Make Mask.

    Change the value of the SExtractor mask in the vignet.

    Parameters
    ----------
    galcat_path : str
        Path to the SExtractor catalogue
    mask_value : float
        New value of the mask

    Returns
    -------
    numpy.array
        Array of the vignets with the new mask value

    """
    vignet = get_original_vignet(galcat_path)

    vignet[np.where(vignet < -1e29)] = mask_value

    return vignet


def save_vignet(vignet, sexcat_path, output_dir, suffix, image_num):
    """Save Vignet.

    Save the vignet to a SExtractor format catalogue.

    Parameters
    ----------
    vignet : numpy.ndarray
        Array containing the vignets to be saved
    sexcat_path : str
        Path to the original SExtractor catalog
    output_dir : str
        Path to the output directory
    suffix : str
        Suffix to use for the output file name.
    image_num : str
        Image numbering.

    """
    output_name = f'{output_dir}/{suffix}_vignet{image_num}.fits'

    file = file_io.FITSCatalogue(
        output_name,
        SEx_catalogue=True,
        open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
    )

    file.save_as_fits(vignet, names=['VIGNET'], sex_cat_path=sexcat_path)
