# -*- coding: utf-8 -*-

"""SEXTRACTOR SCRIPT

This module builds the SExtractor command line.

:Author: Axel Guinot & Martin Kilbinger

"""

import re
from shapepipe.pipeline import file_io as io

import numpy as np
from sqlitedict import SqliteDict
from astropy.io import fits


def get_header_value(image_path, key):
    """Get header value

    This function reads a value from the header image.

    Parameters
    ----------
    image_path: str
        Path to the input image
    key: str
        Key from which the value is requested (has to be float)

    Returns
    -------
    val: float
        Value associated to the key provided

    """

    h = fits.getheader(image_path)

    val = h[key]

    try:
        val = float(val)
    except Exception:
        raise ValueError(
            f'The key {key} does not return a float value. Got {val}'
        )

    return val


def make_post_process(cat_path, f_wcs_path, pos_params, ccd_size):
    """Make post process

    This function will add one hdu by epoch to the SExtractor catalog.
    Only works for tiles.
    The columns will be: NUMBER same as SExtractor NUMBER
                         EXP_NAME name of the single exposure for this epoch
                         CCD_N extension where the object is

    Parameters
    ----------
    cat_path: str
        Path to the outputed SExtractor catalog
    f_wcs_path: str
        Path to the log file containing wcs for all single exp CCDs
    pos_params: list
        World coordinates to use to match the objects.
    ccd_size: list
        Size of a ccd [nx, ny]

    Raises
    ------
    IOError

    """

    cat = io.FITSCatalog(
        cat_path, SEx_catalog=True,
        open_mode=io.BaseCatalog.OpenMode.ReadWrite,
    )
    cat.open()

    f_wcs = SqliteDict(f_wcs_path)
    key_list = list(f_wcs.keys())
    if len(key_list) == 0:
        raise IOError(f'Could not read sql file \'{f_wcs_path}\'')
    n_hdu = len(f_wcs[key_list[0]])

    history = []
    for i in cat.get_data(1)[0][0]:
        if re.split('HISTORY', i)[0] == '':
            history.append(i)

    exp_list = []
    pattern = r'([0-9]*)p\.(.*)'
    for hist in history:
        m = re.search(pattern, hist)
        exp_list.append(m.group(1))

    obj_id = np.copy(cat.get_data()['NUMBER'])

    n_epoch = np.zeros(len(obj_id), dtype='int32')
    for idx, exp in enumerate(exp_list):
        pos_tmp = np.ones(len(obj_id), dtype='int32') * -1
        for j in range(n_hdu):
            w = f_wcs[exp][j]['WCS']
            pix_tmp = w.all_world2pix(cat.get_data()[pos_params[0]],
                                      cat.get_data()[pos_params[1]], 0)
            ind = ((pix_tmp[0] > int(ccd_size[0])) &
                   (pix_tmp[0] < int(ccd_size[1])) &
                   (pix_tmp[1] > int(ccd_size[2])) &
                   (pix_tmp[1] < int(ccd_size[3])))
            pos_tmp[ind] = j
            n_epoch[ind] += 1
        exp_name = np.array([exp_list[idx] for n in range(len(obj_id))])
        a = np.array([(obj_id[ii], exp_name[ii], pos_tmp[ii])
                      for ii in range(len(exp_name))],
                     dtype=[('NUMBER', obj_id.dtype),
                            ('EXP_NAME', exp_name.dtype),
                            ('CCD_N', pos_tmp.dtype)])
        cat.save_as_fits(data=a, ext_name='EPOCH_{}'.format(idx))
        cat.open()

    f_wcs.close()

    cat.add_col('N_EPOCH', n_epoch)

    cat.close()


class sextractor_caller():
    """SExtractor Caller

    This class constructs the command line to call SExtractor based on the
    input files and parameters.

    Parameters
    ----------
    path_input_files: list
        List with all the path for the input files
    path_output_dir: str
        Path for the output directory
    number_string: str
        Pipeline intern numerotation
    path_dot_sex: str
        Path to the ".sex" config file
    path_dot_param: str
        Path to the ".param" config file
    path_dot_conv: str
        Path to the ".conv" kernel file
    use_weight: bool
        Weither a weight is profided for the measurement
    use_flag: bool
        Weither a flag is provided for the measurement
    use_psf: bool
        Weither a psf is provided for the model
    use_detection_image: bool
        Weither a detection image is provided
    use_detection_weight: bool
        Weither a detection weight is provided
    use_zero_point: bool
        Weither to use a zero point from the input image
    zero_point_key: str
        Header key corresponding to the zero point
    use_backgroup: bool
        Weither to use a background value form the input image
    background_key: str
        Header key corresponding to the background value
    check_image: str
        If provided, add SExtractor check image to the output
    output_suffix: str
        If provided, add a suffix to the output files

    """
    def __init__(
            self,
            path_input_files,
            path_output_dir,
            number_string,
            path_dot_sex,
            path_dot_param,
            path_dot_conv,
            use_weight,
            use_flag,
            use_psf,
            use_detection_image,
            use_detection_weight,
            use_zero_point,
            use_background,
            zero_point_key=None,
            background_key=None,
            check_image=None,
            output_suffix=None
    ):

        self.cmd_line = ''
        self._cmd_line_extra = ''

        self._meas_img_path = path_input_files[0]
        self._all_input_path = path_input_files

        self._path_output_dir = path_output_dir
        self._num_str = number_string
        self.path_output_file = self.get_output_name(output_suffix)

        self._path_dot_sex = path_dot_sex
        self._path_dot_param = path_dot_param
        self._path_dot_conv = path_dot_conv

        self.set_input_files(use_weight, use_flag, use_psf,
                             use_detection_image, use_detection_weight)

        # Collect optional arguments for SExtractor
        self.get_zero_point(use_zero_point, zero_point_key)
        self.get_background(use_background, background_key)
        self.get_check_image(check_image)

    def get_output_name(self, output_suffix=None):
        """Get output names

        Construct the output file path.

        Parameters
        ----------
        output_suffix: str
            Suffix to add to the output name, can be None

        Returns
        -------
        output_file_path: str
            Full path of the output file

        """
        if isinstance(output_suffix, type(None)):
            self.suffix = ''
        else:
            if (output_suffix.lower() is not None) & (output_suffix != ''):
                self.suffix = output_suffix + '_'
            else:
                self.suffix = ''

        output_file_name = self.suffix + f'sexcat{self._num_str}.fits'
        output_file_path = f'{self._path_output_dir}/{output_file_name}'

        return output_file_path

    def set_input_files(
            self,
            use_weight,
            use_flag,
            use_psf,
            use_detect_img,
            use_detect_weight
    ):
        """Set input files

        Setup all the input image files.

        Parameters
        ----------
        use_weight: bool
            Weither a weight is profided for the measurement
        use_flag: bool
            Weither a flag is provided for the measurement
        use_psf: bool
            Weither a psf is provided for the model
        use_detect_img: bool
            Weither a detection image is provided
        use_detect_weight: bool
            Weither a detection weight is provided

        Raise
        -----
        ValueError

        """
        extra = 1

        if use_weight:
            weight_image = self._all_input_path[extra]
            extra += 1

        if use_flag:
            self._cmd_line_extra += (
                ' -FLAG_IMAGE '
                + f'{self._all_input_path[extra]}'
            )
            extra += 1

        if use_psf:
            self._cmd_line_extra += (
                ' -PSF_NAME '
                + f'{self._all_input_path[extra]}'
            )
            extra += 1

        # Check for separate files for detection and measurement

        # First, consistency checks
        if use_detect_weight and not use_detect_img:
            raise ValueError(
                'DETECTION_WEIGHT cannot be True '
                + 'if DETECTION_IMAGE is False'
            )
        if use_detect_weight and not use_weight:
            raise ValueError(
                'DETECTION_WEIGHT cannot be True '
                + 'if WEIGHT_FILE is False'
            )

        # Check for separate image file for detection and measurement
        if use_detect_img:
            self._detect_img_path = self._all_input_path[extra]
            extra += 1
        else:
            self._detect_img_path = self._meas_img_path

        # Check for separate weight file corresponding to the detection image.
        # If False, use measurement weight image.
        # Note: This could be changed, and no weight image could be used, but
        # this might lead to more user errors.
        if use_weight:
            if use_detect_weight:
                detect_weight_path = self._all_input_path[extra]
                extra += 1
            else:
                detect_weight_path = weight_image
            self._cmd_line_extra += (
                f' -WEIGHT_IMAGE {detect_weight_path}'
                + f',{weight_image}'
            )
        else:
            self._cmd_line_extra += ' -WEIGHT_TYPE None'

        if extra != len(self._all_input_path):
            raise ValueError(
                'Incoherence between input file number and keys '
                + f'related to extra files: 1 regular + {extra-1} extra '
                + 'files not compatible with total file list '
                + f'length of {len(self._all_input_path)}'
            )

    def get_zero_point(self, use_zp, zp_key=None):
        """Get Zero Point

        Use a zero point from input image header.

        Parameters
        ----------
        use_zp: bool
            If True, add the zero point to the command line
        zp_key: str
            Header key corresponding to the zero point

        """
        if use_zp and isinstance(zp_key, type(None)):
            zp_value = get_header_value(self._meas_img_path, zp_key)
            self._cmd_line_extra += f' -MAG_ZEROPOINT {zp_value}'

    def get_background(self, use_bkg, bkg_key=None):
        """Get Background

        Use a background value from input image header.

        Parameters
        ----------
        use_bkg: bool
            If True, add the background value to the command line
        bkg_key: str
            Header key corresponding to the background value

        """
        if use_bkg and isinstance(bkg_key, type(None)):
            bkg_value = get_header_value(self._meas_img_path, bkg_key)
            self._cmd_line_extra += (
                f' -BACK_TYPE MANUAL -BACK_VALUE {bkg_value}'
            )

    def get_check_image(self, check_image):
        """Get check image

        Handle the check images if any are requested.

        Parameters
        ----------
        check_image: list
            List of SExtractor keys corresponding to check images

        """

        if (len(check_image) == 1) & (check_image[0] == ''):
            check_type = ['NONE']
            check_name = ['none']
        else:
            check_type = []
            check_name = []
            for key in check_image:
                check_type.append(key.upper())
                check_name.append(
                    self._path_output_dir + '/' + self.suffix
                    + key.lower()
                    + self._num_str + '.fits'
                )

        self._cmd_line_extra += (
            f" -CHECKIMAGE_TYPE {','.join(check_type)} "
            + f"-CHECKIMAGE_NAME {','.join(check_name)}"
        )

    def make_command_line(self, exec_path):
        """ Make command line

        Main that construct the command line to run SExtractor

        Parameters
        ----------
        exec_path: str
            Path to SExtractor executable

        Return
        ------
        command_line: str
            Full command line to call SExtractor

        """

        # Base arguments for SExtractor
        command_line_base = (
            f'{exec_path} {self._detect_img_path},{self._meas_img_path} '
            + f'-c {self._path_dot_sex} '
            + f'-PARAMETERS_NAME {self._path_dot_param}'
            + f' -FILTER_NAME {self._path_dot_conv} '
            + f'-CATALOG_NAME {self.path_output_file}'
        )

        command_line = f'{command_line_base} {self._cmd_line_extra}'

        return command_line
