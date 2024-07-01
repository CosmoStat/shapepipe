"""SEXTRACTOR SCRIPT.

This module builds the SExtractor command line.

:Author: Axel Guinot, Martin Kilbinger

"""

import re

import numpy as np
from astropy.io import fits
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io


def get_header_value(image_path, key):
    """Get Header Value.

    This function reads a value from the header image.

    Parameters
    ----------
    image_path: str
        Path to the input image
    key: str
        Key from which the value is requested (has to be a float)

    Returns
    -------
    float
        Value associated to the key provided

    """
    h = fits.getheader(image_path)

    val = h[key]

    try:
        val = float(val)
    except Exception:
        raise ValueError(
            f"The key {key} does not return a float value. Got {val}"
        )

    return val


def make_post_process(cat_path, f_wcs_path, pos_params, ccd_size):
    """Make Post Processing.

    This function will add one HDU for each epoch to the SExtractor catalogue.
    Note that this only works for tiles.

    The columns will be:

    - ``NUMBER``: same as SExtractor NUMBER
    - ``EXP_NAME``: name of the single exposure for this epoch
    - ``CCD_N``: extension where the object was detected

    Parameters
    ----------
    cat_path: str
        Path to the outputed SExtractor catalog
    f_wcs_path: str
        Path to the log file containing WCS for all single exp CCDs
    pos_params: list
        World coordinates to use to match the objects.
    ccd_size: list
        Size of a CCD ``[nx, ny]``

    Raises
    ------
    IOError
        If SQL file not found

    """
    cat = file_io.FITSCatalogue(
        cat_path,
        SEx_catalogue=True,
        open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
    )
    cat.open()

    f_wcs = SqliteDict(f_wcs_path)
    key_list = list(f_wcs.keys())
    if len(key_list) == 0:
        raise IOError(f"Could not read sql file '{f_wcs_path}'")
    n_hdu = len(f_wcs[key_list[0]])

    history = []
    for idx in cat.get_data(1)[0][0]:
        if re.split("HISTORY", idx)[0] == "":
            history.append(idx)

    exp_list = []
    pattern = r"([0-9]*)p\.(.*)"
    for hist in history:
        m = re.search(pattern, hist)
        exp_list.append(m.group(1))

    obj_id = np.copy(cat.get_data()["NUMBER"])

    n_epoch = np.zeros(len(obj_id), dtype="int32")
    for idx, exp in enumerate(exp_list):
        pos_tmp = np.ones(len(obj_id), dtype="int32") * -1
        for idx_j in range(n_hdu):
            if exp not in f_wcs:
                raise KeyError(
                    f"Exposure {exp} used in image {cat_path} but not"
                    + f" found in header file {f_wcs_path}. Make sure this"
                    + " file is complete."
                )
            w = f_wcs[exp][idx_j]["WCS"]
            pix_tmp = w.all_world2pix(
                cat.get_data()[pos_params[0]], cat.get_data()[pos_params[1]], 0
            )
            ind = (
                (pix_tmp[0] > int(ccd_size[0]))
                & (pix_tmp[0] < int(ccd_size[1]))
                & (pix_tmp[1] > int(ccd_size[2]))
                & (pix_tmp[1] < int(ccd_size[3]))
            )
            pos_tmp[ind] = idx_j
            n_epoch[ind] += 1
        exp_name = np.array([exp_list[idx] for n in range(len(obj_id))])
        a = np.array(
            [
                (obj_id[ii], exp_name[ii], pos_tmp[ii])
                for ii in range(len(exp_name))
            ],
            dtype=[
                ("NUMBER", obj_id.dtype),
                ("EXP_NAME", exp_name.dtype),
                ("CCD_N", pos_tmp.dtype),
            ],
        )
        cat.save_as_fits(data=a, ext_name=f"EPOCH_{idx}")
        cat.open()

    f_wcs.close()

    cat.add_col("N_EPOCH", n_epoch)

    cat.close()


class SExtractorCaller:
    """The SExtractor Caller.

    This class constructs the command line to call SExtractor based on the
    input files and parameters.

    Parameters
    ----------
    path_input_files: list
        List with all the paths for the input files
    path_output_dir: str
        Path for the output directory
    number_string: str
        Pipeline internal file numbering
    path_dot_sex: str
        Path to the ``.sex`` config file
    path_dot_param: str
        Path to the ``.param`` config file
    path_dot_conv: str
        Path to the ``.conv`` kernel file
    use_weight: bool
        Specify if a weight is profided for the measurement
    use_flag: bool
        Specify if a flag is provided for the measurement
    use_psf: bool
        Specify if a psf is provided for the model
    use_detection_image: bool
        Specify if a detection image is provided
    use_detection_weight: bool
        Specify if a detection weight is provided
    use_zero_point: bool
        Specify whether or not to use a zero point from the input image
    use_backgroud: bool
        Specify whether or not to use a background value form the input image
    zero_point_key: str, optional
        Header key corresponding to the zero point
    background_key: str, optional
        Header key corresponding to the background value
    check_image: str, optional
        If provided, add SExtractor check image to the output
    output_prefix: str, optional
        If provided, add a prefix to the output file names

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
        output_prefix=None,
    ):

        self.cmd_line = ""
        self._cmd_line_extra = ""

        self._meas_img_path = path_input_files[0]
        self._all_input_path = path_input_files

        self._path_output_dir = path_output_dir
        self._num_str = number_string
        self.path_output_file = self.get_output_name(output_prefix)

        self._path_dot_sex = path_dot_sex
        self._path_dot_param = path_dot_param
        self._path_dot_conv = path_dot_conv

        self.set_input_files(
            use_weight,
            use_flag,
            use_psf,
            use_detection_image,
            use_detection_weight,
        )

        # Collect optional arguments for SExtractor
        self.get_zero_point(use_zero_point, zero_point_key)
        self.get_background(use_background, background_key)
        self.get_check_image(check_image)

    def get_output_name(self, output_prefix=None):
        """Get Output Names.

        Construct the output file path.

        Parameters
        ----------
        output_prefix: str, optional
            Prefix to add to the output name

        Returns
        -------
        str
            Full path of the output file

        """
        if isinstance(output_prefix, type(None)):
            self.prefix = ""
        else:
            if (output_prefix.lower() is not None) & (output_prefix != ""):
                self.prefix = output_prefix + "_"
            else:
                self.prefix = ""

        output_file_name = self.prefix + f"sexcat{self._num_str}.fits"
        output_file_path = f"{self._path_output_dir}/{output_file_name}"

        return output_file_path

    def set_input_files(
        self,
        use_weight,
        use_flag,
        use_psf,
        use_detect_img,
        use_detect_weight,
    ):
        """Set Input Files.

        Set up all of the input image files.

        Parameters
        ----------
        use_weight: bool
            Specify if a weight is profided for the measurement
        use_flag: bool
            Specify if a flag is provided for the measurement
        use_psf: bool
            Specify if a psf is provided for the model
        use_detect_img: bool
            Specify if a detection image is provided
        use_detect_weight: bool
            Specify if a detection weight is provided

        Raises
        ------
        ValueError
            For invalid detection weight

        """
        extra = 1

        if use_weight:
            weight_image = self._all_input_path[extra]
            extra += 1

        if use_flag:
            self._cmd_line_extra += (
                " -FLAG_IMAGE " + f"{self._all_input_path[extra]}"
            )
            extra += 1

        if use_psf:
            self._cmd_line_extra += (
                " -PSF_NAME " + f"{self._all_input_path[extra]}"
            )
            extra += 1

        # Check for separate files for detection and measurement

        # First, consistency checks
        if use_detect_weight and not use_detect_img:
            raise ValueError(
                "DETECTION_WEIGHT cannot be True "
                + "if DETECTION_IMAGE is False"
            )
        if use_detect_weight and not use_weight:
            raise ValueError(
                "DETECTION_WEIGHT cannot be True " + "if WEIGHT_FILE is False"
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
                f" -WEIGHT_IMAGE {detect_weight_path}" + f",{weight_image}"
            )
        else:
            self._cmd_line_extra += " -WEIGHT_TYPE None"

        if extra != len(self._all_input_path):
            raise ValueError(
                "Incoherence between input file number and keys "
                + f"related to extra files: 1 regular + {extra-1} extra "
                + "files not compatible with total file list "
                + f"length of {len(self._all_input_path)}"
            )

    def get_zero_point(self, use_zp, zp_key=None):
        """Get Zero Point.

        Use a zero point from an input image header.

        Parameters
        ----------
        use_zp: bool
            If ``True``, add the zero point to the command line
        zp_key: str
            Header key corresponding to the zero point

        """
        if use_zp and not isinstance(zp_key, type(None)):
            zp_value = get_header_value(self._meas_img_path, zp_key)
            self._cmd_line_extra += f" -MAG_ZEROPOINT {zp_value}"

    def get_background(self, use_bkg, bkg_key=None):
        """Get Background.

        Use a background value from an input image header.

        Parameters
        ----------
        use_bkg: bool
            If ``True``, add the background value to the command line
        bkg_key: str
            Header key corresponding to the background value

        """
        if use_bkg and not isinstance(bkg_key, type(None)):
            bkg_value = get_header_value(self._meas_img_path, bkg_key)
            self._cmd_line_extra += (
                f" -BACK_TYPE MANUAL -BACK_VALUE {bkg_value}"
            )

    def get_check_image(self, check_image):
        """Get Check Image.

        Handle the check images if any are requested.

        Parameters
        ----------
        check_image: list
            List of SExtractor keys corresponding to check images

        """
        if (len(check_image) == 1) & (check_image[0] == ""):
            check_type = ["NONE"]
            check_name = ["none"]
        else:
            check_type = []
            check_name = []
            for key in check_image:
                check_type.append(key.upper())
                check_name.append(
                    self._path_output_dir
                    + "/"
                    + self.prefix
                    + key.lower()
                    + self._num_str
                    + ".fits"
                )

        self._cmd_line_extra += (
            f' -CHECKIMAGE_TYPE {",".join(check_type)} '
            + f'-CHECKIMAGE_NAME {",".join(check_name)}'
        )

    def make_command_line(self, exec_path):
        """Make Command Line.

        This method constructs the command line to run SExtractor.

        Parameters
        ----------
        exec_path: str
            Path to SExtractor executable

        Returns
        -------
        str
            Full command line to call SExtractor

        """
        # Base arguments for SExtractor
        command_line_base = (
            f"{exec_path} {self._detect_img_path},{self._meas_img_path} "
            + f"-c {self._path_dot_sex} "
            + f"-PARAMETERS_NAME {self._path_dot_param} "
            + f"-FILTER_NAME {self._path_dot_conv} "
            + f"-CATALOG_NAME {self.path_output_file}"
        )

        command_line = f"{command_line_base} {self._cmd_line_extra}"

        return command_line

    @staticmethod
    def parse_errors(stderr, stdout):
        """Parse Errors.

        This methoid moves errors from the standard output of SExtractor to
        the standard error.

        Parameters
        ----------
        stderr: str
            String containing the standard error contents
        stdout: str
            String containing the standard output content

        Returns
        -------
        tuple
            Updated standard output and error

        """
        check_error = re.findall("error", stdout.lower())
        check_error2 = re.findall("all done", stdout.lower())

        if check_error == []:
            stderr2 = ""
        else:
            stderr2 = stdout

        if check_error2 == []:
            stderr2 = stdout

        return stdout, stderr2
