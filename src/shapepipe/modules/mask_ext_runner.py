"""MASK EXT RUNNER.

Module runner for ``mask_ext``.

:Author: Cail Daley

"""

from shapepipe.modules.mask_ext_package.mask_ext import MaskExt
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.0",
    file_pattern=["image", "flag"],
    file_ext=[".fits", ".fits"],
    depends=["numpy", "astropy", "healsparse"],
    numbering_scheme="_0",
)
def mask_ext_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Mask Ext Runner.

    Rasterize an external healsparse mask onto the input image pixel grid and
    write the ShapePipe flag file.

    Notes
    -----
    Only the image is strictly required: it supplies the pixel grid and WCS
    onto which the healsparse mask is rasterized. Unlike the ``mask`` module,
    no weight file is used — the healsparse mask already encodes the footprint,
    and missing-data handling is not this module's concern. A second input, an
    external instrument flag file, is consumed only when ``USE_EXT_FLAG`` is
    ``True`` (needed for exposures, where saturation / bleeding / bad columns
    arrive with the data).

    """
    n_inputs = len(input_file_list)
    use_ext_flag = config.getboolean(module_config_sec, "USE_EXT_FLAG") if (
        config.has_option(module_config_sec, "USE_EXT_FLAG")
    ) else False

    if use_ext_flag:
        if n_inputs != 2:
            raise ValueError(
                f"Found {n_inputs} inputs but USE_EXT_FLAG is True, which "
                + 'expects "image" and "flag" in the MASK_EXT_RUNNER section '
                + "of the config file."
            )
        image_path = input_file_list[0]
        ext_flag_name = input_file_list[1]
    else:
        if n_inputs != 1:
            raise ValueError(
                f"Found {n_inputs} inputs but USE_EXT_FLAG is False, which "
                + 'expects only "image" in the MASK_EXT_RUNNER section of the '
                + "config file."
            )
        image_path = input_file_list[0]
        ext_flag_name = None

    # Path to the healsparse mask file
    mask_path = config.getexpanded(module_config_sec, "MASK_PATH")

    # Bit -> flag mapping
    bit_flag_map = MaskExt.parse_bit_flag_map(
        config.get(module_config_sec, "BIT_FLAG_MAP")
    )

    # Flag value for pixels outside the healsparse footprint
    if config.has_option(module_config_sec, "OFF_MAP_FLAG"):
        off_map_flag = config.getint(module_config_sec, "OFF_MAP_FLAG")
    else:
        off_map_flag = 0

    # HDU of the external instrument flag file
    if config.has_option(module_config_sec, "HDU"):
        hdu = config.getint(module_config_sec, "HDU")
    else:
        hdu = 0

    # Output file name prefix
    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(module_config_sec, "PREFIX")
    else:
        prefix = ""

    # Rows per chunk (optional; default derived from image size)
    if config.has_option(module_config_sec, "CHUNK_SIZE"):
        chunk_size = config.getint(module_config_sec, "CHUNK_SIZE")
    else:
        chunk_size = None

    mask_inst = MaskExt(
        image_path,
        mask_path,
        bit_flag_map,
        file_number_string,
        run_dirs["output"],
        w_log,
        off_map_flag=off_map_flag,
        path_external_flag=ext_flag_name,
        image_prefix=prefix.replace(" ", ""),
        outname_base="flag",
        chunk_size=chunk_size,
        hdu=hdu,
    )

    mask_inst.make_mask()

    return None, None
