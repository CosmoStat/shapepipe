"""READ EXTERNAL SEXCAT RUNNER.

Module runner for ``read_ext_sexcat``.

:Author: Martin Kilbinger

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.read_ext_sexcat_package import read_ext_sexcat as rs
from shapepipe.modules.sextractor_package import sextractor_script as ss


@module_runner(
    version="1.0",
    input_module=["get_images_runner"],
    file_pattern=["CFIS_cat", "CFIS_image"],
    file_ext=[".cat", ".fits"],
    depends=["numpy", "astropy"],
)
def read_ext_sexcat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define the Read External SExtractor Catalogue Runner.

    Reads an external ASCII catalogue (SExtractor format), converts it to
    a FITS-LDAC catalogue compatible with downstream ShapePipe modules.
    The inputs are the catalogue and the tile image, then the catalogue's
    segmentation map if SEGMENTATION = True, then the WCS log if
    MAKE_POST_PROCESS = True. With the segmentation map, neighbours' pixels
    in VIGNET are set to -1e30 and the map, relabelled to the catalogue's
    NUMBER, is written as ``seg<number>.fits``. MAKE_POST_PROCESS runs the
    multi-epoch post-processing that adds per-exposure HDUs.
    """
    cat_path, image_path, *extra_inputs = input_file_list
    use_seg = config.has_option(
        module_config_sec, "SEGMENTATION"
    ) and config.getboolean(module_config_sec, "SEGMENTATION")
    seg_path = extra_inputs.pop(0) if use_seg else None

    if config.has_option(module_config_sec, "SUFFIX"):
        suffix = config.get(module_config_sec, "SUFFIX")
    else:
        suffix = "sexcat"

    output_path = f"{run_dirs['output']}/{suffix}{file_number_string}.fits"

    if config.has_option(module_config_sec, "VIGNET_SIZE"):
        stamp_size = config.getint(module_config_sec, "VIGNET_SIZE")
    else:
        stamp_size = 51

    w_log.info(f"Reading external catalogue: {cat_path}")
    w_log.info(f"Reading image header from: {image_path}")

    rs.make_ldac_from_ascii(
        cat_path,
        image_path,
        output_path,
        stamp_size=stamp_size,
        seg_path=seg_path,
        w_log=w_log,
    )

    if config.getboolean(module_config_sec, "MAKE_POST_PROCESS"):
        # The WCS log is supplied as a positional input (via FILE_PATTERN)
        # when post-processing is enabled, not by the decorator default.
        if not extra_inputs:
            raise ValueError(
                "MAKE_POST_PROCESS requires the WCS log file as the last"
                + " input; add 'log_exp_headers' to FILE_PATTERN and"
                + f" FILE_EXT in the [{module_config_sec}] config section."
            )
        f_wcs_path = extra_inputs[0]
        pos_params = config.getlist(module_config_sec, "WORLD_POSITION")
        ccd_size = config.getlist(module_config_sec, "CCD_SIZE")
        w_log.info("Running post-processing")
        ss.make_post_process(output_path, f_wcs_path, pos_params, ccd_size)

    return None, None
