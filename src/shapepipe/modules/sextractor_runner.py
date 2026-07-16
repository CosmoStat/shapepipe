"""SEXTRACTOR RUNNER.

Module runner for ``sextractor``.

:Author: Axel Guinot

"""

import re

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.sextractor_package import sextractor_script as ss
from shapepipe.pipeline.execute import execute


# The trailing log_exp_headers input (merged WCS headers from
# merge_headers_runner) is only consumed when MAKE_POST_PROCESS is True;
# configs without post-processing override FILE_PATTERN/FILE_EXT with the
# first three entries only.
@module_runner(
    version="1.0.1",
    input_module=["mask_runner", "merge_headers_runner"],
    file_pattern=["image", "weight", "flag", "log_exp_headers"],
    file_ext=[".fits", ".fits", ".fits", ".sqlite"],
    executes=["source-extractor"],
    depends=["numpy"],
)
def sextractor_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The SExtractor Runner."""
    # Set the SExtractor executable name
    if config.has_option(module_config_sec, "EXEC_PATH"):
        exec_path = config.getexpanded(module_config_sec, "EXEC_PATH")
    else:
        exec_path = "sex"

    # Get SExtractor config options
    dot_sex = config.getexpanded(module_config_sec, "DOT_SEX_FILE")
    dot_param = config.getexpanded(module_config_sec, "DOT_PARAM_FILE")
    dot_conv = config.getexpanded(module_config_sec, "DOT_CONV_FILE")
    weight_file = config.getboolean(module_config_sec, "WEIGHT_IMAGE")
    flag_file = config.getboolean(module_config_sec, "FLAG_IMAGE")
    psf_file = config.getboolean(module_config_sec, "PSF_FILE")
    detection_image = config.getboolean(module_config_sec, "DETECTION_IMAGE")
    detection_weight = config.getboolean(module_config_sec, "DETECTION_WEIGHT")

    zp_from_header = config.getboolean(module_config_sec, "ZP_FROM_HEADER")
    if zp_from_header:
        zp_key = config.get(module_config_sec, "ZP_KEY")
    else:
        zp_key = None

    bkg_from_header = config.getboolean(module_config_sec, "BKG_FROM_HEADER")
    if bkg_from_header:
        bkg_key = config.get(module_config_sec, "BKG_KEY")
    else:
        bkg_key = None

    if config.has_option(module_config_sec, "CHECKIMAGE"):
        check_image = config.getlist(module_config_sec, "CHECKIMAGE")
    else:
        check_image = [""]

    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(module_config_sec, "PREFIX")
    else:
        prefix = None

    # When post-processing is enabled the sqlite WCS file is the last input;
    # remove it before passing the image files to SExtractorCaller.
    if config.getboolean(module_config_sec, "MAKE_POST_PROCESS"):
        f_wcs_path = input_file_list[-1]
        input_file_list = list(input_file_list[:-1])

    # Create sextractor caller class instance
    ss_inst = ss.SExtractorCaller(
        input_file_list,
        run_dirs["output"],
        file_number_string,
        dot_sex,
        dot_param,
        dot_conv,
        weight_file,
        flag_file,
        psf_file,
        detection_image,
        detection_weight,
        zp_from_header,
        bkg_from_header,
        zero_point_key=zp_key,
        background_key=bkg_key,
        check_image=check_image,
        output_prefix=prefix,
    )

    # Generate sextractor command line
    command_line = ss_inst.make_command_line(exec_path)
    w_log.info(f"Calling command: {command_line}")

    # Execute command line
    stderr, stdout = execute(command_line)

    # Parse SExtractor errors
    stdout, stderr = ss_inst.parse_errors(stderr, stdout)

    # Run sextractor post processing
    if config.getboolean(module_config_sec, "MAKE_POST_PROCESS"):
        pos_params = config.getlist(module_config_sec, "WORLD_POSITION")
        ccd_size = config.getlist(module_config_sec, "CCD_SIZE")
        ss.make_post_process(
            ss_inst.path_output_file,
            f_wcs_path,
            pos_params,
            ccd_size,
            w_log=w_log,
        )

    # Return stdout and stderr
    return stdout, stderr
