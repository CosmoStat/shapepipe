"""SEXTRACTOR RUNNER.

Module runner for ``sextractor``.

:Author: Axel Guinot

"""

import contextlib
import os
import re
import shutil
import tempfile

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.sextractor_package import match_catalogue as mc
from shapepipe.modules.sextractor_package import sextractor_script as ss
from shapepipe.pipeline.execute import execute


# The trailing log_exp_headers input (merged WCS headers from
# merge_headers_runner) is only consumed when MAKE_POST_PROCESS is True;
# configs without post-processing override FILE_PATTERN/FILE_EXT with the
# first three entries only.
@module_runner(
    version="1.0.1",
    input_module=["split_exp_runner", "merge_headers_runner"],
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

    # SEG_VIGNET (optional, environment-expanded boolean): add the
    # SEGMENTATION check image's stamps, on each VIGNET's grid, as the
    # SEG_VIGNET column ngmix's UberSeg blend handling reads. SExtractor then
    # also writes the double-precision positions VIGNET is centred on, which
    # add_seg_vignet reads and drops.
    seg_vignet = config.has_option(
        module_config_sec, "SEG_VIGNET"
    ) and config.getexpandedboolean(module_config_sec, "SEG_VIGNET")
    if seg_vignet:
        if "SEGMENTATION" not in [key.upper() for key in check_image]:
            raise ValueError(
                "SEG_VIGNET needs the SEGMENTATION check image in CHECKIMAGE."
            )
        dot_param = ss.seg_vignet_param_file(
            dot_param,
            f"{run_dirs['tmp']}/seg_vignet{file_number_string}.param",
        )

    # WORK_DIR (optional, environment-expanded): SExtractor writes its
    # catalogue and check images in a fresh directory there, the SEG_VIGNET
    # cut, the join and the post-processing each rewrite the catalogue there,
    # and the files are moved to the run's output directory once complete.
    # On node-local disk this keeps the whole-catalogue rewrites off NFS,
    # where each costs minutes under load. The directory is removed on the way
    # out, also on failure. Without WORK_DIR everything is written in the
    # output directory.
    if config.has_option(module_config_sec, "WORK_DIR"):
        work = tempfile.TemporaryDirectory(
            prefix=f"sp-detect{file_number_string}.",
            dir=config.getexpanded(module_config_sec, "WORK_DIR"),
        )
    else:
        work = contextlib.nullcontext(run_dirs["output"])

    with work as work_dir:
        # Create sextractor caller class instance
        ss_inst = ss.SExtractorCaller(
            input_file_list,
            work_dir,
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

        # SEG_VIGNET is cut before the join, which relabels its stamps with the
        # rows' new NUMBERs.
        if seg_vignet:
            ss.add_seg_vignet(
                ss_inst.path_output_file,
                ss_inst.check_paths["SEGMENTATION"],
                w_log=w_log,
            )

        # MATCH_CATALOGUE (optional, environment-expanded path; empty for none):
        # take membership and NUMBER from that external catalogue of the same
        # image, before the post-processing keys the epoch HDUs on NUMBER.
        match_path = (
            config.getexpanded(module_config_sec, "MATCH_CATALOGUE")
            if config.has_option(module_config_sec, "MATCH_CATALOGUE")
            else ""
        )
        if match_path:
            mc.match_catalogue(
                ss_inst.path_output_file,
                match_path,
                radius=config.getfloat(module_config_sec, "MATCH_RADIUS"),
                min_fraction=config.getfloat(
                    module_config_sec, "MATCH_MIN_FRACTION"
                ),
                tolerated_unpaired=config.getint(
                    module_config_sec, "MATCH_TOLERATED_UNPAIRED"
                ),
                w_log=w_log,
            )

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

        if work_dir != run_dirs["output"]:
            for path in [*ss_inst.check_paths.values(),
                         ss_inst.path_output_file]:
                shutil.move(
                    path,
                    os.path.join(run_dirs["output"], os.path.basename(path)),
                )

    # Return stdout and stderr
    return stdout, stderr
