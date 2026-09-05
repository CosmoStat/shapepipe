"""VIGNET MAKER RUNNER.

Module runner for ``vignetmaker``.

:Authors: Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.vignetmaker_package import vignetmaker as vm

from shapepipe.pipeline.exp_utils import get_exp_output_dirs


@module_runner(
    version="1.1",
    input_module="sextractor_runner",
    file_pattern=["galaxy_selection", "image"],
    file_ext=[".fits", ".fits"],
    depends=["numpy", "astropy", "sf_tools", "sqlitedict"],
)
def vignetmaker_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Vingent Maker Runner."""
    # Get path to galaxy catalogue
    galcat_path = input_file_list[0]

    # Check whether mask is to be created
    if config.getboolean(module_config_sec, "MASKING"):

        # Get the mask value
        mask_value = config.getfloat(module_config_sec, "MASK_VALUE")

        # Create a mask postage stamp
        vignet = vm.make_mask(galcat_path=galcat_path, mask_value=mask_value)

        # Save the mask postage stamp vignet
        vm.save_vignet(
            vignet=vignet,
            sexcat_path=galcat_path,
            output_dir=run_dirs["output"],
            prefix="cat",
            image_num=file_number_string,
        )

    else:

        # Get stamp size
        stamp_size = config.getint(module_config_sec, "STAMP_SIZE") - 1

        # Make sure stamp size is odd number
        if stamp_size % 2 != 0:
            raise ValueError(
                f"Found even value {stamp_size} for postage stame size (STAMP"
                + f"_SIZE entry in config file), must be odd"
            )

        # Set radius
        radius = int(stamp_size / 2)

        # Get position type; allowed are "PIX" and "SPHE"
        pos_type = config.get(module_config_sec, "COORD")

        # Get position column names for x and y / ra and dec
        pos_params = config.getlist(module_config_sec, "POSITION_PARAMS")

        # Get vignet run mode; allowed are "CLASSIC", "MULTI-EPOCH"
        mode = config.get(module_config_sec, "MODE")

        # Create instance of VignetMaker
        vm_inst = vm.VignetMaker(
            galcat_path=galcat_path,
            pos_type=pos_type,
            pos_params=pos_params,
            output_dir=run_dirs["output"],
            image_num=file_number_string,
            w_log=w_log,
        )

        # Process module according to mode
        if mode == "CLASSIC":
            # Classsic mode = Single exposures

            # Get file name prefix
            prefix = config.getlist(module_config_sec, "PREFIX")

            # Make sure number of prefixes match input files
            if len(prefix) != len(input_file_list[1:]):
                raise ValueError(
                    f"The number of prefixes ({len(prefix)}) has to be "
                    + "equal to the number of input file types "
                    + f"({len(input_file_list[1:])})."
                )

            # Process inputs
            vm_inst.process(input_file_list[1:], radius, prefix)

        elif mode == "MULTI-EPOCH":
            # Multi-epoch exposures

            # Locate the runner output dirs through the per-exposure work tree.
            exp_base_dir = config.getexpanded(
                module_config_sec, "ME_IMAGE_EXP_DIR"
            )
            if len(input_file_list) < 3:
                raise ValueError(
                    "ME_IMAGE_EXP_DIR requires the exposure-numbers"
                    + " file as a third input; add 'exp_numbers' to"
                    + " FILE_PATTERN and FILE_EXT in the"
                    + f" [{module_config_sec}] config section."
                )
            exp_numbers_file = input_file_list[2]
            exp_runner_names = config.getlist(
                module_config_sec, "ME_IMAGE_EXP_RUNNERS"
            )
            image_dirs = []
            for runner_name in exp_runner_names:
                dirs = get_exp_output_dirs(
                    exp_base_dir, exp_numbers_file, runner_name, w_log
                )
                image_dirs.append(dirs)

            image_pattern = config.getlist(
                module_config_sec,
                "ME_IMAGE_PATTERN",
            )

            # Get WCS log file path. The WCS log is a positional input (via
            # FILE_PATTERN) in MULTI-EPOCH mode, not a decorator default.
            if len(input_file_list) < 2:
                raise ValueError(
                    "MULTI-EPOCH mode requires the WCS log file as a second"
                    + " input; add 'log_exp_headers' to FILE_PATTERN and"
                    + f" FILE_EXT in the [{module_config_sec}] config"
                    + " section."
                )
            f_wcs_path = input_file_list[1]

            # Process inputs
            vm_inst.process_me(image_dirs, image_pattern, f_wcs_path, radius)

        else:

            raise ValueError(f"Invalid MODE='{mode}'")

    return None, None
