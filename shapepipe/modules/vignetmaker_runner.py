"""VIGNET MAKER RUNNER.

Module runner for ``vignetmaker``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.vignetmaker_package import vignetmaker as vm

from shapepipe.pipeline.run_log import get_last_dir


@module_runner(
    version='1.1',
    input_module='sextractor_runner',
    file_pattern=['galaxy_selection', 'image'],
    file_ext=['.fits', '.fits'],
    depends=['numpy', 'astropy', 'sf_tools', 'sqlitedict'],
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

    # Check if masking should be performed
    # With masking
    if config.getboolean(module_config_sec, 'MASKING'):
        # Fetch the mask value
        mask_value = config.getfloat(module_config_sec, 'MASK_VALUE')
        # Make a mask
        vignet = vm.make_mask(galcat_path=galcat_path, mask_value=mask_value)
        # Save the vignet
        vm.save_vignet(
            vignet=vignet,
            sexcat_path=galcat_path,
            output_dir=run_dirs['output'],
            suffix='cat',
            image_num=file_number_string,
        )

    # Without masking
    else:
        # Fetch stamp size
        stamp_size = config.getint(module_config_sec, 'STAMP_SIZE') - 1
        # Check stamp size
        if stamp_size % 2 != 0:
            raise ValueError('The STAMP_SIZE must be odd')
        # Set radius
        radius = int(stamp_size / 2)

        # Fetch position type and values
        pos_type = config.get(module_config_sec, 'COORD')
        pos_params = config.getlist(module_config_sec, 'POSITION_PARAMS')
        # Fetch vignet run mode
        mode = config.get(module_config_sec, 'MODE')

        # Create instance of VignetMaker
        vm_inst = vm.VignetMaker(
            galcat_path=galcat_path,
            pos_type=pos_type,
            pos_params=pos_params,
            output_dir=run_dirs['output'],
            image_num=file_number_string,
        )

        # Run in CLASSIC mode
        if mode == 'CLASSIC':
            # Fetch suffix
            suffix = config.getlist(module_config_sec, 'SUFFIX')
            # Check suffix
            if len(suffix) != len(input_file_list[1:]):
                raise ValueError(
                    f'The number of suffixes ({len(suffix)}) has to be '
                    + 'equal to the number of input file types '
                    + f'({len(input_file_list[1:])}).'
                )

            # Process inputs
            vm_inst.process(input_file_list[1:], radius, suffix)

        # Run in MULTI-EPOCH mode
        elif mode == 'MULTI-EPOCH':
            # Fetch image directory and patterns
            modules = config.getlist(module_config_sec, 'ME_IMAGE_DIR')
            image_dir = []
            for module in modules:
                last_dir = get_last_dir(run_dirs['run_log'], module)
                image_dir.append(last_dir)
            image_pattern = config.getlist(
                module_config_sec,
                'ME_IMAGE_PATTERN',
            )
            # Fetch WCS log path
            f_wcs_path = config.getexpanded(module_config_sec, 'ME_LOG_WCS')

            # Process inputs
            vm_inst.process_me(image_dir, image_pattern, f_wcs_path, radius)

        # Invalid mode
        else:
            # Raise error for invalid run mode
            raise ValueError(f'Invalid MODE=\'{mode}\'')

    # No return objects
    return None, None
