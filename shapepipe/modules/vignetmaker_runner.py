# -*- coding: utf-8 -*-

"""VIGNET MAKER RUNNER

Module runner for ``vignetmaker``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.vignetmaker_package import vignetmaker as vm


@module_runner(
    input_module='sextractor_runner',
    file_pattern=['galaxy_selection', 'image'],
    file_ext=['.fits', '.fits'],
    depends=['numpy', 'astropy', 'sf_tools', 'sqlitedict']
)
def vignetmaker_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    w_log
):

    # Get path to galaxy catalogue
    galcat_path = input_file_list[0]

    # Check if masking should be performed
    # With masking
    if config.getboolean('VIGNETMAKER_RUNNER', 'MASKING'):
        # Fetch the mask value
        mask_value = config.getfloat('VIGNETMAKER_RUNNER', 'MASK_VALUE')
        # Make a mask
        vignet = vm.make_mask(galcat_path, mask_value)
        # Save the vignet
        vm.save_vignet(
            vignet,
            galcat_path,
            run_dirs['output'],
            'cat',
            file_number_string,
        )

    # Without masking
    else:

        # Fetch stamp size
        stamp_size = config.getint('VIGNETMAKER_RUNNER', 'STAMP_SIZE') - 1
        # Check stamp size
        if stamp_size % 2 != 0:
            raise ValueError('The STAMP_SIZE must be odd')
        # Set radius
        rad = int(stamp_size / 2)

        # Fetch position type and values
        pos_type = config.get('VIGNETMAKER_RUNNER', 'COORD')
        pos_params = config.getlist('VIGNETMAKER_RUNNER', 'POSITION_PARAMS')
        # Fetch vignet run mode
        mode = config.get('VIGNETMAKER_RUNNER', 'MODE')

        # Create instance of VignetMaker
        vm_inst = vm.VignetMaker(
            galcat_path,
            pos_type,
            pos_params,
            run_dirs['output'],
            file_number_string,
        )

        # Run in CLASSIC mode
        if mode == 'CLASSIC':
            # Fetch suffix
            suffix = config.getlist('VIGNETMAKER_RUNNER', 'SUFFIX')
            # Check suffix
            if len(suffix) != len(input_file_list[1:]):
                raise ValueError(
                    'Number of suffixes ({0}) has to be equal to '
                    + 'the number of input file type ({1})'
                    + ''.format(len(suffix), len(input_file_list[1:]))
                )

            # Process inputs
            vm_inst.process(input_file_list[1:], rad, suffix)

        # Run in MULTI-EPOCH mode
        elif mode == 'MULTI-EPOCH':
            # Fetch image directory and patterns
            image_dir = config.getlist('VIGNETMAKER_RUNNER', 'ME_IMAGE_DIR')
            image_pattern = config.getlist(
                'VIGNETMAKER_RUNNER',
                'ME_IMAGE_PATTERN',
            )
            # Fetch WCS log path
            f_wcs_path = config.getexpanded('VIGNETMAKER_RUNNER', 'ME_LOG_WCS')

            # Process inputs
            vm_inst.process_me(image_dir, image_pattern, f_wcs_path, rad)

        # Invalid mode
        else:
            # Raise error for invalid run mode
            raise ValueError('Invalid MODE=\'{}\''.format(mode))

    # No return objects
    return None, None
