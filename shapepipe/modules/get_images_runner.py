# -*- coding: utf-8 -*-

"""GET IMAGES RUNNER

Module runner for ``get_images``

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2019 - 2021

:Package: ShapePipe

"""


from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.get_images_package import get_images as gi


@module_runner(
    version='1.0',
    depends=['numpy'],
    run_method='serial',
    numbering_scheme='_0')
def get_images_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    w_log
):

    # Read config file section

    # Copy/download method
    retrieve_method = config.get('GET_IMAGES_RUNNER', 'RETRIEVE')
    retrieve_ok = ['vos', 'symlink']
    if retrieve_method not in retrieve_ok:
        raise ValueError(
            f'key RETRIEVE={retrieve_method} is invalid, must be in {retrieve_ok}'
        )

    if config.has_option('GET_IMAGES_RUNNER', 'RETRIEVE_OPTIONS'):
        retrieve_options = config.getexpanded(
            'GET_IMAGES_RUNNER',
            'RETRIEVE_OPTIONS'
        )
    else:
        retrieve_options = None

    # Paths
    input_dir = config.getlist('GET_IMAGES_RUNNER', 'INPUT_PATH')
    input_file_pattern = config.getlist(
        'GET_IMAGES_RUNNER',
        'INPUT_FILE_PATTERN'
    )
    input_file_ext = config.getlist(
        'GET_IMAGES_RUNNER',
        'INPUT_FILE_EXT'
    )
    output_file_pattern = config.getlist(
        'GET_IMAGES_RUNNER',
        'OUTPUT_FILE_PATTERN'
    )

    input_numbering = config.get('GET_IMAGES_RUNNER', 'INPUT_NUMBERING')

    if config.has_option('GET_IMAGES_RUNNER', 'OUTPUT_PATH'):
        output_dir = config.getexpanded('GET_IMAGES_RUNNER', 'OUTPUT_PATH')
    else:
        output_dir = run_dirs['output']

    # Flags to check for already retrieved files
    if config.has_option('GET_IMAGES_RUNNER', 'CHECK_EXISTING_DIR'):
        check_existing_dir = config.getexpanded(
            'GET_IMAGES_RUNNER',
            'CHECK_EXISTING_DIR'
        )
        if config.has_option('GET_IMAGES_RUNNER', 'N_EXPECTED'):
            n_expected = config.getint('GET_IMAGES_RUNNER', 'N_EXPECTED')
        else:
            n_expected = 1
    else:
        check_existing_dir = None
        n_expected = None

    inst = gi.GetImages(
        retrieve_method,
        retrieve_options,
        input_file_list,
        input_numbering,
        input_file_pattern,
        input_file_ext,
        output_file_pattern,
        w_log,
        check_existing_dir=check_existing_dir,
        n_expected=n_expected
    )

    inst.process(input_dir, output_dir)

    return None, None
