# -*- coding: utf-8 -*-

"""MATCH EXTERNAL RUNNER

This module matches an external catalogue to a ShapePipe (SExtractor) catalog

:Author: Martin Kilbinger, Xavier Jimenez

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.match_external_package import match_external


@module_runner(
    version='1.1',
    input_module='sextractor_runner',
    file_pattern='tile_sexcat',
    file_ext='.fits',
    depends=['numpy', 'astropy'],
    run_method='parallel',
)
def match_external_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):

    # Get processing tolerance
    tolerance = config.getfloat(module_config_sec, 'TOLERANCE')

    # Internal data
    col_match = config.getlist(module_config_sec, 'COL_MATCH')
    if config.has_option(module_config_sec, 'HDU'):
        hdu_no = config.getint(module_config_sec, 'HDU')
    else:
        hdu_no = 2

    # Set run mode
    mode = config.get(module_config_sec, 'MODE')
    valid_modes = ('CLASSIC', 'MULTI-EPOCH')
    if mode not in valid_modes:
        raise ValueError(
            f'mode \'{mode}\' is invalid, must be one of {valid_modes}.'
        )

    # External data
    external_cat_path = config.getexpanded(
        module_config_sec,
        'EXTERNAL_CAT_PATH',
    )
    external_col_match = config.getlist(
        module_config_sec,
        'EXTERNAL_COL_MATCH',
    )

    # TODO: optional or 'none', 'all'
    # Also TODO: change column name if already present in internal cat
    external_col_copy = config.getlist(module_config_sec, 'EXTERNAL_COL_COPY')

    if config.has_option(module_config_sec, 'EXTERNAL_HDU'):
        external_hdu_no = config.getint(module_config_sec, 'EXTERNAL_HDU')
    else:
        external_hdu_no = 1

    # Output
    if config.has_option(module_config_sec, 'OUTPUT_FILE_PATTERN'):
        output_file_pattern = config.get(
            module_config_sec,
            'OUTPUT_FILE_PATTERN',
        )
    else:
        output_file_pattern = 'cat_matched'

    if config.has_option(module_config_sec, 'MARK_NON_MATCHED'):
        mark_non_matched = config.getfloat(
            module_config_sec,
            'MARK_NON_MATCHED',
        )
    else:
        mark_non_matched = None

    # Set output file path
    file_ext = 'fits'
    output_path = (
        f'{run_dirs['output']}/{output_file_pattern}{file_number_string}.'
        + f'{file_ext}'
    )

    # Create instance of MatchCats
    mc_inst = MatchCats(
        input_file_list,
        output_path,
        w_log,
        tolerance,
        col_match,
        hdu_no,
        mode,
        external_cat_path,
        external_col_match,
        external_col_copy,
        external_hdu_no=external_hdu_no,
        mark_non_matched=mark_non_matched,
    )

    # Process inputs
    mc_inst.process()

    # No return objects
    return None, None
