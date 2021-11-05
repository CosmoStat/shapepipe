# -*- coding: utf-8 -*-

"""PASTE CAT RUNNER

Pipeline runner for the PasteCat package.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Axel Guinot

:Date: 10/2020

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.pastecat_package.pastecat_script import PasteCat


@module_runner(
    version='1.1',
    input_module='sextractor_runner',
    file_pattern='tile_sexcat',
    file_ext='.fits',
    depends=['numpy', 'astropy'],
    run_method='parallel'
)
def paste_cat_runner(input_file_list, run_dirs, file_number_string,
                     config, module_config_sec, w_log):

    if config.has_option('PASTE_CAT_RUNNER', 'CHECK_COL_NAME'):
        check_col_name = config.get('PASTE_CAT_RUNNER', 'CHECK_COL_NAME')
    else:
        check_col_name = None

    if config.has_option('PASTE_CAT_RUNNER', 'HDU'):
        tmp = config.getlist('PASTE_CAT_RUNNER', 'HDU')
        hdu_no = [int(i) for i in tmp]
        if len(hdu_no) != len(input_file_list):
            raise IndexError(
                f'Different lengths for input file list'
                +' ({len(input_file_list)}) and HDU ({len(hdu_no)})'
            )
    else:
        hdu_no = None

    if config.has_option('PASTE_CAT_RUNNER', 'OUTPUT_FILE_PATTERN'):
        output_file_pattern = config.get(
            'PASTE_CAT_RUNNER',
            'OUTPUT_FILE_PATTERN'
        )
    else:
        output_file_pattern = 'cat_pasted'

    if config.has_option('PASTE_CAT_RUNNER', 'EXT_NAME'):
        ext_name_list = config.getlist('PASTE_CAT_RUNNER', 'EXT_NAME')
        if len(ext_name_list) != len(input_file_list):
            raise ValueError(
        f'Input file list length ({len(input_file_list)})'
        +' and EXT_NAME list ({len(ext_name_list)})'
        + 'need to be equal'
            )
    else:
        ext_name_list = None

    file_ext = 'fits'

    output_path = '{}/{}{}.{}'.format(
        run_dirs['output'],
        output_file_pattern,
        file_number_string,
        file_ext
    )

    inst = PasteCat(
        input_file_list,
        output_path,
        w_log,
        ext_name = ext_name_list,
        check_col_name = check_col_name,
        hdu_no=hdu_no
    )

    inst.process()

    return None, None
