# -*- coding: utf-8 -*-

"""PYTHON MODULE EXAMPLE

This module defines methods for an example Python module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.python_example_package import python_example as pe


@module_runner(
    version='1.1',
    file_pattern=['numbers', 'letters'],
    file_ext='.txt',
    depends=[
        'numpy',
        'astropy',
        'galsim',
        'joblib',
        'mccd',
        'ngmix',
        'pandas',
        'pysap',
        'scipy',
        'sf_tools',
        'sip_tpv',
        'sqlitedict',
        'treecorr',
    ],
    run_method='parallel'
)
def python_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):

    output_file_name = (
        f'{run_dirs["output"]}/pyex_output{file_number_string}.cat'
    )
    message = config.get(module_config_sec, 'MESSAGE')

    inst = pe.PythonExample(0)
    inst.read_files(*input_file_list)
    inst.write_file(output_file_name, message)

    return inst.content, None
