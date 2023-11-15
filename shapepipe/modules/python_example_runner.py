"""PYTHON MODULE EXAMPLE.

This module defines methods for an example Python module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.python_example_package import python_example


@module_runner(
    version='1.1',
    file_pattern=['numbers', 'letters'],
    file_ext='.txt',
    depends=[
        'numpy',
        'astropy',
        'galsim',
        'joblib',
        'ngmix',
        'pandas',
        'scipy',
        'sf_tools',
        'sip_tpv',
        'sqlitedict',
        'treecorr',
    ],
    run_method='parallel',
)
def python_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Python Example Runner."""
    # Set output file name
    output_file_name = (
        f'{run_dirs["output"]}/pyex_output{file_number_string}.cat'
    )

    # Retrieve log message from config file
    message = config.get(module_config_sec, 'MESSAGE')

    # Create an instance of the Python example class
    py_ex_inst = python_example.PythonExample(0)

    # Read input files
    py_ex_inst.read_files(*input_file_list)

    # Write output files
    py_ex_inst.write_file(output_file_name, message)

    # Return file content and no stderr
    return py_ex_inst.content, None
