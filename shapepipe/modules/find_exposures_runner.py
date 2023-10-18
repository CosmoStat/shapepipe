"""FIND_EXPOSURES RUNNER.

Module runner for ``find_exposures``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.find_exposures_package import find_exposures
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version='1.1',
    file_pattern=['image'],
    file_ext='.fits',
    depends=['numpy', 'astropy', 'sip_tpv'],
    numbering_scheme='_0',
)
def find_exposures_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Find Exposures Runner."""
    # Get input file name of tile
    input_file_name = input_file_list[0]

    # Create output ascii file name
    output_path = f'{run_dirs["output"]}/exp_numbers{file_number_string}.txt'

    # Give column number for exposure name in fits header
    colnum = config.getint(module_config_sec, 'COLNUM')

    # Give the prefix of exposures
    exp_prefix = config.get(module_config_sec, 'EXP_PREFIX')
    # Create find exposures class instance
    find_exp_inst = find_exposures.FindExposures(
        input_file_name,
        output_path,
        w_log,
        colnum,
        exp_prefix,
    )

    # Run processing
    find_exp_inst.process()

    # No return objects
    return None, None
