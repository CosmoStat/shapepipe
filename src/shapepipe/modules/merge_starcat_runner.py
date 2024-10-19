"""MERGE STARCAT RUNNER.

Module runner for ``merge_starcat``.

:Author: Tobias Liaudat, Morgan Schmitz, Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.merge_starcat_package import merge_starcat
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version='1.1',
    input_module=['mccd_fit_val_runner'],
    file_pattern=['validation_psf'],
    file_ext=['.fits'],
    numbering_scheme='-0000000',
    depends=['numpy', 'astropy'],
    run_method='serial',
)
def merge_starcat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Merge Star Catalogues Runner."""
    # Read config file options
    psf_model = config.get(module_config_sec, 'PSF_MODEL')
    allowed_psf_models = ('psfex', 'mccd', 'setools')
    if psf_model not in allowed_psf_models:
        raise ValueError(
            f'Invalid config entry PSF_MODEL={psf_model} found, '
            + f'needs to be one of {allowed_psf_models}'
        )

    # Set output directory
    output_dir = run_dirs['output']

    # Set merge class to use
    if psf_model == 'mccd':
        MSC = merge_starcat.MergeStarCatMCCD
    elif psf_model == 'psfex':
        MSC = merge_starcat.MergeStarCatPSFEX
    elif psf_model == 'setools':
        MSC = merge_starcat.MergeStarCatSetools

    # Create instance of merge class
    merge_inst = MSC(input_file_list, output_dir, w_log)

    # Run processing
    merge_inst.process()

    # No return objects
    return None, None
