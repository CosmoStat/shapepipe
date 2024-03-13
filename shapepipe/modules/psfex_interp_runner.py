"""PSFEX_INTERP RUNNER.

Module runner for ``psfex_interp``.

:Author: Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.psfex_interp_package import psfex_interp

from shapepipe.pipeline.run_log import get_last_dir, get_all_dirs


@module_runner(
    version='1.1',
    input_module=['psfex_runner', 'setools_runner'],
    file_pattern=['star_selection', 'galaxy_selection'],
    file_ext=['.psf', '.fits'],
    depends=['numpy', 'astropy', 'galsim', 'sqlitedict'],
)
def psfex_interp_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The PSFEx Interpolation Runner."""
    # Fetch interpolation run mode
    mode = config.get(module_config_sec, 'MODE')
    # Fetch parameter values
    pos_params = config.getlist(module_config_sec, 'POSITION_PARAMS')
    get_shapes = config.getboolean(module_config_sec, 'GET_SHAPES')
    star_thresh = config.getint(module_config_sec, 'STAR_THRESH')
    chi2_thresh = config.getint(module_config_sec, 'CHI2_THRESH')

    # Run in CLASSIC mode
    if mode == 'CLASSIC':

        # Set input paths
        psfcat_path, galcat_path = input_file_list

        # Create instance of PSFExInterpolator
        pi_inst = psfex_interp.PSFExInterpolator(
            psfcat_path,
            galcat_path,
            run_dirs['output'],
            file_number_string,
            w_log,
            pos_params,
            get_shapes,
            star_thresh,
            chi2_thresh,
        )

        # Process inputs
        pi_inst.process()

    # Run in MULTI-EPOCH mode
    elif mode == 'MULTI-EPOCH':

        # Fetch multi-epoch parameters
        module = config.getexpanded(
            module_config_sec,
            'ME_DOT_PSF_DIR',
        )
        module_name = module.split(":")[-1]
        if "last" in module:
            dot_psf_dirs = [get_last_dir(run_dirs['run_log'], module_name)]
        elif "all" in module:
            dot_psf_dirs = get_all_dirs(run_dirs["run_log"], module_name)
        else:
            raise ValueError(
                "Expected qualifier 'last:' or 'all' before module"
               + f" '{module}' in config entry 'ME_DOT_PSF_DIR'")

        dot_psf_pattern = config.get(
            module_config_sec,
            'ME_DOT_PSF_PATTERN',
        )
        f_wcs_path = config.getexpanded(module_config_sec, 'ME_LOG_WCS')

        # Set input paths
        galcat_path = input_file_list[0]

        # Create instance of PSFExInterpolator
        psfex_interp_inst = psfex_interp.PSFExInterpolator(
            None,
            galcat_path,
            run_dirs['output'],
            file_number_string,
            w_log,
            pos_params,
            get_shapes,
            star_thresh,
            chi2_thresh,
        )

        # Process inputs multi-epoch
        psfex_interp_inst.process_me(dot_psf_dirs, dot_psf_pattern, f_wcs_path)

    # Run in VALIDATION mode
    elif mode == 'VALIDATION':

        # Set input paths
        psfcat_path, galcat_path, psfex_cat_path = input_file_list

        # Create instance of PSFExInterpolator
        psfex_interp_inst = psfex_interp.PSFExInterpolator(
            psfcat_path,
            galcat_path,
            run_dirs['output'],
            file_number_string,
            w_log,
            pos_params,
            get_shapes,
            star_thresh,
            chi2_thresh,
        )

        # Process inputs validation
        psfex_interp_inst.process_validation(psfex_cat_path)

    else:
        # Raise error for invalid run mode
        ValueError('MODE has to be in : [CLASSIC, MULTI-EPOCH, VALIDATION]')

    # No return objects
    return None, None
