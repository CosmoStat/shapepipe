"""MCCD FIT RUNNER.

This file is the pipeline fit runner for the MCCD package.

:Author: Tobias Liaudat

"""

import mccd

from shapepipe.modules.mccd_package import shapepipe_auxiliary_mccd as aux_mccd
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    input_module=["mccd_preprocessing_runner"],
    file_pattern=["train_star_selection"],
    file_ext=[".fits"],
    numbering_scheme="-0000000",
    depends=["numpy", "mccd", "galsim"],
    run_method="parallel",
)
def mccd_fit_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The MCCD Fit Runner."""
    # Recover the MCCD config file and its params
    config_file_path = config.getexpanded(module_config_sec, "CONFIG_PATH")
    mccd_mode = config.get(module_config_sec, "MODE")
    verbose = config.getboolean(module_config_sec, "VERBOSE")

    # Parse MCCD config file
    mccd_parser = mccd.auxiliary_fun.MCCDParamsParser(config_file_path)
    mccd_parser.parse_document()

    # Prepare inputs to run the main fit function
    trainstar_path = input_file_list[0]
    output_dir = run_dirs["output"] + "/"
    saving_name = "fitted_model"

    if mccd_mode == "FIT":
        aux_mccd.mccd_fit_pipeline(
            trainstar_path=trainstar_path,
            file_number_string=file_number_string,
            mccd_parser=mccd_parser,
            output_dir=output_dir,
            verbose=verbose,
            saving_name=saving_name,
            w_log=w_log,
        )

    else:
        raise ValueError(
            'mccd_fit_runner should be called when the MODE is "FIT".'
        )

    # No return objects
    return None, None
