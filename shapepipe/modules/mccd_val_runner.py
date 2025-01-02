"""MCCD VAL RUNNER.

This file is the pipeline validation runner for the MCCD package.

:Author: Tobias Liaudat

"""

import mccd

from shapepipe.modules.mccd_package import shapepipe_auxiliary_mccd as aux_mccd
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    input_module=["mccd_preprocessing_runner"],
    file_pattern=["fitted_model", "test_star_selection"],
    file_ext=[".npy", ".fits"],
    numbering_scheme="-0000000",
    depends=["numpy", "mccd", "galsim"],
    run_method="parallel",
)
def mccd_val_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The MCCD Validation Runner."""
    # Recover the MCCD config file and its params
    config_file_path = config.getexpanded(module_config_sec, "CONFIG_PATH")
    mccd_mode = config.get(module_config_sec, "MODE")
    verbose = config.getboolean(module_config_sec, "VERBOSE")

    # Parse MCCD config file
    mccd_parser = mccd.auxiliary_fun.MCCDParamsParser(config_file_path)
    mccd_parser.parse_document()

    # Prepare inputs to run the main fit function
    output_dir = run_dirs["output"] + "/"
    val_saving_name = "validation_psf"

    # Extract the MCCD model path
    mccd_model_path = input_file_list[0]
    # Validation stars are in the second position of the list
    teststar_path = input_file_list[1]

    if mccd_mode == "VALIDATION":

        aux_mccd.mccd_validation_pipeline(
            teststar_path=teststar_path,
            mccd_model_path=mccd_model_path,
            mccd_parser=mccd_parser,
            output_dir=output_dir,
            file_number_string=file_number_string,
            w_log=w_log,
            val_saving_name=val_saving_name,
        )

    else:
        raise ValueError(
            "The mccd_val_runner should only be called when the MODE is "
            + '"VALIDATION".'
        )

    # No return objects
    return None, None
