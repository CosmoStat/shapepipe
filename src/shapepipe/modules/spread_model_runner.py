"""SPREAD MODEL RUNNER.

Module runner for ``spread_model``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.spread_model_package.spread_model import SpreadModel


@module_runner(
    version="1.1",
    input_module=[
        "sextractor_runner",
        "psfex_interp_runner",
        "vignetmaker_runner",
    ],
    file_pattern=["sexcat", "galaxy_psf", "weight_vign"],
    file_ext=[".fits", ".sqlite", ".fits"],
    depends=["numpy", "galsim"],
    run_method="parallel",
)
def spread_model_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Spread Model Runner."""
    # Get input files
    sex_cat_path, psf_cat_path, weight_cat_path = input_file_list

    # Get file prefix (optional)
    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(module_config_sec, "PREFIX")
        if (prefix.lower() != "none") & (prefix != ""):
            prefix = prefix + "_"
        else:
            prefix = ""
    else:
        prefix = ""

    # Get pixel scale and output mode
    pixel_scale = config.getfloat(module_config_sec, "PIXEL_SCALE")
    output_mode = config.get(module_config_sec, "OUTPUT_MODE")

    # Set output file path
    file_name = f"{prefix}sexcat_sm{file_number_string}.fits"
    output_path = f'{run_dirs["output"]}/{file_name}'

    # Create spread model class instance
    sm_inst = SpreadModel(
        sex_cat_path,
        psf_cat_path,
        weight_cat_path,
        output_path,
        pixel_scale,
        output_mode,
    )

    # Process spread model computation
    sm_inst.process()

    # No return objects
    return None, None
