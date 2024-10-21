"""SPLIT EXP RUNNER.

Module runner for ``split_exp``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.split_exp_package.split_exp import SplitExposures


@module_runner(
    version="1.1",
    input_module="get_images_runner",
    file_pattern=["image", "weight", "flag"],
    file_ext=[".fz", ".fz", ".fz"],
    depends=["numpy", "astropy", "sip_tpv"],
    run_method="parallel",
)
def split_exp_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Split Exposures Runner."""
    # Get file suffix
    output_suffix = config.getlist(module_config_sec, "OUTPUT_SUFFIX")
    # Get the number of HDUs
    n_hdu = config.getint(module_config_sec, "N_HDU")

    # Create split exposures class instance
    split_inst = SplitExposures(
        input_file_list,
        run_dirs["output"],
        file_number_string,
        output_suffix,
        n_hdu,
    )

    # Process splitting
    split_inst.process()

    # No return objects
    return None, None
