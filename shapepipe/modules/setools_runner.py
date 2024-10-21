"""SETOOLS RUNNER.

Module runner for ``setools``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.setools_package.setools import SETools


@module_runner(
    input_module="sextractor_runner",
    version="1.1",
    file_pattern=["sexcat"],
    file_ext=[".fits"],
    depends=["numpy", "matplotlib"],
)
def setools_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The SETools Runner."""
    # Get path to setools configuration file
    config_file = config.getexpanded(module_config_sec, "SETOOLS_CONFIG_PATH")

    # Create instance of SETools
    se_inst = SETools(
        input_file_list[0],
        run_dirs["output"],
        file_number_string,
        config_file,
    )

    # Process inputs
    se_inst.process(w_log)

    # No return objects
    return None, None
