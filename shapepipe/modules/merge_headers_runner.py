"""MERGE HEADERS RUNNER.

Module runner for ``merge_headers``.

:Author: Axel Guinot

"""

from shapepipe.modules.merge_headers_package.merge_headers import merge_headers
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    input_module="split_exp_runner",
    version="1.1",
    file_pattern=["headers"],
    file_ext=[".npy"],
    depends=["numpy", "sqlitedict"],
    run_method="serial",
)
def merge_headers_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Merge Headers Runner."""
    # Set output directory
    output_dir = run_dirs["output"]
    if config.has_option(module_config_sec, "OUTPUT_PATH"):
        output_dir = config.getexpanded(module_config_sec, "OUTPUT_PATH")

    # Log output directory
    w_log.info(f"output_dir = {output_dir}")

    # Merge header files
    merge_headers(input_file_list, output_dir)

    w_log.info(f"Merged {len(input_file_list)} input file headers")

    # No return objects
    return None, None
