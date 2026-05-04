"""MERGE HEADERS RUNNER.

Module runner for ``merge_headers``.

:Authors: Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.merge_headers_package.merge_headers import merge_headers
from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline.exp_utils import get_exp_output_files


@module_runner(
    input_module=["split_exp_runner", "find_exposures_runner"],
    version="1.2",
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
    output_dir = run_dirs["output"]
    w_log.info(f"output_dir = {output_dir}")

    if config.has_option(module_config_sec, "EXP_BASE_DIR"):
        # Tile-level mode: input is an exp_numbers txt file; collect header
        # files from each per-exposure work directory via get_exp_output_files.
        exp_base_dir = config.getexpanded(module_config_sec, "EXP_BASE_DIR")
        exp_numbers_file = input_file_list[0][0]
        w_log.info(
            f"Tile-level merge: collecting headers from {exp_base_dir} "
            f"using {exp_numbers_file}"
        )
        headers_file_list = get_exp_output_files(
            exp_base_dir,
            exp_numbers_file,
            "split_exp_runner",
            "headers",
            ".npy",
            w_log=w_log,
        )
        merge_headers(headers_file_list, output_dir, file_number_string)
        w_log.info(f"Merged {len(headers_file_list)} exposure header files")
    else:
        # Per-exposure mode: input_file_list already contains the header files.
        merge_headers(input_file_list, output_dir, file_number_string)
        w_log.info(f"Merged {len(input_file_list)} input file headers")

    # No return objects
    return None, None
