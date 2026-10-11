"""MERGE HEADERS RUNNER.

Module runner for ``merge_headers``.

:Authors: Axel Guinot, Martin Kilbinger

"""

import os
import re

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

    # The input is an exp_numbers txt file (find_exposures_runner); the header
    # files themselves are collected from each per-exposure work directory
    # under EXP_BASE_DIR via get_exp_output_files.
    exp_base_dir = config.getexpanded(module_config_sec, "EXP_BASE_DIR")

    # In serial mode several tiles' exp_numbers files can arrive in a
    # single call; merge each tile into its own per-tile sqlite file.
    for exp_numbers_file in (item[0] for item in input_file_list):
        w_log.info(
            f"Collecting headers from {exp_base_dir} "
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
        # Extract tile number from the exp_numbers filename, e.g.
        # "exp_numbers-284.272-1.000.txt" -> "-284.272-1.000"
        base = os.path.splitext(os.path.basename(exp_numbers_file))[0]
        tile_number = re.sub(r"^exp_numbers", "", base)
        merge_headers(headers_file_list, output_dir, tile_number)
        w_log.info(f"Merged {len(headers_file_list)} exposure header files")

    # No return objects
    return None, None
