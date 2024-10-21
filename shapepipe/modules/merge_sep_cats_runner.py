"""MERGE SEP CATS RUNNER.

Module runner for ``merge_sep_cats``

:Authors: Morgan Schmitz, Axel Guinot,
    Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.merge_sep_cats_package.merge_sep_cats import MergeSep
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    input_module="ngmix_runner",
    file_pattern=["ngmix"],
    file_ext=[".fits"],
    depends=["numpy"],
)
def merge_sep_cats_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Merge SEP Catalogues Runner."""
    # Get config entries
    n_split_max = config.getint(module_config_sec, "N_SPLIT_MAX")

    file_pattern = config.getlist(module_config_sec, "FILE_PATTERN")
    file_ext = config.getlist(module_config_sec, "FILE_EXT")

    if config.has_option(module_config_sec, "WARNING"):
        warning = config.get(module_config_sec, "WARNING")
    else:
        warning = "error"

    # Create merge sep cat class instance
    merge_sep_inst = MergeSep(
        input_file_list,
        file_number_string,
        file_pattern,
        file_ext,
        run_dirs["output"],
        n_split_max,
        warning,
        w_log,
    )

    # Run processing
    merge_sep_inst.process()

    # No return objects
    return None, None
