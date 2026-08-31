"""MASK_QUERY RUNNER.

Module runner for ``mask_query``.

:Author: Claude Fable 5, for PR #847

"""

from shapepipe.modules.mask_query_package.mask_query import MaskQuery
from shapepipe.modules.module_decorator import module_runner
from shapepipe.utilities import mask_query as mask_query_util


@module_runner(
    version="1.0",
    input_module="sextractor_runner",
    file_pattern=["sexcat"],
    file_ext=[".fits"],
    depends=["numpy", "healsparse"],
)
def mask_query_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Mask Query Runner."""
    sexcat_path = input_file_list[0]

    # Get file prefix (optional)
    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(module_config_sec, "PREFIX")
        if (prefix.lower() != "none") & (prefix != ""):
            prefix = prefix + "_"
        else:
            prefix = ""
    else:
        prefix = ""

    mask_paths = mask_query_util.parse_map_paths(
        config.getexpanded(module_config_sec, "MASK_PATHS")
    )
    if not mask_paths:
        raise ValueError(
            f"[{module_config_sec}] MASK_PATHS is empty; the module has"
            + " nothing to query."
        )

    # Any nonzero map value flags unless a bit selection is given
    if config.has_option(module_config_sec, "MASK_BITS"):
        bits = config.getint(module_config_sec, "MASK_BITS")
    else:
        bits = None

    output_path = (
        f'{run_dirs["output"]}/{prefix}sexcat_ext{file_number_string}.fits'
    )

    mq_inst = MaskQuery(
        sexcat_path,
        output_path,
        mask_paths,
        bits=bits,
        w_log=w_log,
    )
    n_flagged = mq_inst.process()

    w_log.info(f"FLAG_EXT nonzero for {n_flagged} objects")

    # No return objects
    return None, None
