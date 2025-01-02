"""MATCH EXTERNAL RUNNER.

Module runner for ``match_external``.

:Author: Martin Kilbinger, Xavier Jimenez

"""

from shapepipe.modules.match_external_package.match_external import MatchCats
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    input_module="sextractor_runner",
    file_pattern="tile_sexcat",
    file_ext=".fits",
    depends=["numpy", "astropy"],
    run_method="parallel",
)
def match_external_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Match External Runner."""
    # Get processing tolerance
    tolerance = config.getfloat(module_config_sec, "TOLERANCE")

    # Internal data
    col_match = config.getlist(module_config_sec, "COL_MATCH")
    if config.has_option(module_config_sec, "HDU"):
        hdu_no = config.getint(module_config_sec, "HDU")
    else:
        hdu_no = 2

    # Set run mode
    mode = config.get(module_config_sec, "MODE")
    valid_modes = ("CLASSIC", "MULTI-EPOCH")
    if mode not in valid_modes:
        raise ValueError(
            f"mode '{mode}' is invalid, must be one of {valid_modes}."
        )

    # External data
    external_cat_path = config.getexpanded(
        module_config_sec,
        "EXTERNAL_CAT_PATH",
    )
    external_col_match = config.getlist(
        module_config_sec,
        "EXTERNAL_COL_MATCH",
    )

    # TODO: optional or 'none', 'all'
    # Also TODO: change column name if already present in internal cat
    external_col_copy = config.getlist(module_config_sec, "EXTERNAL_COL_COPY")

    if config.has_option(module_config_sec, "EXTERNAL_HDU"):
        external_hdu_no = config.getint(module_config_sec, "EXTERNAL_HDU")
    else:
        external_hdu_no = 1

    # Output
    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(
            module_config_sec,
            "PREFIX",
        )
    else:
        prefix = "cat_matched"

    if config.has_option(module_config_sec, "MARK_NON_MATCHED"):
        mark_non_matched = config.getfloat(
            module_config_sec,
            "MARK_NON_MATCHED",
        )
    else:
        mark_non_matched = None

    if config.has_option(module_config_sec, "OUTPUT_DISTANCE"):
        output_distance = config.getboolean(
            module_config_sec,
            "OUTPUT_DISTANCE",
        )
    else:
        output_distance = False

    # Set output file path
    file_ext = "fits"
    output_path = (
        f'{run_dirs["output"]}/{prefix}{file_number_string}.' + f"{file_ext}"
    )

    # Create instance of MatchCats
    match_cats_inst = MatchCats(
        input_file_list,
        output_path,
        w_log,
        tolerance,
        col_match,
        hdu_no,
        mode,
        external_cat_path,
        external_col_match,
        external_col_copy,
        external_hdu_no=external_hdu_no,
        mark_non_matched=mark_non_matched,
        output_distance=output_distance,
    )

    # Process inputs
    match_cats_inst.process()

    # No return objects
    return None, None
