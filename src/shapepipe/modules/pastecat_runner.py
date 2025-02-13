"""PASTECAT RUNNER.

Module runner for ``pastecat``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.pastecat_package.pastecat import PasteCat


@module_runner(
    version="1.1",
    input_module="sextractor_runner",
    file_pattern="tile_sexcat",
    file_ext=".fits",
    depends=["numpy", "astropy"],
    run_method="parallel",
)
def paste_cat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Paste Catalogue Runner."""
    # Get config options
    if config.has_option(module_config_sec, "CHECK_COL_NAME"):
        check_col_name = config.get(module_config_sec, "CHECK_COL_NAME")
    else:
        check_col_name = None

    if config.has_option(module_config_sec, "HDU"):
        tmp = config.getlist(module_config_sec, "HDU")
        hdu_no = [int(idx) for idx in tmp]
        if len(hdu_no) != len(input_file_list):
            raise IndexError(
                "Different lengths for input file list "
                + f"({len(input_file_list)}) and HDU ({len(hdu_no)})."
            )
    else:
        hdu_no = None

    if config.has_option(module_config_sec, "PREFIX"):
        prefix = config.get(module_config_sec, "PREFIX")
    else:
        prefix = "cat_pasted"

    if config.has_option(module_config_sec, "EXT_NAME"):
        ext_name_list = config.getlist(module_config_sec, "EXT_NAME")
        if len(ext_name_list) != len(input_file_list):
            raise ValueError(
                f"Input file list length ({len(input_file_list)}) "
                + f"and EXT_NAME list ({len(ext_name_list)})"
                + "need to be equal."
            )
    else:
        ext_name_list = None

    # Set file extension
    file_ext = "fits"

    # Set output path
    output_path = (
        f'{run_dirs["output"]}/{prefix}' + f"{file_number_string}.{file_ext}"
    )

    # Create rand cat class instance
    paste_cat_inst = PasteCat(
        input_file_list,
        output_path,
        w_log,
        ext_name=ext_name_list,
        check_col_name=check_col_name,
        hdu_no=hdu_no,
    )

    # Run processing
    paste_cat_inst.process()

    # No return objects
    return None, None
