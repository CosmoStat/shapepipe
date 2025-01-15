"""UNCOMPRESS FITS RUNNER.

Module runner for ``uncompress_fits``.

:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.uncompress_fits_package import uncompress_fits


@module_runner(
    version="1.1",
    file_pattern=["image"],
    file_ext=[".fits"],
    numbering_scheme="_0",
)
def uncompress_fits_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Uncompress Fits Runner."""
    # Get HDU number of input image data
    if config.has_option(module_config_sec, "HDU_DATA"):
        data_hdu = config.getint(module_config_sec, "HDU_DATA")
    else:
        data_hdu = 0

    # Get output patterns
    output_pattern_list = config.getlist(module_config_sec, "OUTPUT_PATTERN")

    # Check consistency of input and output list lengths
    if len(input_file_list) != len(output_pattern_list):
        raise ValueError(
            f"Lists INPUT_PATH ({len(input_file_list)})"
            + f" and OUTPUT_PATTERN ({len(output_pattern_list)})"
            + "need to be of equal length."
        )

    # Create instance of uncompress
    uncompress_inst = uncompress_fits.Uncompress(
        input_file_list,
        output_pattern_list,
        run_dirs["output"],
        file_number_string,
        data_hdu,
    )

    # Process uncompressing
    uncompress_inst.process()

    # No return objects
    return None, None
