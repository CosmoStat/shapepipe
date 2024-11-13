"""RANDOM CAT RUNNER.

Module runner for ``random_cat``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.random_cat_package.random_cat import RandomCat


@module_runner(
    version="1.1",
    file_pattern=["image", "pipeline_flag"],
    file_ext=[".fits", "fits"],
    depends=["astropy"],
    numbering_scheme="_0",
)
def random_cat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Random Catalogue Runner."""
    # Get input file names of image and mask
    input_image_name = input_file_list[0]
    input_mask_name = input_file_list[1]

    # Set output file name
    if config.has_option(module_config_sec, "OUTPUT_FILE_PATTERN"):
        output_file_pattern = config.get(
            module_config_sec, "OUTPUT_FILE_PATTERN"
        )
    else:
        output_file_pattern = "random_cat"

    # Get number of random objects requested on output
    n_rand = config.getfloat(module_config_sec, "N_RANDOM")

    # Flag whether n_rand is total (DENSITY=False, default)
    # or per square degree (DENSITY=True)
    if config.has_option(module_config_sec, "DENSITY"):
        density = config.getboolean(module_config_sec, "DENSITY")
    else:
        density = False

    # Get healpix output options
    save_mask_as_healpix = config.getboolean(
        module_config_sec, "SAVE_MASK_AS_HEALPIX"
    )
    if save_mask_as_healpix:
        healpix_options = {}
        for option_trunc in ['FILE_BASE', 'OUT_NSIDE']:
            option = f'HEALPIX_OUT_{option_trunc}'
            healpix_options[option_trunc] = config.get(
                module_config_sec, option
            )
    # Create rand cat class instance
    rand_cat_inst = RandomCat(
        input_image_name,
        input_mask_name,
        run_dirs["output"],
        file_number_string,
        output_file_pattern,
        n_rand,
        density,
        w_log,
        healpix_options,
    )

    # Run processing
    rand_cat_inst.process()

    # No return objects
    return None, None
