"""MAKE CATALOGUE RUNNER.

Module runner for ``make_cat``.

:Authors: Axel Guinot, Martin Kilbinger

"""

import os

from shapepipe.modules.make_cat_package import make_cat
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    input_module=[
        "sextractor_runner",
        "spread_model_runner",
        "psfex_interp_runner",
        "ngmix_runner",
    ],
    file_pattern=[
        "tile_sexcat",
        "sexcat",
        "galaxy_psf",
        "ngmix",
    ],
    file_ext=[".fits", ".fits", ".sqlite", ".fits"],
    depends=["numpy", "sqlitedict"],
)
def make_cat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Make Catalogue Runner."""
    # Set input file paths
    if len(input_file_list) == 3:
        # No spread model input
        (
            tile_sexcat_path,
            galaxy_psf_path,
            shape1_cat_path,
        ) = input_file_list
    else:
        raise IndexError(
            f"Invalid number of input files {len(input_file_List)}, expected 3."
        )

    # Fetch shape measurement type
    shape_type_list = config.getlist(
        module_config_sec,
        "SHAPE_MEASUREMENT_TYPE",
    )
    for shape_type in shape_type_list:
        if shape_type.lower() not in ["ngmix", "galsim"]:
            raise ValueError(
                "SHAPE_MEASUREMENT_TYPE must be in [ngmix, galsim]"
            )

    # Fetch PSF data option
    if config.has_option(module_config_sec, "SAVE_PSF_DATA"):
        save_psf = config.getboolean(module_config_sec, "SAVE_PSF_DATA")
    else:
        save_psf = False

    # Set final output file
    final_cat_file = make_cat.prepare_final_cat_file(
        run_dirs["output"],
        file_number_string,
    )

    # Save SExtractor data
    w_log.info("Save SExtractor data")
    n_obj = make_cat.save_sextractor_data(final_cat_file, tile_sexcat_path)
    cat_size_sextractor = n_obj

    # Save shape data
    sc_inst = make_cat.SaveCatalogue(final_cat_file, cat_size_sextractor, w_log)
    w_log.info("Save shape measurement data")
    for shape_type in shape_type_list:
        w_log.info(f"Save {shape_type.lower()} data")
        cat_path = (
            shape2_cat_path if shape_type == "galsim" else shape1_cat_path
        )
        err_msg = sc_inst.process(shape_type.lower(), cat_path)


        # If error message: delete (incomplete) output file and raise error
        if err_msg is not None:
            os.remove(
                make_cat.get_output_name(
                    run_dirs["output"],
                    file_number_string,
                )
            )
            #raise ValueError(err_msg)
            w_log.info(err_msg)

    if save_psf:
        err_msg = sc_inst.process("psf", galaxy_psf_path)

    return None, None
