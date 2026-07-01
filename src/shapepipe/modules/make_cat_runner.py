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
        "sexcat_sm",
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
        sexcat_sm_path = None
    else:
        # With spread model input
        (
            tile_sexcat_path,
            sexcat_sm_path,
            galaxy_psf_path,
            shape1_cat_path,
        ) = input_file_list[0:4]
        if len(input_file_list) == 5:
            # With second shape catalogue input
            shape2_cat_path = input_file_list[4]

    # Fetch classification options
    do_classif = config.getboolean(
        module_config_sec,
        "SM_DO_CLASSIFICATION",
    )
    if do_classif:
        star_thresh = config.getfloat(module_config_sec, "SM_STAR_THRESH")
        gal_thresh = config.getfloat(module_config_sec, "SM_GAL_THRESH")
    else:
        star_thresh = None
        gal_thresh = None

    # Fetch shape measurement type
    shape_type_list = config.getlist(
        module_config_sec,
        "SHAPE_MEASUREMENT_TYPE",
    )
    for shape_type in shape_type_list:
        if shape_type.lower() != "ngmix":
            raise ValueError("SHAPE_MEASUREMENT_TYPE must be ngmix")

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

    # Save spread-model data
    if sexcat_sm_path is None:
        w_log.info("No sm cat input, setting spread model to 99")
    else:
        w_log.info("Save spread-model data")
        cat_size_sm = make_cat.save_sm_data(
            final_cat_file,
            sexcat_sm_path,
            do_classif,
            star_thresh,
            gal_thresh,
            n_obj=n_obj
        )

        if cat_size_sextractor != cat_size_sm:
            w_log(
                f"Warnign: SExtractor catalogue {tile_sexcat_path} has different size"
                + f" ({cat_size_sextractor} than spread_model catalogue"
                + f" {sexcat_sm_path} ({cat_size_sm})"
            )

    # Save shape data
    sc_inst = make_cat.SaveCatalogue(final_cat_file, cat_size_sextractor, w_log)
    w_log.info("Save shape measurement data")
    for shape_type in shape_type_list:
        w_log.info(f"Save {shape_type.lower()} data")
        err_msg = sc_inst.process(shape_type.lower(), shape1_cat_path)

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
