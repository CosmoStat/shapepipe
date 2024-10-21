"""MAKE CATALOGUE RUNNER.

Module runner for ``make_cat``.

:Author: Axel Guinot

"""

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
    (
        tile_sexcat_path,
        sexcat_sm_path,
        galaxy_psf_path,
        shape1_cat_path,
    ) = input_file_list[0:4]
    if len(input_file_list) == 5:
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
    make_cat.save_sextractor_data(final_cat_file, tile_sexcat_path)

    # Save spread-model data
    w_log.info("Save spread-model data")
    make_cat.save_sm_data(
        final_cat_file,
        sexcat_sm_path,
        do_classif,
        star_thresh,
        gal_thresh,
    )

    # Save shape data
    sc_inst = make_cat.SaveCatalogue(final_cat_file)
    w_log.info("Save shape measurement data")
    for shape_type in shape_type_list:
        w_log.info(f"Save {shape_type.lower()} data")
        cat_path = (
            shape2_cat_path if shape_type == "galsim" else shape1_cat_path
        )
        sc_inst.process(shape_type.lower(), cat_path)
    if save_psf:
        sc_inst.process("psf", galaxy_psf_path)

    # No return objects
    return None, None
