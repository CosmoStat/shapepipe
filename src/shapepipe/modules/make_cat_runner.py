"""MAKE CATALOGUE RUNNER.

Module runner for ``make_cat``.

:Authors: Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.make_cat_package import make_cat
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="2.0",
    input_module=[
        "sextractor_runner",
        "psfex_interp_runner",
        "ngmix_runner",
    ],
    file_pattern=[
        "tile_sexcat",
        "galaxy_psf",
        "ngmix",
    ],
    file_ext=[".fits", ".sqlite", ".fits"],
    depends=["numpy", "h5py", "sqlitedict"],
)
def make_cat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Make Catalogue Runner.

    @sc [decision:catalogue_assembly.star_galaxy_classification]

    @sc [decision:masking.sky_mask_application]
    """
    tile_sexcat_path, galaxy_psf_path, shape1_cat_path = input_file_list

    # Fetch shape measurement type
    shape_type_list = config.getlist(
        module_config_sec,
        "SHAPE_MEASUREMENT_TYPE",
    )
    for shape_type in shape_type_list:
        if shape_type.lower() != "ngmix":
            raise ValueError("SHAPE_MEASUREMENT_TYPE must be [ngmix]")

    # Fetch PSF data option
    if config.has_option(module_config_sec, "SAVE_PSF_DATA"):
        save_psf = config.getboolean(module_config_sec, "SAVE_PSF_DATA")
    else:
        save_psf = False
    if config.has_option(module_config_sec, "N_EPOCH_SLOTS"):
        n_epoch_slots = config.getint(module_config_sec, "N_EPOCH_SLOTS")
    else:
        n_epoch_slots = None

    # The catalogue is assembled in memory, column by column, and written once.
    w_log.info("Save SExtractor data")
    final_cat = make_cat.read_sextractor_data(tile_sexcat_path)

    sc_inst = make_cat.SaveCatalogue(
        final_cat, len(final_cat["NUMBER"]), w_log
    )
    w_log.info("Save shape measurement data")
    for shape_type in shape_type_list:
        w_log.info(f"Save {shape_type.lower()} data")
        err_msg = sc_inst.process(shape_type.lower(), shape1_cat_path)
        # An incomplete catalogue is never written.
        if err_msg is not None:
            w_log.info(err_msg)
            return None, None

    if save_psf:
        sc_inst.process("psf", galaxy_psf_path, n_epoch_slots=n_epoch_slots)

    # Optional per-band external healsparse mask lookup (UNIONS-WL/spherex#38):
    # add one MASK_<BAND> column per band, queried at each object's world
    # position. Absent config is a strict no-op.
    if config.has_option(module_config_sec, "MASK_EXT_PATHS"):
        band_paths = make_cat.parse_mask_ext_paths(
            config.getexpanded(module_config_sec, "MASK_EXT_PATHS")
        )
        w_log.info("Save external mask data")
        make_cat.save_mask_ext_data(final_cat, band_paths, w_log)

    make_cat.write_final_cat(
        make_cat.get_output_name(run_dirs["output"], file_number_string),
        final_cat,
    )

    return None, None
