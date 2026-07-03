"""NGMIX RUNNER.

Module runner for ``ngmix``.

:Author: Axel Guinot

"""

from sqlitedict import SqliteDict

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.ngmix_package.ngmix import Ngmix


@module_runner(
    version="0.0.1",
    input_module=[
        "sextractor_runner",
        "psfex_interp_runner",
        "vignetmaker_runner",
        "merge_headers_runner",
    ],
    file_pattern=[
        "tile_sexcat",
        "image",
        "exp_background",
        "galaxy_psf",
        "weight",
        "flag",
        "log_exp_headers",
    ],
    file_ext=[".fits", ".sqlite", ".sqlite", ".sqlite", ".sqlite", ".sqlite", ".sqlite"],
    depends=["numpy", "ngmix", "galsim", "sqlitedict", "astropy"],
)
def ngmix_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Ngmix Runner."""
    # Read config file entries

    # Photometric zero point
    zero_point = config.getfloat(module_config_sec, "MAG_ZP")

    # Pixel scale
    pixel_scale = config.getfloat(module_config_sec, "PIXEL_SCALE")

    # Input directory to check for already retrieved files
    if config.has_option(module_config_sec, "CHECK_EXISTING_DIR"):
        check_existing_dir = config.getexpanded(
            module_config_sec,
            "CHECK_EXISTING_DIR",
        )
    else:
        check_existing_dir = None

    # Batch save option
    if config.has_option(module_config_sec, "SAVE_BATCH"):
        save_batch = config.getint(
            module_config_sec,
            "SAVE_BATCH",
        )
    else:
        # No batch saving
        save_batch = -1

    # First and last galaxy ID to process
    id_obj_min = config.getint(module_config_sec, "ID_OBJ_MIN")
    id_obj_max = config.getint(module_config_sec, "ID_OBJ_MAX")

    # Background subtraction (disable for image sims where background is absent)
    bkg_sub = config.getboolean(module_config_sec, "BKG_SUB", fallback=True)

    # wcs path and vignet slice depend on whether background vignet is present
    wcs_idx = 6 if bkg_sub else 5
    f_wcs_path = input_file_list[wcs_idx]

    # Check PSF vignets first: if all are empty dicts {}, the exposures for this
    # tile are absent from the PSF dictionary and no shape measurement is possible.
    # This check must come before reading image vignets to avoid a C-level malloc
    # crash that occurs when large numbers of numpy arrays are allocated then freed.
    psf_idx = 3 if bkg_sub else 2
    psf_vignet_path = input_file_list[psf_idx]
    with SqliteDict(psf_vignet_path) as db:
        psf_keys = list(db.keys())
        n_empty_psf = sum(1 for k in psf_keys if len(db[k]) == 0)
    if psf_keys and n_empty_psf == len(psf_keys):
        w_log.warning(
            f"All {len(psf_keys)} PSF vignet entries are empty in "
            f"{psf_vignet_path} — no PSF coverage for this tile. Skipping ngmix."
        )
        return None, None

    # Check that image vignets are not all empty before initialising ngmix
    image_vignet_path = input_file_list[1]
    with SqliteDict(image_vignet_path) as db:
        keys = list(db.keys())
        n_empty = sum(1 for k in keys if db[k] == "empty")
    if keys and n_empty == len(keys):
        w_log.warning(
            f"All {len(keys)} image vignets are 'empty' in {image_vignet_path} "
            "— no valid CCD coverage for this tile. Skipping ngmix."
        )
        return None, None

    # Initialise class instance
    ngmix_inst = Ngmix(
        input_file_list[:wcs_idx],
        run_dirs["output"],
        file_number_string,
        zero_point,
        pixel_scale,
        f_wcs_path,
        w_log,
        check_existing_dir=check_existing_dir,
        save_batch=save_batch,
        id_obj_min=id_obj_min,
        id_obj_max=id_obj_max,
        bkg_sub=bkg_sub,
    )

    # Process ngmix shape measurement and metacalibration
    w_log.info("ngmix processing start")
    ngmix_inst.process()
    w_log.info("ngmix end")

    # No return objects
    return None, None
