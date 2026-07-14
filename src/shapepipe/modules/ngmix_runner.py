"""NGMIX RUNNER.

Module runner for ``ngmix``.

:Author: Axel Guinot

"""

import os

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

    # Path to merged single-exposure single-HDU headers
    f_wcs_path = input_file_list[6]

    if config.has_option(module_config_sec, "BKG_RMS_VIGNET_PATH"):
        bkg_rms_vignet_path = config.getexpanded(
            module_config_sec,
            "BKG_RMS_VIGNET_PATH",
        ).format(file_number_string=file_number_string)
        if not os.path.exists(bkg_rms_vignet_path):
            raise FileNotFoundError(
                f"Background RMS vignet file not found: {bkg_rms_vignet_path}"
            )
        input_file_list = input_file_list[:6] + [bkg_rms_vignet_path]
    else:
        input_file_list = input_file_list[:6]

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

    # Centroid source for the galaxy Jacobian origin: "wcs" (default — the
    # catalog sky position projected through the WCS, trusting the astrometry)
    # or "hsm" (legacy HSM adaptive-moment centroid, being phased out: noisy
    # for stars and flagged as incorrect by Fabian — see #767).
    if config.has_option(module_config_sec, "CENTROID_SOURCE"):
        centroid_source = config.get(module_config_sec, "CENTROID_SOURCE")
    else:
        centroid_source = "wcs"

    # Seed the per-object RNG from sky position instead of per tile, so
    # metacal's fixnoise counter-noise (and the fit guesses) cancel across
    # Pujol image-simulation shear branches (ngmix#796). Default False leaves
    # the production path byte-identical.
    if config.has_option(module_config_sec, "SEED_FROM_POSITION"):
        seed_from_position = config.getboolean(
            module_config_sec, "SEED_FROM_POSITION"
        )
    else:
        seed_from_position = False

    # Metacal reconvolution-kernel scheme (metacal_pars['psf']): "fitgauss"
    # (default; fit a Gaussian to the PSF and round it), "gauss" (fixed round
    # Gaussian sized from the PSF), "dilate" (dilate the original PSF), or
    # "azgauss" (ngmix >= 2.4.1; noise-robust variant of "gauss").
    if config.has_option(module_config_sec, "METACAL_PSF"):
        metacal_psf = config.get(module_config_sec, "METACAL_PSF")
    else:
        metacal_psf = "fitgauss"

    # Initialise class instance
    ngmix_inst = Ngmix(
        input_file_list,
        run_dirs["output"],
        file_number_string,
        zero_point,
        pixel_scale,
        f_wcs_path,
        w_log,
        save_batch=save_batch,
        id_obj_min=id_obj_min,
        id_obj_max=id_obj_max,
        centroid_source=centroid_source,
        seed_from_position=seed_from_position,
        metacal_psf=metacal_psf,
    )

    # Process ngmix shape measurement and metacalibration
    w_log.info("ngmix processing start")
    ngmix_inst.process()
    w_log.info("ngmix end")

    # No return objects
    return None, None
