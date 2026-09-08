"""NGMIX RUNNER.

Module runner for ``ngmix``.

:Author: Axel Guinot

"""

import os

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

    # Pixel scale -- optional override. When absent (or non-positive) it is
    # derived from the image WCS inside Ngmix, so it cannot drift from the
    # pixels. Only the centroid-prior width and noise window use it.
    if config.has_option(module_config_sec, "PIXEL_SCALE"):
        pixel_scale = config.getfloat(module_config_sec, "PIXEL_SCALE")
    else:
        pixel_scale = None

    # Background subtraction: disable for image sims, where there is no
    # background image to subtract. The background vignet occupies one input
    # slot, so the f_wcs headers -- and the optional background-rms vignet --
    # shift by one when it is absent.
    bkg_sub = config.getboolean(module_config_sec, "BKG_SUB", fallback=True)
    wcs_idx = 6 if bkg_sub else 5

    # Path to merged single-exposure single-HDU headers
    f_wcs_path = input_file_list[wcs_idx]

    # Optional background-rms vignet -> ngmix inverse-variance weight map (#779);
    # all-or-nothing. Appended as the final positional input after the wcs slice.
    if config.has_option(module_config_sec, "BKG_RMS_VIGNET_PATH"):
        bkg_rms_vignet_path = config.getexpanded(
            module_config_sec,
            "BKG_RMS_VIGNET_PATH",
        ).format(file_number_string=file_number_string)
        if not os.path.exists(bkg_rms_vignet_path):
            raise FileNotFoundError(
                f"Background RMS vignet file not found: {bkg_rms_vignet_path}"
            )
        input_file_list = input_file_list[:wcs_idx] + [bkg_rms_vignet_path]
    else:
        input_file_list = input_file_list[:wcs_idx]

    # SEG_VIGNET_PATH (optional): coadd-frame SExtractor segmentation vignets
    # (a CLASSIC-mode vignetmaker output), row-aligned to the tile catalogue.
    # Required for BLEND_HANDLING = uberseg; when set, the file must exist for
    # every tile (missing file -> error). Read on Tile_cat, not via Vignet, so
    # it is threaded to Ngmix as its own argument rather than into
    # input_file_list.
    if config.has_option(module_config_sec, "SEG_VIGNET_PATH"):
        seg_vignet_path = config.getexpanded(
            module_config_sec,
            "SEG_VIGNET_PATH",
        ).format(file_number_string=file_number_string)
        if not os.path.exists(seg_vignet_path):
            raise FileNotFoundError(
                f"Segmentation vignet file not found: {seg_vignet_path}"
            )
    else:
        seg_vignet_path = None

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

    # Centroid source for the galaxy Jacobian origin: "wcs" (default -- the
    # catalog sky position projected through the WCS, trusting the astrometry)
    # or "hsm" (legacy HSM adaptive-moment centroid, being phased out: noisy
    # for stars and flagged as incorrect by Fabian -- see #767).
    if config.has_option(module_config_sec, "CENTROID_SOURCE"):
        centroid_source = config.get(module_config_sec, "CENTROID_SOURCE")
    else:
        centroid_source = "wcs"

    # Neighbour treatment: "noisefill" (default, historical) replaces a
    # neighbour's pixels with a noise realisation; "uberseg" hard-masks
    # (weight -> 0) every pixel closer to a neighbour than to the central
    # object, from the segmentation map. See the ngmix module docstrings.
    if config.has_option(module_config_sec, "BLEND_HANDLING"):
        blend_handling = config.get(module_config_sec, "BLEND_HANDLING")
    else:
        blend_handling = "noisefill"

    # DILATE_NEIGHBOUR (optional): binary-dilation iterations enlarging the
    # uberseg neighbour mask, to absorb the few-pixel coadd-vs-epoch seg-overlay
    # offset. Ignored unless BLEND_HANDLING = uberseg. Default 1 (~one pixel).
    if config.has_option(module_config_sec, "DILATE_NEIGHBOUR"):
        dilate_neighbour = config.getint(module_config_sec, "DILATE_NEIGHBOUR")
    else:
        dilate_neighbour = 1

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
        bkg_sub=bkg_sub,
        centroid_source=centroid_source,
        blend_handling=blend_handling,
        seg_cat_path=seg_vignet_path,
        dilate_neighbour=dilate_neighbour,
        seed_from_position=seed_from_position,
        metacal_psf=metacal_psf,
    )

    # Process ngmix shape measurement and metacalibration
    w_log.info("ngmix processing start")
    ngmix_inst.process()
    w_log.info("ngmix end")

    # No return objects
    return None, None
