"""STAR COLOUR RUNNER.

Module runner for ``star_colour``: between ``setools_runner`` and
``psfex_runner``, drop training stars an external purity verdict rejects and
give every star a ``GI_COLOR`` for a colour-dependent PSF model.

Config ([STAR_COLOUR_RUNNER]):
    COLOUR_CATALOGUE  FITS table with RA, DEC, the colour and the rejection flag
    COLOUR_COLUMN     colour column name (default GI)
    REJECT_COLUMN     rejection column name (default REMOVE)
    MATCH_RADIUS      arcsec (default 0.5)

:Author: Claude (Opus 5.5) for Cail Daley

"""

import os

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.star_colour_package.star_colour import StarColour


@module_runner(
    input_module="setools_runner",
    version="1.0",
    file_pattern=["star_split_ratio_80", "star_split_ratio_20"],
    file_ext=[".fits", ".fits"],
    depends=["numpy", "astropy", "scipy"],
)
def star_colour_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Star Colour Runner."""
    def opt(key, default):
        return config.getexpanded(module_config_sec, key) if config.has_option(module_config_sec, key) else default

    sc = StarColour(
        config.getexpanded(module_config_sec, "COLOUR_CATALOGUE"),
        colour_col=opt("COLOUR_COLUMN", "GI"),
        reject_col=opt("REJECT_COLUMN", "REMOVE"),
        radius_arcsec=float(opt("MATCH_RADIUS", 0.5)),
    )
    train, reserve = input_file_list
    out = run_dirs["output"]
    sc.process(train, reserve, f"{out}/{os.path.basename(train)}", f"{out}/{os.path.basename(reserve)}", w_log)
    return None, None
