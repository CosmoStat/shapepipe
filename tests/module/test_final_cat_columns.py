"""The final-catalogue columns each tile_detection mode requests.

``final_cat_merge`` asks every tile catalogue for the columns of
``workflow/config/cfis/final_cat.param`` and stops on any one missing. That
list is written against SExtractor-mode tiles; under ``tile_detection:
unions_catalogue`` the detection columns are whatever the UNIONS per-tile
catalogue carries, copied through ``read_ext_sexcat``. So the catalogue-mode
request is the param list less ``merge_final_cat.SEXTRACTOR_ONLY_COLUMNS``.

Here the catalogue-mode detection columns are not written down but derived:
the converter runs on a catalogue with the UNIONS DR6 header (copied verbatim
from ``vos:cfis/tiles_DR6/CFIS.202.301.r.cat``) and its output is what a tile
can carry. The requested SExtractor columns are those of the tile SExtractor
parameter file. Every other requested column comes from stages that run the
same way in both modes (post-processing, ngmix, make_cat), so the detection
columns are the whole of the difference between them.
"""

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from shapepipe.modules.read_ext_sexcat_package import read_ext_sexcat as rs

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
CFIS = REPO_ROOT / "workflow" / "config" / "cfis"

# The header of a UNIONS DR6 per-tile catalogue and its first object, moved to
# pixel (4, 4) so that its stamp falls on the small test image.
DR6_CATALOGUE = """\
#   1 NUMBER                 Running object number
#   2 X_IMAGE                Object position along x                                    [pixel]
#   3 Y_IMAGE                Object position along y                                    [pixel]
#   4 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#   5 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#   6 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
#   7 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
#   8 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
#   9 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
#  10 MAG_APER               Fixed aperture magnitude vector                            [mag]
#  11 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#  12 A_WORLD                Profile RMS along major axis (world units)                 [deg]
#  13 ERRA_WORLD             World RMS position error along major axis                  [deg]
#  14 B_WORLD                Profile RMS along minor axis (world units)                 [deg]
#  15 ERRB_WORLD             World RMS position error along minor axis                  [deg]
#  16 THETA_J2000            Position angle (east of north) (J2000)                     [deg]
#  17 ERRTHETA_J2000         J2000 error ellipse pos. angle (east of north)             [deg]
#  18 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
#  19 MU_MAX                 Peak surface brightness above background                   [mag * arcsec**(-2)]
#  20 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
#  21 FLAGS                  Extraction flags
         1      4.0000      4.0000 205.3679556 +60.2576198  14.8792   0.0002  14.8792   0.0002  14.9100   0.0002  0.000389709  2.80088e-07 0.0003095117 1.974931e-07  -1.48  -0.72      5707  16.3802      4.474   0
"""


def _load_script(name, path):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


merge_final_cat = _load_script("merge_final_cat",
                               SCRIPTS / "merge_final_cat.py")
create_final_cat = _load_script(
    "create_final_cat",
    REPO_ROOT / "scripts" / "python" / "create_final_cat.py")


@pytest.fixture(scope="module")
def param_list():
    return create_final_cat.read_param_file(str(CFIS / "final_cat.param"))


@pytest.fixture(scope="module")
def sextractor_columns():
    """Column names the tile SExtractor run writes (vector sizes stripped)."""
    names = set()
    for line in (CFIS / "default_noimaflags.param").read_text().splitlines():
        entry = line.split("#")[0].strip()
        if entry:
            names.add(re.sub(r"\(.*\)$", "", entry))
    return names


@pytest.fixture(scope="module")
def catalogue_columns(tmp_path_factory):
    """Detection columns of a tile under tile_detection: unions_catalogue."""
    tmp = tmp_path_factory.mktemp("dr6")
    cat, img, out = (tmp / "CFIS.202.301.r.cat", tmp / "CFIS.202.301.r.fits",
                     tmp / "sexcat-202-301.fits")
    cat.write_text(DR6_CATALOGUE)
    fits.PrimaryHDU(np.zeros((8, 8), dtype=np.float32)).writeto(img)
    rs.make_ldac_from_ascii(str(cat), str(img), str(out), stamp_size=3)
    with fits.open(out) as hdul:
        return set(hdul["LDAC_OBJECTS"].columns.names) - {"VIGNET"}


def test_sextractor_only_is_exactly_what_the_catalogue_lacks(
        param_list, sextractor_columns, catalogue_columns):
    """Every requested SExtractor column the catalogue lacks, and no other."""
    lacking = {c for c in param_list
               if c in sextractor_columns and c not in catalogue_columns}
    assert set(merge_final_cat.SEXTRACTOR_ONLY_COLUMNS) == lacking


def test_catalogue_mode_request_is_carried(
        param_list, sextractor_columns, catalogue_columns):
    """Each requested detection column is one the converter writes."""
    requested = merge_final_cat.requested_columns(param_list,
                                                  "unions_catalogue")
    detection = [c for c in requested if c in sextractor_columns]
    assert detection and set(detection) <= catalogue_columns
    # Only detection columns are dropped, and the order is the param file's.
    assert requested == [c for c in param_list
                         if c not in merge_final_cat.SEXTRACTOR_ONLY_COLUMNS]
    assert {"MAG_AUTO", "MAGERR_AUTO", "FLUX_RADIUS"} <= set(requested)


def test_sextractor_mode_requests_the_param_file():
    for input_type in ("cfis", "cfis_image_sims"):
        params = create_final_cat.read_param_file(
            str(REPO_ROOT / "workflow" / "config" / input_type
                / "final_cat.param"))
        assert merge_final_cat.requested_columns(
            params, "sextractor") == params


def test_unknown_mode_is_refused(param_list):
    with pytest.raises(ValueError):
        merge_final_cat.requested_columns(param_list, "sextractr")
