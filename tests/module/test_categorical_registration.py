"""Native-grid categories from independent tile and exposure astrometry."""

from configparser import ConfigParser
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from astropy.wcs import Sip, WCS
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package import ngmix as n
from shapepipe.modules.vignetmaker_package.vignetmaker import (
    VignetMaker,
    get_stamps,
    register_categories,
)

pytestmark = pytest.mark.decision("preparation.categorical_registration")


def _wcs(crpix=(41.0, 43.0), angle=0.0, distortion=False):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [188.0, 30.7]
    w.wcs.crpix = crpix
    theta = np.deg2rad(angle)
    w.wcs.cd = (
        0.18577
        / 3600
        * np.array(
            [
                [-np.cos(theta), np.sin(theta)],
                [np.sin(theta), np.cos(theta)],
            ]
        )
    )
    if distortion:
        a, b = np.zeros((4, 4)), np.zeros((4, 4))
        a[2, 0], a[1, 1], b[0, 2] = 8e-4, -5e-4, -7e-4
        w.sip = Sip(a, b, None, None, w.wcs.crpix)
        w.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
    return w


def _segmentation():
    seg = np.zeros((89, 91), np.int32)
    # An asymmetric target, a departed detection, and a nearby neighbour.
    seg[35:46, 33:38] = 101
    seg[40:44, 38:46] = 101
    seg[45:49, 34:36] = 101
    seg[30:37, 46:50] = -1
    seg[43:52, 48:56] = 205
    return seg


def _oracle(seg, tw, ew, centre, shape):
    """Scalar high-level WCS oracle; no shared mapping/extraction helper."""
    expected = np.zeros(shape, np.int32)
    off = np.zeros(shape, bool)
    for row, col in np.ndindex(shape):
        sky = ew.pixel_to_world(
            col + centre[1] - shape[1] // 2, row + centre[0] - shape[0] // 2
        )
        tx, ty = tw.world_to_pixel(sky)
        x, y = int(np.floor(tx + 0.5)), int(np.floor(ty + 0.5))
        off[row, col] = not (0 <= x < seg.shape[1] and 0 <= y < seg.shape[0])
        if not off[row, col]:
            expected[row, col] = seg[y, x]
    return expected, off


@pytest.mark.parametrize(
    "crpix,angle,distortion",
    [
        ((43.0, 40.0), 0.0, False),  # integer translation
        ((40.37, 42.61), 0.0, False),  # independent subpixel origins
        ((41.24, 42.83), 17.0, False),  # residual rotation
        ((40.81, 43.22), 181.8, True),  # upside-down CCD plus distortion
    ],
)
def test_independent_wcs_registration(crpix, angle, distortion):
    seg, tw, ew = _segmentation(), _wcs(), _wcs(crpix, angle, distortion)
    positions = np.array([[41.29, 40.38]])
    image = np.arange(90 * 90).reshape(90, 90)
    stamps, ip, offsets = get_stamps(image, positions, 12)
    got = register_categories(seg, tw, ew, ip[0], stamps.shape[1:], 101)
    expected, off = _oracle(seg, tw, ew, ip[0], stamps.shape[1:])
    npt.assert_array_equal(got["SEG_VIGNET"], expected)
    npt.assert_array_equal(got["OFF_TILE"], off)
    npt.assert_array_equal(
        got["NEIGHBOUR"], (expected != 0) & (expected != 101)
    )
    npt.assert_allclose(ip[0] + offsets[0], positions[0])
    assert np.any(expected == 101) and np.any(got["NEIGHBOUR"])


def test_rotated_tile_edge_is_not_a_border_run_heuristic():
    tw, ew = _wcs(), _wcs((42.3, 43.2), 27.0)
    seg = _segmentation()
    centre = [2, 3]
    got = register_categories(seg, tw, ew, centre, (19, 19), 101)
    expected, off = _oracle(seg, tw, ew, centre, (19, 19))
    npt.assert_array_equal(got["SEG_VIGNET"], expected)
    npt.assert_array_equal(got["OFF_TILE"], off)
    assert np.any(off) and np.any(~off)
    # At least one off-tile cell has both in-tile row and column companions.
    assert np.any(
        off & (~off).any(axis=1)[:, None] & (~off).any(axis=0)[None, :]
    )
    assert not np.any(got["NEIGHBOUR"] & off)


def test_half_pixel_ties_go_upward():
    # An exact coordinate stub isolates floating-point tie semantics from
    # spherical WCS round-trip precision, exercised in the tests above.
    identity = SimpleNamespace(all_pix2world=lambda xy, origin: xy)
    shift = SimpleNamespace(all_world2pix=lambda xy, origin: xy + 0.5)
    seg = np.arange(25, dtype=np.int32).reshape(5, 5)
    got = register_categories(seg, shift, identity, [2, 2], (3, 3), 12)
    npt.assert_array_equal(got["SEG_VIGNET"], seg[2:5, 2:5])


def test_full_source_covers_pixels_beyond_original_vignet():
    # A target shifted 30 px from the tile barycentre still samples real
    # labels; a 51-px tile cutout padded with sky cannot represent this.
    tw, ew = _wcs(), _wcs((11.0, 43.0))
    seg = _segmentation()
    seg[40:46, 65:72] = 101
    got = register_categories(seg, tw, ew, [42, 40], (51, 51), 101)
    expected, off = _oracle(seg, tw, ew, [42, 40], (51, 51))
    npt.assert_array_equal(got["SEG_VIGNET"], expected)
    npt.assert_array_equal(got["OFF_TILE"], off)
    assert got["SEG_VIGNET"][25, 25] == 101


def _prepare(epoch, tw, ew, obj_id, tile_seg, weight=None, flag=None):
    shape = epoch["VIGNET"].shape
    name = "2254827-27"
    per = lambda a: {str(obj_id): {name: {"VIGNET": a}}}
    vignet = SimpleNamespace(
        bkg_vign_cat=None,
        bkg_rms_vign_cat=None,
        flag_vign_cat=per(np.zeros(shape) if flag is None else flag),
        weight_vign_cat=per(np.ones(shape) if weight is None else weight),
        f_wcs_file={
            "2254827": {
                27: {
                    "WCS": ew,
                    "header": fits.Header({"FSCALE": 1.0}).tostring(),
                }
            }
        },
    )
    sky = ew.pixel_to_world(epoch["INT_POS"][1], epoch["INT_POS"][0])
    tile = SimpleNamespace(
        vign=np.zeros((1, *shape)),
        seg=tile_seg[None],
        ra=[sky.ra.deg],
        dec=[sky.dec.deg],
    )
    return n.prepare_postage_stamps(
        vignet,
        obj_id,
        0,
        tile,
        bkg_sub=False,
        psf_obj={name: {"VIGNET": np.ones(shape)}},
        gal_obj={name: epoch},
        blend_handling="uberseg",
    )


def test_ngmix_consumes_registered_maps_without_flipping(monkeypatch):
    seg, tw, ew = _segmentation(), _wcs(), _wcs((40.37, 42.61), 181.2, True)
    centre, shape = np.array([42, 40]), (51, 51)
    epoch = dict(
        VIGNET=np.random.default_rng(3).normal(size=shape),
        INT_POS=centre,
        OFFSET=np.array([0.2, -0.3]),
    )
    epoch.update(register_categories(seg, tw, ew, centre, shape, 101))
    expected, _ = _oracle(seg, tw, ew, centre, shape)
    monkeypatch.setattr(
        n.Ngmix,
        "MegaCamFlip",
        lambda *a: pytest.fail("native maps must not be flipped"),
    )
    stamp = _prepare(epoch, tw, ew, 101, seg[17:68, 15:66])
    assert len(stamp.gals) == 1
    npt.assert_array_equal(stamp.segs[0], expected)
    npt.assert_array_equal(stamp.neighbours[0], epoch["NEIGHBOUR"])
    npt.assert_array_equal(stamp.offsets[0], epoch["OFFSET"])


def test_off_tile_reaches_fraction_cut_before_diagnostics():
    tw, ew = _wcs(), _wcs((41.2, 43.4), 12.0)
    shape, centre = (51, 51), [1, 2]
    epoch = dict(
        VIGNET=np.random.default_rng(2).normal(size=shape), INT_POS=centre
    )
    epoch.update(
        register_categories(_segmentation(), tw, ew, centre, shape, 101)
    )
    assert epoch["OFF_TILE"].mean() > 0.1
    stamp = _prepare(epoch, tw, ew, 101, np.zeros(shape, np.int32))
    assert stamp.epoch_cuts["masked_fraction"] == 1
    assert not stamp.gals
    assert stamp.n_epoch_interp == 0


def test_retained_188307_object_5422_epoch_2254827_27(monkeypatch):
    path = Path(__file__).parents[1] / "data/registration_188307.npz"
    with np.load(path) as f:
        tw = WCS(fits.Header.fromstring(str(f["tile_header"])))
        ew = WCS(fits.Header.fromstring(str(f["epoch_header"])))
        epoch = dict(
            VIGNET=f["image"], INT_POS=f["int_pos"], OFFSET=f["offset"]
        )
        epoch.update(
            register_categories(
                f["segmentation"], tw, ew, f["int_pos"], (51, 51), 5422
            )
        )
        npt.assert_array_equal(epoch["SEG_VIGNET"], f["expected_seg"])
        assert (
            (epoch["SEG_VIGNET"] == 5422) ^ (f["old_seg"] == 5422)
        ).sum() == 10
        assert f["old_seg"][25, 27] == 0
        assert epoch["SEG_VIGNET"][25, 27] == 5422
        assert f["weight"][25, 27] == 0 and f["flag"][25, 27] == 3
        # Isolate the footprint diagnostic; production's fixed-radius veto
        # rejects this particular epoch in either registration.
        monkeypatch.setattr(n, "central_defect_vetoes", lambda *args: False)
        stamp = _prepare(
            epoch, tw, ew, 5422, f["old_seg"], f["weight"], f["flag"]
        )
        assert len(stamp.gals) == 1
        assert stamp.imaflags_iso == 3
        assert stamp.n_epoch_interp_iso == 1


@pytest.mark.parametrize(
    "relative",
    [
        "workflow/config/cfis/config_tile_PiViVi_psfex.ini",
        "workflow/config/cfis/config_tile_PiViVi_mccd.ini",
        "workflow/config/cfis_image_sims/config_tile_PiViVi_fake.ini",
        "example/cfis_image_sims/config_tile_PiViVi_canfar_sx.ini",
    ],
)
def test_multi_epoch_configs_supply_full_segmentation(relative):
    config = ConfigParser(interpolation=None)
    config.read(Path(__file__).parents[2] / relative)
    sections = [
        s
        for s in config.sections()
        if s.startswith("VIGNETMAKER")
        and config.get(s, "MODE") == "MULTI-EPOCH"
    ]
    assert sections
    for section in sections:
        assert (
            config.get(section, "FILE_PATTERN").split(",")[-1].strip()
            == "segmentation"
        )
        assert (
            config.get(section, "FILE_EXT").split(",")[-1].strip() == ".fits"
        )


def test_simulation_example_produces_registration_labels():
    root = Path(__file__).parents[2] / "example/cfis_image_sims"
    config = ConfigParser(interpolation=None)
    config.read(root / "config_tile_Sx.ini")
    assert config.getboolean("SEXTRACTOR_RUNNER", "SEG_VIGNET")
    assert "SEGMENTATION" in config.get("SEXTRACTOR_RUNNER", "CHECKIMAGE")


def test_vignetmaker_writes_joined_categories_with_image_grid(tmp_path):
    seg, tw, ew = _segmentation(), _wcs(), _wcs((40.37, 42.61), 1.7, True)
    raw = np.where(seg == 101, 7, np.where(seg != 0, 8, 0)).astype(np.int32)
    sky = ew.pixel_to_world(40.27, 42.34)
    columns = [
        fits.Column(name="NUMBER", format="J", array=[101]),
        fits.Column(name="SEG_LABEL", format="J", array=[7]),
        fits.Column(name="N_EPOCH", format="J", array=[1]),
        fits.Column(name="XWIN_WORLD", format="D", array=[sky.ra.deg]),
        fits.Column(name="YWIN_WORLD", format="D", array=[sky.dec.deg]),
    ]
    cat = tmp_path / "cat.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(),
            fits.BinTableHDU.from_columns(columns),
            fits.BinTableHDU.from_columns(
                [
                    fits.Column(
                        name="EXP_NAME", format="8A", array=["2254827"]
                    ),
                    fits.Column(name="CCD_N", format="J", array=[27]),
                ],
                name="EPOCH_0",
            ),
        ]
    ).writeto(cat)
    seg_path, image_path = (
        tmp_path / "segmentation.fits",
        tmp_path / "image-2254827-27.fits",
    )
    fits.writeto(seg_path, raw, tw.to_header(relax=True))
    image = np.random.default_rng(4).normal(size=(89, 91))
    fits.writeto(image_path, image)
    headers = tmp_path / "headers.sqlite"
    with SqliteDict(str(headers)) as db:
        db["2254827"] = {27: {"WCS": ew}}
        db.commit()
    maker = VignetMaker(
        str(cat),
        "SPHE",
        ["XWIN_WORLD", "YWIN_WORLD"],
        str(tmp_path),
        "-188-307",
        SimpleNamespace(info=lambda msg: None),
    )
    maker.process_me(
        [[str(tmp_path)]], ["image"], str(headers), 25, str(seg_path)
    )
    with SqliteDict(
        str(tmp_path / "image_vignet-188-307.sqlite"), flag="r"
    ) as db:
        epoch = db["101"]["2254827-27"]
    stamps, ip, offset = get_stamps(image, np.array([[42.34, 40.27]]), 25)
    npt.assert_array_equal(epoch["VIGNET"], stamps[0])
    npt.assert_array_equal(epoch["INT_POS"], ip[0])
    npt.assert_allclose(epoch["OFFSET"], offset[0], atol=1e-7)
    joined = np.where(raw == 7, 101, np.where(raw != 0, -1, 0))
    expected, off = _oracle(joined, tw, ew, ip[0], (51, 51))
    npt.assert_array_equal(epoch["SEG_VIGNET"], expected)
    npt.assert_array_equal(epoch["OFF_TILE"], off)
    assert np.any(epoch["SEG_VIGNET"] == -1)
    npt.assert_array_equal(epoch["NEIGHBOUR"], expected == -1)
