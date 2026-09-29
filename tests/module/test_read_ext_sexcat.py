"""UNIT TESTS FOR MODULE PACKAGE: READ_EXT_SEXCAT.

Drives ``make_ldac_from_ascii`` on a synthetic ASCII SExtractor-format
catalogue and a synthetic tile image, and checks the FITS-LDAC it writes is
what the tile chain downstream of ``tile_detect`` reads: the LDAC_IMHEAD
extension carrying the tile header, the SExtractor column aliases, one
``VIGNET`` stamp per object cut from the image, and the input ``NUMBER``
kept as is. It follows the catalogue through
``make_cat.save_sextractor_data``, which builds ``TILE_UNIQUE_ID``. The rest
covers the segmentation map: relabelling it to the catalogue's ``NUMBER`` and
setting neighbours' ``VIGNET`` pixels to -1e30, as SExtractor does, which is
all ngmix reads to mask neighbours.
"""

from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

from shapepipe.modules.make_cat_package import make_cat
from shapepipe.modules.read_ext_sexcat_package import read_ext_sexcat as rs

NX, NY = 40, 30
STAMP = 5
# (NUMBER, X_IMAGE, Y_IMAGE): an interior object, one on the left edge, one
# in the top-right corner.
OBJECTS = [(1, 10.0, 12.0), (2, 1.0, 20.0), (7, 40.0, 30.0)]


def _write_ascii_cat(path):
    lines = [
        "#   1 NUMBER          Running object number",
        "#   2 X_IMAGE         Object position along x        [pixel]",
        "#   3 Y_IMAGE         Object position along y        [pixel]",
        "#   4 ALPHA_J2000     Right ascension of barycenter  [deg]",
        "#   5 DELTA_J2000     Declination of barycenter      [deg]",
        "#   6 MAG_AUTO        Kron-like elliptical aperture magnitude [mag]",
    ]
    for num, x, y in OBJECTS:
        lines.append(f"{num} {x} {y} {150.0 + num} {30.0 + num} {20.0 + num}")
    path.write_text("\n".join(lines) + "\n")


def _write_image(path):
    # pixel value = 1000*row + column (0-based), so every stamp pixel names
    # where it came from.
    data = (np.arange(NY)[:, None] * 1000 + np.arange(NX)[None, :]).astype(
        np.float32
    )
    hdu = fits.PrimaryHDU(data)
    hdu.header["HISTORY"] = "input image 2605805p.fits"
    hdu.header["TILEKEY"] = "kept"
    hdu.writeto(path, overwrite=True)


@pytest.fixture
def ldac(tmp_path):
    cat = tmp_path / "CFIS_cat-301-279.cat"
    img = tmp_path / "CFIS_image-301-279.fits"
    out = tmp_path / "sexcat-301-279.fits"
    _write_ascii_cat(cat)
    _write_image(img)
    rs.make_ldac_from_ascii(
        str(cat), str(img), str(out), stamp_size=STAMP
    )
    return out


def test_ldac_layout_and_header(ldac):
    with fits.open(ldac) as hdul:
        assert [h.name for h in hdul] == ["PRIMARY", "LDAC_IMHEAD", "LDAC_OBJECTS"]
        cards = hdul["LDAC_IMHEAD"].data[0][0]
        assert isinstance(cards, str) or cards.ndim == 1
        text = "".join(cards) if not isinstance(cards, str) else cards
        assert "TILEKEY" in text and "2605805p" in text


def test_number_is_kept_and_aliases_added(ldac):
    with fits.open(ldac) as hdul:
        data = hdul["LDAC_OBJECTS"].data
    # The input NUMBER (gapped here: 1, 2, 7) is the object's identity and is
    # copied unchanged; the converter builds no ID of its own.
    npt.assert_array_equal(data["NUMBER"], [o[0] for o in OBJECTS])
    assert "TILE_UNIQUE_ID" not in data.names
    npt.assert_array_equal(data["XWIN_IMAGE"], data["X_IMAGE"])
    npt.assert_array_equal(data["YWIN_IMAGE"], data["Y_IMAGE"])
    npt.assert_array_equal(data["XWIN_WORLD"], data["ALPHA_J2000"])
    npt.assert_array_equal(data["YWIN_WORLD"], data["DELTA_J2000"])


def test_vignets_are_cut_from_the_image_and_padded_as_sextractor(ldac):
    with fits.open(ldac) as hdul:
        vignets = hdul["LDAC_OBJECTS"].data["VIGNET"]
    assert vignets.shape == (len(OBJECTS), STAMP, STAMP)

    # Interior object at (10, 12), 1-based: centre pixel is (row 11, col 9).
    assert vignets[0, STAMP // 2, STAMP // 2] == 11 * 1000 + 9
    assert vignets[0, 0, 0] == 9 * 1000 + 7

    # Left edge, x = 1: the two columns left of the image are -1e30, the value
    # SExtractor writes off the image.
    assert (vignets[1, :, :2] == rs.BIG).all()
    assert vignets[1, STAMP // 2, STAMP // 2] == 19 * 1000 + 0

    # Top-right corner: only the lower-left quadrant of the stamp is in the
    # image.
    assert (vignets[2, STAMP // 2 + 1:, :] == rs.BIG).all()
    assert (vignets[2, :, STAMP // 2 + 1:] == rs.BIG).all()
    assert vignets[2, STAMP // 2, STAMP // 2] == 29 * 1000 + 39


def test_tile_unique_id_reaches_the_final_catalogue(ldac, tmp_path):
    """make_cat builds the ID from the tile and the catalogue's own NUMBER."""
    final = make_cat.prepare_final_cat_file(str(tmp_path), "-301-279")
    n_obj = make_cat.save_sextractor_data(final, str(ldac))
    assert n_obj == len(OBJECTS)
    with fits.open(tmp_path / "final_cat-301-279.fits") as hdul:
        data = hdul["RESULTS"].data
    assert "VIGNET" not in data.names
    npt.assert_array_equal(
        data["TILE_UNIQUE_ID"], 301279 * 10**6 + np.array([1, 2, 7])
    )
    npt.assert_allclose(data["TILE_ID"], 301.279)


# --- the segmentation map -------------------------------------------------


def _seg_map():
    """Two footprints, labelled 7 and 9, on a 20x20 sky."""
    seg = np.zeros((20, 20), dtype=np.int32)
    seg[2:6, 2:6] = 7
    seg[12:18, 12:18] = 9
    return seg


@pytest.mark.decision("detection.catalogue_neighbour_marking")
class TestRelabel:
    """The relabelled map carries each object's NUMBER on its own footprint.

    Segmentation labels are not the catalogue's NUMBER; the claim is by the
    pixel under the object's position, which is all both maps share.
    """

    def test_each_object_owns_the_footprint_it_sits_in(self):
        seg = _seg_map()
        # Positions in an order that is NOT the label order, so a
        # relabelling that merely renumbered would fail.
        out, counts = rs.relabel_seg(
            seg, number=np.array([1, 2]),
            x_image=np.array([15.0, 4.0]), y_image=np.array([15.0, 4.0]))
        assert counts["matched"] == 2
        assert set(np.unique(out[seg == 9])) == {1}
        assert set(np.unique(out[seg == 7])) == {2}
        assert np.all(out[seg == 0] == 0)

    def test_an_unclaimed_footprint_becomes_a_neighbour(self):
        seg = _seg_map()
        out, _ = rs.relabel_seg(
            seg, number=np.array([1]),
            x_image=np.array([4.0]), y_image=np.array([4.0]))
        assert set(np.unique(out[seg == 9])) == {rs.NEIGHBOUR_LABEL}
        assert set(np.unique(out[seg == 7])) == {1}

    def test_an_object_on_sky_gets_a_disc(self):
        seg = _seg_map()
        out, counts = rs.relabel_seg(
            seg, number=np.array([1, 5]),
            x_image=np.array([4.0, 10.0]), y_image=np.array([4.0, 10.0]),
            fallback_radius=2)
        assert counts["unclaimed"] == 1
        assert out[9, 9] == 5
        assert np.count_nonzero(out == 5) == np.count_nonzero(
            np.add.outer(np.arange(-2, 3) ** 2, np.arange(-2, 3) ** 2) <= 4)

    def test_two_objects_in_one_footprint_both_keep_a_centre(self):
        seg = _seg_map()
        out, counts = rs.relabel_seg(
            seg, number=np.array([1, 2]),
            x_image=np.array([14.0, 16.0]), y_image=np.array([14.0, 16.0]),
            fallback_radius=1)
        assert counts == dict(matched=1, unclaimed=0, shared=1, off_image=0,
                              shared_pixel=0)
        assert out[13, 13] == 1 and out[15, 15] == 2
        assert np.count_nonzero(out == 1) > np.count_nonzero(out == 2)

    def test_every_object_is_self_somewhere(self):
        rng = np.random.default_rng(0)
        seg = np.zeros((60, 60), dtype=np.int32)
        for label in range(1, 12):
            row, col = rng.integers(0, 55, size=2)
            seg[row:row + 5, col:col + 5] = label
        number = np.arange(1, 31)
        x_image = rng.uniform(1, 60, size=30)
        y_image = rng.uniform(1, 60, size=30)
        out, counts = rs.relabel_seg(seg, number, x_image, y_image)
        assert sum(counts[k] for k in
                   ("matched", "unclaimed", "shared", "off_image")) == 30
        for num in number:
            assert np.any(out == num), f"object {num} has no self pixels"
        assert set(np.unique(out)) <= set(number) | {0, rs.NEIGHBOUR_LABEL}

    def test_two_objects_on_one_pixel(self):
        seg = _seg_map()
        out, counts = rs.relabel_seg(
            seg, number=np.array([4, 6]), x_image=np.array([10.0, 10.2]),
            y_image=np.array([10.0, 10.1]), fallback_radius=1)
        assert counts["shared_pixel"] == 1
        assert out[9, 9] == 4
        assert np.any(out == 6)

    @pytest.mark.parametrize("x, y", [(0.4, 5.0), (5.0, 21.0)])
    def test_a_position_off_the_image_is_counted(self, x, y):
        out, counts = rs.relabel_seg(
            _seg_map(), number=np.array([1]), x_image=np.array([x]),
            y_image=np.array([y]))
        assert counts["off_image"] == 1
        assert not np.any(out == 1)


# Objects on _seg_map(): 1 in footprint 7, 2 on sky between the footprints;
# footprint 9 is claimed by nobody.
SEG_OBJECTS = [(1, 4.0, 4.0), (2, 10.0, 9.0)]
SEG_STAMP = 21


def _marked(number, x, y, seg=None):
    image = np.ones(_seg_map().shape, np.float32)
    seg = _seg_map() if seg is None else seg
    relabelled, _ = rs.relabel_seg(seg, number, x, y)
    return rs._extract_vignets(image, x, y, SEG_STAMP, seg=relabelled,
                               number=number)


def _stamp_of(array, x, y, fill):
    """The SEG_STAMP stamp of ``array`` centred on 1-based (x, y)."""
    half = SEG_STAMP // 2
    padded = np.pad(array, half, constant_values=fill)
    row, col = int(round(y)) - 1, int(round(x)) - 1
    return padded[row:row + SEG_STAMP, col:col + SEG_STAMP]


@pytest.mark.decision("detection.catalogue_neighbour_marking")
def test_neighbour_footprints_become_big_and_nothing_else():
    number = np.array([o[0] for o in SEG_OBJECTS])
    x = np.array([o[1] for o in SEG_OBJECTS])
    y = np.array([o[2] for o in SEG_OBJECTS])
    vignets = _marked(number, x, y)
    seg = _seg_map()

    # Object 1 in footprint 7: footprint 9 (unclaimed) is a neighbour; its
    # own footprint and the sky keep their image values.
    s = _stamp_of(seg, x[0], y[0], fill=-99)
    v = vignets[0]
    assert (v[s == 9] == rs.BIG).all()
    assert (v[s == 7] == 1).all()
    assert (v[s == -99] == rs.BIG).all()  # off the image
    # Object 2's disc is a catalogue object's pixels, so a neighbour too; the
    # rest of the sky is untouched.
    relabelled, _ = rs.relabel_seg(seg, number, x, y)
    disc = _stamp_of(relabelled, x[0], y[0], fill=-99) == 2
    assert disc.any() and (v[disc] == rs.BIG).all()
    assert (v[(s == 0) & ~disc] == 1).all()

    # Object 2 on sky: every footprint in its stamp is a neighbour; the sky,
    # its own centre included, keeps its image values.
    s = _stamp_of(seg, x[1], y[1], fill=-99)
    v = vignets[1]
    assert (v[(s == 7) | (s == 9)] == rs.BIG).all()
    assert (v[s == 0] == 1).all()
    assert v[SEG_STAMP // 2, SEG_STAMP // 2] == 1


@pytest.mark.decision("detection.catalogue_neighbour_marking")
def test_a_claimed_neighbour_is_masked_by_its_number():
    """With both footprints claimed, each object masks the other's."""
    number, x, y = np.array([3, 8]), np.array([4.0, 15.0]), np.array([4.0, 15.0])
    vignets = _marked(number, x, y)
    seg = _seg_map()
    for i, (own, other) in enumerate([(7, 9), (9, 7)]):
        s = _stamp_of(seg, x[i], y[i], fill=-99)
        assert (vignets[i][s == other] == rs.BIG).all()
        assert (vignets[i][s == own] == 1).all()


@pytest.mark.decision("detection.catalogue_neighbour_marking")
def test_converter_writes_the_relabelled_map_and_marks(tmp_path):
    """End to end from a compressed map, as fetched from vos."""
    cat = tmp_path / "CFIS_cat-301-279.cat"
    img = tmp_path / "CFIS_image-301-279.fits"
    seg_in = tmp_path / "CFIS_seg-301-279.fitsfz"
    seg_out = tmp_path / "seg-301-279.fits"
    out = tmp_path / "sexcat-301-279.fits"
    lines = ["#   1 NUMBER", "#   2 X_IMAGE", "#   3 Y_IMAGE",
             "#   4 ALPHA_J2000", "#   5 DELTA_J2000"]
    lines += [f"{n} {x} {y} 150.0 30.0" for n, x, y in SEG_OBJECTS]
    cat.write_text("\n".join(lines) + "\n")
    fits.PrimaryHDU(np.ones((20, 20), np.float32)).writeto(img)
    header = fits.Header()
    header["CRPIX1"] = 10.0
    fits.HDUList([fits.PrimaryHDU(),
                  fits.CompImageHDU(_seg_map(), header=header)]).writeto(seg_in)

    rs.make_ldac_from_ascii(str(cat), str(img), str(out), stamp_size=SEG_STAMP,
                            seg_path=str(seg_in), seg_output_path=str(seg_out))

    with fits.open(seg_out) as hdul:
        relabelled = hdul[0].data
        assert hdul[0].header["CRPIX1"] == 10.0
    seg = _seg_map()
    assert set(np.unique(relabelled[seg == 7])) == {1}
    assert set(np.unique(relabelled[seg == 9])) == {rs.NEIGHBOUR_LABEL}
    with fits.open(out) as hdul:
        v = hdul["LDAC_OBJECTS"].data["VIGNET"][0]
    s = _stamp_of(seg, 4.0, 4.0, fill=-99)
    assert (v[s == 9] == rs.BIG).all() and (v[s == 7] == 1).all()

    fits.HDUList([fits.PrimaryHDU(),
                  fits.CompImageHDU(_seg_map()[:10])]).writeto(
        seg_in, overwrite=True)
    with pytest.raises(ValueError, match="one grid"):
        rs.make_ldac_from_ascii(str(cat), str(img), str(out),
                                stamp_size=SEG_STAMP, seg_path=str(seg_in),
                                seg_output_path=str(seg_out))


DR6_PATCH = Path(__file__).parent / "data" / "dr6_202.301_seg_patch.fits"


@pytest.mark.decision("detection.catalogue_neighbour_marking")
def test_dr6_marks_every_neighbour_pixel_and_no_own_pixel():
    """On a crowded 200x200 patch of the real 202.301 map and catalogue.

    The raw labels are not NUMBER, so this checks the whole chain on real
    data: 1-based positions, the centre-pixel claim, and the stamp geometry.
    Every stamp fully on the patch is checked against the RAW map: pixels of
    footprints other than the one under the object are all -1e30, its own
    footprint and the sky are untouched. (SExtractor's own VIGNET, on a
    SExtractor seg map of an image-sim tile, marks 93% of neighbour-label
    pixels and 0.14% of own pixels: check_sex_vignet.py.)
    """
    with fits.open(DR6_PATCH) as hdul:
        seg = hdul["SEG"].data
        objects = hdul["OBJECTS"].data
    number = np.array(objects["NUMBER"])
    x, y = np.array(objects["X_IMAGE"]), np.array(objects["Y_IMAGE"])
    relabelled, counts = rs.relabel_seg(seg, number, x, y)
    assert counts["matched"] == len(number)

    stamp, half = 51, 25
    image = np.ones(seg.shape, np.float32)
    vignets = rs._extract_vignets(image, x, y, stamp, seg=relabelled,
                                  number=number)
    col, row = np.rint(x).astype(int) - 1, np.rint(y).astype(int) - 1
    full = ((col >= half) & (col < seg.shape[1] - half)
            & (row >= half) & (row < seg.shape[0] - half))
    assert full.sum() >= 10

    marked = {"neighbour": [0, 0], "own": [0, 0], "sky": [0, 0]}
    for i in np.flatnonzero(full):
        s = seg[row[i] - half:row[i] + half + 1, col[i] - half:col[i] + half + 1]
        big = vignets[i] == rs.BIG
        own = s[half, half]
        for kind, where in (("neighbour", (s != 0) & (s != own)),
                            ("own", s == own), ("sky", s == 0)):
            marked[kind][0] += (big & where).sum()
            marked[kind][1] += where.sum()
    assert marked["neighbour"][1] > 1000
    assert marked["neighbour"][0] == marked["neighbour"][1]
    assert marked["own"][0] == 0
    assert marked["sky"][0] == 0
