"""Joining the tile SExtractor catalogue to the UNIONS catalogue."""

import numpy as np
import pytest
from astropy.io import fits

from shapepipe.modules.sextractor_package import match_catalogue as mc
from shapepipe.pipeline.config import CustomParser

HEADER = """\
#   1 NUMBER                 Running object number
#   2 X_IMAGE                Object position along x                                    [pixel]
#   3 Y_IMAGE                Object position along y                                    [pixel]
#   4 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
"""


def _external(path, number, x, y):
    rows = "".join(f"{n:10d} {xi:11.4f} {yi:11.4f}  24.0000\n"
                   for n, xi, yi in zip(number, x, y))
    path.write_text(HEADER + rows)
    return str(path)


def _sexcat(path, x, y, seg_vignets=None):
    """A FITS-LDAC sexcat: NUMBER 1..n at (x, y), 3x3 stamps."""
    n = len(x)
    number = np.arange(1, n + 1, dtype=np.int32)
    cols = [
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="X_IMAGE", format="E", array=np.asarray(x)),
        fits.Column(name="Y_IMAGE", format="E", array=np.asarray(y)),
        fits.Column(name="XWIN_IMAGE", format="D",
                    array=np.asarray(x) + 0.01),
        fits.Column(name="VIGNET", format="9E", dim="(3,3)",
                    array=np.arange(n * 9, dtype=np.float32).reshape(n, 9)),
    ]
    if seg_vignets is not None:
        cols.append(fits.Column(name="SEG_VIGNET", format="9J", dim="(3,3)",
                                array=seg_vignets.reshape(n, 9)))
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="80A",
                     array=np.array(["SIMPLE  =                    T"]))],
        name="LDAC_IMHEAD")
    objects = fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(path)
    return str(path)


def _objects(path):
    with fits.open(path) as hdul:
        assert [h.name for h in hdul] == ["PRIMARY", "LDAC_IMHEAD",
                                          "LDAC_OBJECTS"]
        return hdul["LDAC_OBJECTS"].data.copy()


X = np.array([10.0, 50.0, 90.0, 130.0, 170.0])
Y = np.array([20.0, 60.0, 100.0, 140.0, 180.0])


def test_rows_take_the_external_number_one_to_one(tmp_path):
    """Every row pairs with its object, whatever the external order."""
    cat = _sexcat(tmp_path / "sexcat.fits", X, Y)
    order = [3, 0, 4, 1, 2]
    numbers = [1001, 1002, 1003, 1004, 1005]
    # Shuffled, offset by 1e-4 px, plus one object SExtractor did not detect.
    ext = _external(tmp_path / "ext.cat", numbers + [1006],
                    list(X[order] + 1e-4) + [300.0],
                    list(Y[order] - 1e-4) + [300.0])
    counts = mc.match_catalogue(cat, ext, min_fraction=0.8)
    data = _objects(cat)
    expected = np.empty(5, int)
    expected[order] = numbers
    assert list(data["NUMBER"]) == list(expected)
    assert counts == dict(n_sextractor=5, n_external=6, n_paired=5,
                          n_dropped=0, n_external_only=1)
    # Every measurement is SExtractor's, row for row.
    np.testing.assert_array_equal(data["XWIN_IMAGE"], X + 0.01)
    np.testing.assert_array_equal(data["VIGNET"][2].ravel(),
                                  np.arange(18, 27))


def test_unpaired_rows_leave_the_catalogue(tmp_path):
    """A detection with no object within the radius is dropped."""
    cat = _sexcat(tmp_path / "sexcat.fits", X, Y)
    ext = _external(tmp_path / "ext.cat", [7, 8, 9, 10],
                    X[[0, 1, 2, 4]], Y[[0, 1, 2, 4]] + 0.9)
    counts = mc.match_catalogue(cat, ext, radius=1.0, min_fraction=0.5)
    data = _objects(cat)
    assert list(data["NUMBER"]) == [7, 8, 9, 10]
    np.testing.assert_allclose(data["X_IMAGE"], X[[0, 1, 2, 4]])
    assert counts["n_dropped"] == 1 and counts["n_external_only"] == 0


def test_pairs_are_mutual_nearest_neighbours(tmp_path):
    """Two detections near one object: only the nearer one takes it."""
    x = np.array([10.0, 10.6, 90.0])
    y = np.array([10.0, 10.0, 90.0])
    cat = _sexcat(tmp_path / "sexcat.fits", x, y)
    ext = _external(tmp_path / "ext.cat", [5, 6], [10.1, 90.0], [10.0, 90.0])
    mc.match_catalogue(cat, ext, min_fraction=0.5)
    data = _objects(cat)
    assert list(data["NUMBER"]) == [5, 6]
    np.testing.assert_allclose(data["X_IMAGE"], [10.0, 90.0])


def test_pairs_lie_within_the_radius(tmp_path):
    """Mutual nearest neighbours 2 px apart, beyond the radius, do not pair."""
    x = np.append(X, 300.0)
    y = np.append(Y, 300.0)
    cat = _sexcat(tmp_path / "sexcat.fits", x, y)
    ext = _external(tmp_path / "ext.cat", [1, 2, 3, 4, 5, 6],
                    np.append(X, 302.0), np.append(Y, 300.0))
    counts = mc.match_catalogue(cat, ext, radius=1.0, min_fraction=0.8)
    assert list(_objects(cat)["NUMBER"]) == [1, 2, 3, 4, 5]
    assert counts["n_dropped"] == 1 and counts["n_external_only"] == 1


def test_too_few_catalogue_objects_paired_stop_the_run(tmp_path):
    """Fewer detections than catalogue objects (a drifted detection
    configuration) fails, though every detection pairs."""
    cat = _sexcat(tmp_path / "sexcat.fits", X[:4], Y[:4])
    ext = _external(tmp_path / "ext.cat", [1, 2, 3, 4, 5], X, Y)
    with pytest.raises(ValueError, match="fewer objects"):
        mc.match_catalogue(cat, ext, min_fraction=0.99)
    assert list(_objects(cat)["NUMBER"]) == [1, 2, 3, 4]
    counts = mc.match_catalogue(cat, ext, min_fraction=0.8)
    assert counts["n_paired"] == 4 and counts["n_external_only"] == 1


def test_too_few_pairs_stop_the_run(tmp_path):
    """Below MATCH_MIN_FRACTION (other pixels: the DR5 image) the run fails."""
    cat = _sexcat(tmp_path / "sexcat.fits", X, Y)
    ext = _external(tmp_path / "ext.cat", [1, 2, 3, 4], X[:4], Y[:4])
    with pytest.raises(ValueError, match="DR6"):
        mc.match_catalogue(cat, ext, min_fraction=0.99)
    # Nothing was written.
    assert list(_objects(cat)["NUMBER"]) == [1, 2, 3, 4, 5]
    mc.match_catalogue(cat, ext, min_fraction=0.8)
    assert list(_objects(cat)["NUMBER"]) == [1, 2, 3, 4]


def test_seg_vignet_is_relabelled(tmp_path):
    """SEG_VIGNET labels follow NUMBER; a dropped row's footprint is -1."""
    x, y = X[:3], Y[:3]
    seg = np.array([
        [[0, 1, 1], [0, 1, 2], [0, 0, 3]],
        [[2, 2, 0], [1, 2, 0], [0, 0, 0]],
        [[3, 3, 3], [0, 3, 3], [2, 0, 0]],
    ], dtype=np.int32)
    cat = _sexcat(tmp_path / "sexcat.fits", x, y, seg)
    # Row 2 (NUMBER 2) has no partner; rows 1 and 3 become 40 and 70.
    ext = _external(tmp_path / "ext.cat", [70, 40], [x[2], x[0]],
                    [y[2], y[0]])
    mc.match_catalogue(cat, ext, min_fraction=0.5)
    data = _objects(cat)
    assert list(data["NUMBER"]) == [40, 70]
    lut = {0: 0, 1: 40, 2: mc.UNMATCHED_LABEL, 3: 70}
    expected = np.vectorize(lut.get)(seg[[0, 2]])
    np.testing.assert_array_equal(
        np.asarray(data["SEG_VIGNET"]).reshape(2, 3, 3), expected)
    # Each object's own centre label is its new NUMBER, as UberSeg needs.
    centre = np.asarray(data["SEG_VIGNET"]).reshape(2, 3, 3)[:, 1, 1]
    assert list(centre) == [40, 70]


def test_relabel_rejects_unknown_negative_labels():
    with pytest.raises(ValueError):
        mc.relabel(np.array([[-1, 1]]), [1], [5])


def test_match_catalogue_key_expands_to_empty_when_unset(monkeypatch):
    """config_tile_Sx.ini's ${SP_MATCH_CATALOGUE:-} means "no join" unset."""
    parser = CustomParser()
    parser.read_string("[S]\nMATCH = ${SP_MATCH_CATALOGUE:-}\n"
                       "BARE = $SP_MATCH_CATALOGUE\n")
    monkeypatch.delenv("SP_MATCH_CATALOGUE", raising=False)
    assert parser.getexpanded("S", "MATCH") == ""
    with pytest.raises(ValueError):
        parser.getexpanded("S", "BARE")
    monkeypatch.setenv("SP_MATCH_CATALOGUE", "/a/CFIS_cat-186-307.cat")
    assert parser.getexpanded("S", "MATCH") == "/a/CFIS_cat-186-307.cat"
