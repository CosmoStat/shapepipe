"""UNIT TESTS FOR SEXTRACTOR MODE'S SEG_VIGNET.

``sextractor_script.add_seg_vignet`` cuts the SEGMENTATION check image on the
grid SExtractor cut each object's VIGNET on and adds it to the sexcat as the
int32 ``SEG_VIGNET`` column. SExtractor 2.25.0 centres VIGNET on
``(int)(mx + 0.49999)`` of its 0-based double-precision barycentre
(``src/analyse.c``), which the catalogue's float32 X_IMAGE / Y_IMAGE cannot
resolve near half pixels, so the runner asks SExtractor for X_IMAGE_DBL /
Y_IMAGE_DBL too (``seg_vignet_param_file``), and ``add_seg_vignet`` centres
on them and then drops them.

The expected centres below are written out by hand from that rule, and the
expected stamps are cut by plain slicing, not by the module's helpers.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits

from shapepipe.modules.sextractor_package import sextractor_script as ss

NX, NY = 60, 50
STAMP = 7
BIG = np.float32(-1e30)

# (NUMBER, X_IMAGE_DBL, Y_IMAGE_DBL, VIGNET centre col, row (0-based)).
# Object 3's x is an exact half with an odd floor (np.rint goes to 32 - 1),
# object 4's y has a fractional part in [0.5, 0.50001) (np.rint goes up), and
# object 6's x is just past 0.50001, where float32 falls back below it.
OBJECTS = [
    (1, 20.2, 15.7, 19, 15),
    (2, 1.3, 48.9, 0, 48),
    (3, 32.5, 30.1, 31, 29),
    (4, 44.8, 22.5000099, 44, 21),
    (5, 59.9, 1.2, 59, 0),
    (6, 31.5000101, 10.2, 31, 9),
]
DBL = ["X_IMAGE_DBL", "Y_IMAGE_DBL"]


@pytest.mark.parametrize("pos, centre", [
    (32.5, 31),          # exact half, odd floor: down, not to even
    (31.5, 30),          # exact half, even floor
    (31.49999, 30),
    (31.5000099, 30),    # fraction 0.5000099 < 0.50001: down
    (31.5000101, 31),    # fraction 0.5000101 > 0.50001: up
    (31.4999, 30),
    (20.2, 19),
    (20.8, 20),
    (0.6, 0),            # mx = -0.4 truncates to 0
])
def test_vignet_centre_is_sextractors_rule(pos, centre):
    assert ss.vignet_centre([pos]).tolist() == [centre]


def _cut(array, col, row, fill):
    """The STAMP stamp of ``array`` with 0-based (row, col) at its centre."""
    half = STAMP // 2
    padded = np.pad(array, half, constant_values=fill)
    return padded[row:row + STAMP, col:col + STAMP]


def _scene(tmp_path, flat=False):
    """A sexcat whose VIGNETs are cut on the hand-written centres.

    ``flat`` makes the image uniform, so no pixel comparison can tell a tie
    object's two candidate centres apart.
    """
    rng = np.random.default_rng(3)
    image = (np.full((NY, NX), 100, np.float32) if flat
             else rng.normal(100, 5, (NY, NX)).astype(np.float32))
    seg = np.zeros((NY, NX), np.int32)
    for number, _, _, col, row in OBJECTS:
        seg[max(row - 2, 0):row + 3, max(col - 2, 0):col + 3] = number
    number = np.array([o[0] for o in OBJECTS])
    x = np.array([o[1] for o in OBJECTS])
    y = np.array([o[2] for o in OBJECTS])

    vignets = np.array([_cut(image, o[3], o[4], BIG) for o in OBJECTS])
    seg_true = np.array([_cut(seg, o[3], o[4], 0) for o in OBJECTS])
    vignets[(seg_true != 0) & (seg_true != number[:, None, None])] = BIG

    paths = {name: str(tmp_path / f"{name}-001-001.fits")
             for name in ("segmentation", "sexcat")}
    fits.PrimaryHDU(seg).writeto(paths["segmentation"])
    objects = fits.BinTableHDU.from_columns([
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="X_IMAGE", format="E", array=x.astype(np.float32)),
        fits.Column(name="Y_IMAGE", format="E", array=y.astype(np.float32)),
        fits.Column(name="VIGNET", format=f"{STAMP * STAMP}E",
                    array=vignets.reshape(len(number), -1),
                    dim=f"({STAMP},{STAMP})"),
        fits.Column(name="THETA_J2000", format="E",
                    array=np.zeros(len(number))),
        fits.Column(name="X_IMAGE_DBL", format="D", array=x),
        fits.Column(name="Y_IMAGE_DBL", format="D", array=y),
    ], name="LDAC_OBJECTS")
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="80A",
                     array=np.array(["HISTORY x"]))], name="LDAC_IMHEAD")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(
        paths["sexcat"])
    return paths, seg_true


def _add(paths):
    ss.add_seg_vignet(paths["sexcat"], paths["segmentation"])
    with fits.open(paths["sexcat"]) as hdul:
        return [h.name for h in hdul], hdul["LDAC_OBJECTS"].data.copy()


def test_neither_rint_nor_float32_positions_give_the_centres():
    """The premise: np.rint of the double positions misses objects 3 and 4,
    and SExtractor's rule on the float32 positions misses object 6."""
    x = np.array([o[1] for o in OBJECTS])
    y = np.array([o[2] for o in OBJECTS])
    col = np.array([o[3] for o in OBJECTS])
    row = np.array([o[4] for o in OBJECTS])
    rint_miss = (np.rint(x - 1) != col) | (np.rint(y - 1) != row)
    assert np.flatnonzero(rint_miss).tolist() == [2, 3]
    x32 = x.astype(np.float32).astype(np.float64)
    y32 = y.astype(np.float32).astype(np.float64)
    f32_miss = ((ss.vignet_centre(x32) != col)
                | (ss.vignet_centre(y32) != row))
    assert np.flatnonzero(f32_miss).tolist() == [5]


@pytest.mark.parametrize("flat", [False, True], ids=["noisy", "flat"])
def test_seg_vignet_is_registered_with_vignet(tmp_path, flat):
    """SEG_VIGNET is the check image on VIGNET's grid for every object, the
    half-pixel ties included, even where the image cannot tell the two
    candidate centres apart: int32, VIGNET's shape, 0 off the image."""
    paths, seg_true = _scene(tmp_path, flat=flat)
    names, data = _add(paths)
    assert names == ["PRIMARY", "LDAC_IMHEAD", "LDAC_OBJECTS"]
    seg_vignets = data["SEG_VIGNET"]
    assert seg_vignets.dtype.kind == "i" and seg_vignets.dtype.itemsize == 4
    assert seg_vignets.shape == data["VIGNET"].shape
    npt.assert_array_equal(seg_vignets, seg_true)
    centre = STAMP // 2
    npt.assert_array_equal(seg_vignets[:, centre, centre], data["NUMBER"])


def test_seg_vignet_replaces_the_double_positions(tmp_path):
    """The written catalogue has the columns SExtractor writes without the
    double positions, plus SEG_VIGNET, and every other HDU and column as
    they were."""
    paths, _ = _scene(tmp_path)
    with fits.open(paths["sexcat"]) as hdul:
        before = hdul["LDAC_OBJECTS"].data.copy()
        imhead = hdul["LDAC_IMHEAD"].data.copy()
    _, after = _add(paths)
    kept = [name for name in before.names if name not in DBL]
    assert after.names == kept + ["SEG_VIGNET"]
    for name in kept:
        npt.assert_array_equal(after[name], before[name])
    with fits.open(paths["sexcat"]) as hdul:
        npt.assert_array_equal(hdul["LDAC_IMHEAD"].data, imhead)


def test_seg_vignet_needs_the_double_positions(tmp_path):
    paths, _ = _scene(tmp_path)
    with fits.open(paths["sexcat"]) as hdul:
        objects = hdul["LDAC_OBJECTS"]
        cols = [c for c in objects.columns if c.name not in DBL]
        hdus = [hdul[0].copy(), hdul[1].copy(),
                fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS")]
    fits.HDUList(hdus).writeto(paths["sexcat"], overwrite=True)
    with pytest.raises(ValueError, match="X_IMAGE_DBL"):
        _add(paths)


def test_param_file_adds_the_double_positions(tmp_path):
    """The runner's parameter file is the configured one plus the two
    double-precision positions, appended after its last parameter."""
    dot_param = tmp_path / "default.param"
    dot_param.write_text("NUMBER  # running number\nX_IMAGE\nVIGNET(51,51)")
    out = ss.seg_vignet_param_file(str(dot_param), str(tmp_path / "out.param"))
    lines = [line.split("#")[0].strip()
             for line in open(out).read().splitlines()]
    assert [line for line in lines if line] == [
        "NUMBER", "X_IMAGE", "VIGNET(51,51)"] + DBL
    assert dot_param.read_text().count("DBL") == 0


def test_sextractor_caller_names_its_check_images(tmp_path):
    """The runner finds the SEGMENTATION check image by type."""
    caller = ss.SExtractorCaller(
        [str(tmp_path / "image-001-001.fits"),
         str(tmp_path / "weight-001-001.fits")],
        str(tmp_path), "-001-001", "d.sex", "d.param", "d.conv",
        True, False, False, False, False, False, False,
        check_image=["BACKGROUND", "SEGMENTATION"],
    )
    assert caller.check_paths == {
        "BACKGROUND": f"{tmp_path}/background-001-001.fits",
        "SEGMENTATION": f"{tmp_path}/segmentation-001-001.fits",
    }


# --- SEG_VIGNET through the join to the UNIONS catalogue -------------------
#
# UberSeg (uberseg_mask, seg_has_neighbour, Ngmix._check_central_seg_label)
# compares each stamp's labels with the row's NUMBER. The check image is
# labelled with SExtractor's NUMBERs; the join replaces NUMBER with the
# catalogue's, and relabels every stamp through the whole SExtractor -> UNIONS
# map, the footprints of rows it drops included (UNMATCHED_LABEL). The scene
# makes old and new numbers collide: each paired row's new NUMBER is another
# object's SExtractor NUMBER, one of them the dropped detection's.

CROWD_STAMP = 9
# (SExtractor NUMBER, x, y, UNIONS NUMBER or None for no partner)
CROWD = [
    (1, 10.0, 10.0, 2),
    (2, 13.0, 10.0, None),
    (3, 10.0, 13.0, 1),
    (4, 22.0, 22.0, 3),
]
EXT_ONLY = (7, 5.0, 25.0)
EXT_HEADER = (
    "#   1 NUMBER                 Running object number\n"
    "#   2 X_IMAGE                Object position along x  [pixel]\n"
    "#   3 Y_IMAGE                Object position along y  [pixel]\n"
)


def _crowd_seg():
    seg = np.zeros((30, 30), np.int32)
    for number, x, y, _ in CROWD:
        col, row = int(x) - 1, int(y) - 1
        seg[row - 1:row + 2, col - 1:col + 2] = number
    return seg


def _write_crowd(sexcat_path, seg_path):
    """SExtractor's outputs for CROWD: the sexcat (with the double positions
    seg_vignet_param_file asks for) and the SEGMENTATION check image."""
    number = np.array([o[0] for o in CROWD], np.int32)
    x = np.array([o[1] for o in CROWD])
    y = np.array([o[2] for o in CROWD])
    n = len(number)
    size = CROWD_STAMP * CROWD_STAMP
    objects = fits.BinTableHDU.from_columns([
        fits.Column(name="NUMBER", format="J", array=number),
        fits.Column(name="X_IMAGE", format="E", array=x.astype(np.float32)),
        fits.Column(name="Y_IMAGE", format="E", array=y.astype(np.float32)),
        fits.Column(name="VIGNET", format=f"{size}E",
                    array=np.ones((n, size), np.float32),
                    dim=f"({CROWD_STAMP},{CROWD_STAMP})"),
        fits.Column(name="X_IMAGE_DBL", format="D", array=x),
        fits.Column(name="Y_IMAGE_DBL", format="D", array=y),
    ], name="LDAC_OBJECTS")
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="80A",
                     array=np.array(["HISTORY x"]))], name="LDAC_IMHEAD")
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(sexcat_path)
    fits.PrimaryHDU(_crowd_seg()).writeto(seg_path)


def _write_crowd_external(path):
    rows = [(new, x, y) for _, x, y, new in CROWD if new is not None]
    rows.append(EXT_ONLY)
    path.write_text(EXT_HEADER + "".join(
        f"{n:10d} {x:11.4f} {y:11.4f}\n" for n, x, y in rows))
    return str(path)


def _old_stamps():
    """The check image's stamps in SExtractor's numbering, row by row."""
    col = np.array([int(o[1]) - 1 for o in CROWD])
    row = np.array([int(o[2]) - 1 for o in CROWD])
    return ss.cut_stamps(_crowd_seg(), col, row, CROWD_STAMP, 0)


def _run_runner(tmp_path, monkeypatch, match):
    """sextractor_runner as the workflow runs it under blend_handling:
    uberseg, with SExtractor replaced by the CROWD outputs."""
    import logging

    from shapepipe.modules import sextractor_runner as runner_module
    from shapepipe.pipeline.config import CustomParser

    out = tmp_path / "output"
    tmp = tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    num = "-001-001"

    def fake_execute(command_line):
        _write_crowd(out / f"sexcat{num}.fits", out / f"segmentation{num}.fits")
        return "", "All done"

    monkeypatch.setattr(runner_module, "execute", fake_execute)
    dot_param = tmp_path / "default.param"
    dot_param.write_text("NUMBER\nX_IMAGE\nY_IMAGE\nVIGNET(9,9)\n")
    ext = _write_crowd_external(tmp_path / "CFIS_cat-001-001.cat")
    config = CustomParser()
    config["SEXTRACTOR_RUNNER"] = {
        "EXEC_PATH": "source-extractor",
        "DOT_SEX_FILE": "d.sex", "DOT_PARAM_FILE": str(dot_param),
        "DOT_CONV_FILE": "d.conv", "WEIGHT_IMAGE": "True",
        "FLAG_IMAGE": "False", "PSF_FILE": "False",
        "DETECTION_IMAGE": "False", "DETECTION_WEIGHT": "False",
        "ZP_FROM_HEADER": "False", "BKG_FROM_HEADER": "False",
        "CHECKIMAGE": "BACKGROUND, SEGMENTATION",
        "SEG_VIGNET": "${SP_SEG_VIGNET:-False}",
        "MATCH_CATALOGUE": "${SP_MATCH_CATALOGUE:-}",
        "MATCH_RADIUS": "1.0", "MATCH_MIN_FRACTION": "0.98",
        "MATCH_TOLERATED_UNPAIRED": "20", "MAKE_POST_PROCESS": "False",
    }
    monkeypatch.setenv("SP_SEG_VIGNET", "True")
    monkeypatch.setenv("SP_MATCH_CATALOGUE", ext if match else "")
    runner_module.sextractor_runner(
        [str(tmp_path / f"image{num}.fits"),
         str(tmp_path / f"weight{num}.fits")],
        {"output": str(out), "tmp": str(tmp)}, num, config,
        "SEXTRACTOR_RUNNER", logging.getLogger("test"),
    )
    with fits.open(out / f"sexcat{num}.fits") as hdul:
        return hdul["LDAC_OBJECTS"].data.copy()


def test_runner_seg_vignet_carries_the_final_number(tmp_path, monkeypatch):
    """Through the runner, each row's own footprint carries its final
    (UNIONS) NUMBER and nothing else does, though every new NUMBER is some
    other object's SExtractor NUMBER; the dropped detection's footprint is
    UNMATCHED_LABEL. Fails if the join runs before add_seg_vignet, or relabels
    only the paired labels."""
    from shapepipe.modules.sextractor_package import match_catalogue as mc

    data = _run_runner(tmp_path, monkeypatch, match=True)
    paired = [i for i, o in enumerate(CROWD) if o[3] is not None]
    old = _old_stamps()[paired]
    old_number = np.array([CROWD[i][0] for i in paired])
    new_number = np.array([CROWD[i][3] for i in paired])
    npt.assert_array_equal(data["NUMBER"], new_number)
    seg = data["SEG_VIGNET"]
    assert "X_IMAGE_DBL" not in data.names
    centre = CROWD_STAMP // 2
    npt.assert_array_equal(seg[:, centre, centre], data["NUMBER"])
    dropped = CROWD[1][0]
    for i in range(len(paired)):
        npt.assert_array_equal(seg[i] == new_number[i],
                               old[i] == old_number[i])
        npt.assert_array_equal(seg[i] == mc.UNMATCHED_LABEL,
                               old[i] == dropped)
        npt.assert_array_equal(seg[i] == 0, old[i] == 0)
    # Rows 0 and 1 overlap the dropped detection, whose SExtractor label (2)
    # is row 0's new NUMBER.
    assert (seg[0] == mc.UNMATCHED_LABEL).any()
    assert (old[0] == new_number[0]).any()


def test_runner_seg_vignet_without_a_join_keeps_sextractor_numbers(
    tmp_path, monkeypatch,
):
    """Image simulations (empty MATCH_CATALOGUE): every row stays, and the
    stamps are the check image's, labelled with SExtractor's NUMBER."""
    data = _run_runner(tmp_path, monkeypatch, match=False)
    npt.assert_array_equal(data["NUMBER"], [o[0] for o in CROWD])
    npt.assert_array_equal(data["SEG_VIGNET"], _old_stamps())


def test_uberseg_sees_the_same_neighbours_after_the_join(tmp_path,
                                                         monkeypatch):
    """UberSeg only asks own-versus-other, so the relabelled stamps give the
    neighbour mask and flag the SExtractor-labelled ones give."""
    from shapepipe.modules.ngmix_package.ngmix import (
        seg_has_neighbour,
        uberseg_mask,
    )

    data = _run_runner(tmp_path, monkeypatch, match=True)
    paired = [i for i, o in enumerate(CROWD) if o[3] is not None]
    old = _old_stamps()[paired]
    for i, j in enumerate(paired):
        new_seg, number = data["SEG_VIGNET"][i], data["NUMBER"][i]
        npt.assert_array_equal(
            uberseg_mask(new_seg, number),
            uberseg_mask(old[i], CROWD[j][0]))
        assert (seg_has_neighbour(new_seg, number)
                == seg_has_neighbour(old[i], CROWD[j][0]))
    assert seg_has_neighbour(data["SEG_VIGNET"][0], data["NUMBER"][0])
    assert not seg_has_neighbour(data["SEG_VIGNET"][2], data["NUMBER"][2])


def test_a_partial_relabel_would_collide():
    """The negative control: relabelling only the paired labels leaves the
    dropped detection's SExtractor label 2 in row 0's stamp, equal to row 0's
    new NUMBER, so UberSeg would take the neighbour for the object."""
    old = _old_stamps()[0]
    partial = old.copy()
    for number, _, _, new in CROWD:
        if new is not None:
            partial[old == number] = new
    assert ((partial == 2) & (old != 1)).any()


def _join_with_numbers(case, numbers):
    """Join CROWD to a catalogue whose three objects carry ``numbers``
    (written as given), and return the joined NUMBER column."""
    from shapepipe.modules.sextractor_package import match_catalogue as mc

    case.mkdir()
    _write_crowd(case / "sexcat.fits", case / "seg.fits")
    rows = zip(numbers, [10.0, 10.0, 22.0], [10.0, 13.0, 22.0])
    ext = case / "ext.cat"
    ext.write_text(EXT_HEADER + "".join(
        f"{n:>10} {x:11.4f} {y:11.4f}\n" for n, x, y in rows))
    mc.match_catalogue(str(case / "sexcat.fits"), str(ext))
    with fits.open(case / "sexcat.fits") as hdul:
        return hdul["LDAC_OBJECTS"].data["NUMBER"].tolist()


@pytest.mark.parametrize("numbers", [
    [5, 5, 6],                   # repeated
    [0, 5, 6],                   # not positive
    [-3, 5, 6],
    ["5.5", 7, 6],               # not an integer
    ["5.0", 7, 6],               # read as a float
    [2**31, 5, 6],               # past int32
    [2**32 - 1, 5, 6],           # would narrow to -1, UNMATCHED_LABEL
])
def test_join_rejects_numbers_that_could_collide(tmp_path, numbers):
    """Unique catalogue NUMBERs that the int32 columns hold exactly are what
    keep a relabelled stamp's own footprint apart from every other one."""
    with pytest.raises(ValueError, match="NUMBER"):
        _join_with_numbers(tmp_path / "case", numbers)


def test_join_accepts_the_largest_int32_number(tmp_path):
    from shapepipe.modules.sextractor_package import match_catalogue as mc

    assert _join_with_numbers(tmp_path / "case", [1, mc.MAX_NUMBER, 6]) == [
        1, mc.MAX_NUMBER, 6]


def test_seg_vignet_stays_out_of_the_final_catalogue(tmp_path):
    """make_cat drops SEG_VIGNET with VIGNET: the stamps are ngmix inputs,
    not catalogue columns."""
    from shapepipe.modules.make_cat_package import make_cat

    _write_crowd(tmp_path / "sexcat-001-001.fits", tmp_path / "seg.fits")
    ss.add_seg_vignet(str(tmp_path / "sexcat-001-001.fits"),
                      str(tmp_path / "seg.fits"))
    names = list(
        make_cat.read_sextractor_data(str(tmp_path / "sexcat-001-001.fits"))
    )
    assert "VIGNET" not in names and "SEG_VIGNET" not in names
    assert "NUMBER" in names
