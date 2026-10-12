"""The per-exposure fragment and the campaign sum (``exp_maps.py``,
``merge_exposure_maps.py``).

The fixture is a synthetic store of two CCDs on a rotated TAN WCS at MegaCam's
0.187"/pixel, so a sky pixel covers ~74 CCD pixels as on the sky. Its flag
image carries a one-pixel bad column, a saturated blob, isolated pixels on the
CCD edges, a lattice of hot pixels, and a flagged overscan border outside
DATASEC, which must count for nothing.

Needs healsparse, hpgeom and astropy, so it runs inside the container and
skips outside.
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = [pytest.mark.unions,
              pytest.mark.decision("masking.defect_map_from_flags")]

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
EXP = "2079612"
NY, NX = 320, 240                 # DATASEC
PAD_X, PAD_Y = 8, 6               # overscan columns either side, rows on top
SCALE_DEG = 0.187 / 3600.0
ROTATION_DEG = 23.0


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.remove(str(SCRIPTS))


@pytest.fixture(scope="module")
def maps():
    for dep in ("healsparse", "hpgeom", "astropy"):
        pytest.importorskip(dep)
    return _load("exp_maps")


def _wcs(crval):
    from astropy.wcs import WCS

    theta = np.deg2rad(ROTATION_DEG)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = list(crval)
    wcs.wcs.crpix = [PAD_X + NX / 2 + 0.5, NY / 2 + 0.5]
    wcs.wcs.cd = SCALE_DEG * np.array([[-np.cos(theta), np.sin(theta)],
                                       [np.sin(theta), np.cos(theta)]])
    return wcs


def _flags():
    """The DATASEC flags; ``_split`` adds the overscan border."""
    flags = np.zeros((NY, NX), dtype=np.int16)
    flags[:, 101] = 1                                       # bad column
    yy, xx = np.mgrid[:NY, :NX]
    flags[(yy - 210) ** 2 + (xx - 60) ** 2 <= 7 ** 2] = 2   # saturated blob
    for row, col in [(0, 0), (NY - 1, NX - 1), (40, NX - 1), (77, 33)]:
        flags[row, col] = 8                                 # isolated pixels
    flags[5::11, 3::11] = 8                                 # sparse hot pixels
    return flags


def _split(flags):
    """The full split: DATASEC flags inside a border flagged 3."""
    full = np.full((NY + PAD_Y, NX + 2 * PAD_X), 3, dtype=np.int16)
    full[:NY, PAD_X:PAD_X + NX] = flags
    return full


def _store(tmp_path, ccds, crvals, valid):
    """A split dir with image (header only) and flag splits for ``ccds``, and
    an exp_persist manifest naming the ``valid`` CCDs' validation_psf files."""
    from astropy.io import fits

    split = tmp_path / "store" / "output/run_sp_exp_Sp/split_exp_runner/output"
    split.mkdir(parents=True)
    for ccd, crval in zip(ccds, crvals):
        header = _wcs(crval).to_header()
        header["DATASEC"] = f"[{PAD_X + 1}:{PAD_X + NX},1:{NY}]"
        fits.PrimaryHDU(header=header).writeto(split / f"image-{EXP}-{ccd}.fits")
        fits.PrimaryHDU(data=_split(_flags())).writeto(
            split / f"flag-{EXP}-{ccd}.fits")
    manifest = tmp_path / "exp_persist.json"
    manifest.write_text(json.dumps({"files": [
        {"name": f"validation_psf-{EXP}-{c}.fits"} for c in valid]
        + [{"name": f"psf-{EXP}-0.psf"}]}))
    return tmp_path / "store", manifest


def _run(maps, tmp_path, monkeypatch, **store):
    import healsparse as hsp

    exp_dir, persist = _store(tmp_path, **store)
    fragment = tmp_path / "maps" / f"maps-{EXP}.hsp"
    monkeypatch.setattr(sys, "argv", [
        "exp_maps.py", "--exp-dir", str(exp_dir), "--exp", EXP,
        "--persist-manifest", str(persist), "--fragment", str(fragment),
        "--manifest", str(tmp_path / "exp_maps.json")])
    maps.main()
    return hsp.HealSparseMap.read(str(fragment)), json.loads(
        (tmp_path / "exp_maps.json").read_text())


def _sky_pixels(maps, crval, x, y):
    """Sky pixels of DATASEC pixel positions (0-based within DATASEC)."""
    import hpgeom as hpg

    ra, dec = _wcs(crval).pixel_to_world_values(np.asarray(x) + PAD_X, y)
    return hpg.angle_to_pixel(maps.NSIDE, ra, dec, nest=True)


def test_fragment_counts_every_flagged_pixel(maps, tmp_path, monkeypatch):
    """Each flagged DATASEC pixel is counted once, in the sky pixel holding
    its centre; the overscan border counts for nothing."""
    crval = (150.3, 31.7)
    fragment, record = _run(maps, tmp_path, monkeypatch,
                            ccds=[0], crvals=[crval], valid=[0])
    rows, cols = np.nonzero(_flags())
    expected = dict(zip(*np.unique(_sky_pixels(maps, crval, cols, rows),
                                   return_counts=True)))
    pixels = fragment.valid_pixels
    values = fragment.get_values_pix(pixels).astype(int)
    got = {p: v - 1 for p, v in zip(pixels, values) if v > 1}
    # Centres on the very edge of DATASEC can fall in a sky pixel whose own
    # centre is outside the coverage polygon; those are dropped.
    assert set(got) <= set(expected)
    assert all(got[p] == expected[p] for p in got)
    assert sum(got.values()) > 0.97 * rows.size
    assert record["ccds"]["0"]["flagged_pixels"] == rows.size


def test_coverage_is_the_imaging_area_of_valid_psf_ccds(maps, tmp_path,
                                                        monkeypatch):
    """Interior pixel centres are covered, the overscan is not, and a CCD
    without a PSF model contributes nothing."""
    crvals = [(150.3, 31.7), (150.5, 31.7)]
    fragment, record = _run(maps, tmp_path, monkeypatch,
                            ccds=[0, 1], crvals=crvals, valid=[0])
    yy, xx = np.mgrid[13:NY - 13:7, 13:NX - 13:7]   # > one sky pixel in
    inside = _sky_pixels(maps, crvals[0], xx.ravel(), yy.ravel())
    assert (fragment.get_values_pix(inside) > 0).all()
    overscan = _sky_pixels(maps, crvals[0], np.full(50, -PAD_X - 4.0),
                           np.linspace(20, NY - 20, 50))
    assert (fragment.get_values_pix(overscan) == 0).all()
    other = _sky_pixels(maps, crvals[1], xx.ravel(), yy.ravel())
    assert (fragment.get_values_pix(other) == 0).all()
    assert list(record["ccds"]) == ["0"]


def test_merge_counts_exposures_and_flagged_pixels(maps, tmp_path,
                                                   monkeypatch):
    """Two fragments overlapping on half their pixels: nexp is 2 in the
    overlap, nflagged sums the counts, and a fragment-less exposure is
    skipped."""
    import healsparse as hsp

    merge = _load("merge_exposure_maps")
    products = tmp_path / "products"
    a = np.arange(1000, dtype=np.int64) + 10**9
    for exp, pix in (("2000001", a), ("2000002", a + 500)):
        frag = hsp.HealSparseMap.make_empty(maps.NSIDE_COVERAGE, maps.NSIDE,
                                            np.uint8)
        frag[pix] = np.ones(pix.size, np.uint8)
        frag[pix[500:510]] = np.full(10, 1 + 7, np.uint8)
        path = merge.fragment_path(products, exp)
        path.parent.mkdir(parents=True)
        frag.write(str(path))
    monkeypatch.setattr(merge.build_index, "campaign_exposures",
                        lambda *_: ["2000001", "2000002", "2000003"])
    out = {k: tmp_path / f"{k}.hsp" for k in ("nexp", "nflagged")}
    monkeypatch.setattr(sys, "argv", [
        "merge_exposure_maps.py", "--products-dir", str(products),
        "--tile-list", "-", "--index-db", "-",
        "--nexp", str(out["nexp"]), "--nflagged", str(out["nflagged"])])
    merge.main()

    nexp = hsp.HealSparseMap.read(str(out["nexp"]))
    nflagged = hsp.HealSparseMap.read(str(out["nflagged"]))
    assert nexp.get_values_pix(a[:500]).tolist() == [1] * 500
    assert nexp.get_values_pix(a[500:]).tolist() == [2] * 500
    assert nexp.get_values_pix(a[500:] + 500).tolist() == [1] * 500
    assert nflagged.get_values_pix(a[500:510]).tolist() == [7] * 10
    assert nflagged.get_values_pix(a[500:510] + 500).tolist() == [7] * 10
    assert nflagged.valid_pixels.size == 20
