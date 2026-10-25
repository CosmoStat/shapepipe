"""The per-exposure fragment and the campaign sum (``exp_maps.py``,
``merge_exposure_maps.py``).

The fixture is a synthetic CCD grid on a rotated TAN WCS at MegaCam's
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
              pytest.mark.decision("masking.defect_map_from_flags"),
              pytest.mark.decision("masking.nexp_map_valid_psf_ccds"),
              pytest.mark.decision("masking.pointing_square")]

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
    values = (fragment.get_values_pix(pixels) >> maps.FLAG_SHIFT).astype(int)
    got = {p: v for p, v in zip(pixels, values) if v > 0}
    # Centres on the very edge of DATASEC can fall in a sky pixel whose own
    # centre is outside the coverage polygon; those are dropped.
    assert set(got) <= set(expected)
    assert all(got[p] == expected[p] for p in got)
    assert sum(got.values()) > 0.97 * rows.size
    assert record["ccds"]["0"]["flagged_pixels"] == rows.size


def test_coverage_is_the_imaging_area_of_valid_psf_ccds(maps, tmp_path,
                                                        monkeypatch):
    """Validated CCD coverage excludes overscan; unvalidated CCDs are taken."""
    crvals = [(150.3, 31.7), (150.5, 31.7)]
    fragment, record = _run(maps, tmp_path, monkeypatch,
                            ccds=[0, 1], crvals=crvals, valid=[0])
    yy, xx = np.mgrid[13:NY - 13:7, 13:NX - 13:7]   # > one sky pixel in
    inside = _sky_pixels(maps, crvals[0], xx.ravel(), yy.ravel())
    assert (fragment.get_values_pix(inside) & 4 != 0).all()
    overscan = _sky_pixels(maps, crvals[0], np.full(50, -PAD_X - 4.0),
                           np.linspace(20, NY - 20, 50))
    assert (fragment.get_values_pix(overscan) & 2 == 0).all()
    other = _sky_pixels(maps, crvals[1], xx.ravel(), yy.ravel())
    values = fragment.get_values_pix(other)
    assert (values & 2 != 0).all()
    assert (values & 4 == 0).all()
    assert (values >> maps.FLAG_SHIFT == 0).all()
    assert list(record["ccds"]) == ["0", "1"]


def test_merge_counts_exposures_and_flagged_pixels(maps, tmp_path,
                                                   monkeypatch):
    """Two fragments overlapping on half their pixels: nexp is 2 in the
    overlap and nflagged sums the counts."""
    import healsparse as hsp

    merge = _load("merge_exposure_maps")
    products = tmp_path / "products"
    a = np.arange(1000, dtype=np.int64) + 10**9
    for exp, pix in (("2000001", a), ("2000002", a + 500)):
        frag = hsp.HealSparseMap.make_empty(maps.NSIDE_COVERAGE, maps.NSIDE,
                                            maps.FRAGMENT_DTYPE)
        values = np.full(pix.size, 7, maps.FRAGMENT_DTYPE)
        values[500:510] |= 7 << maps.FLAG_SHIFT
        frag[pix] = values
        path = merge.fragment_path(products, exp)
        path.parent.mkdir(parents=True)
        frag.write(str(path))
    monkeypatch.setattr(merge.build_index, "campaign_exposures",
                        lambda *_: ["2000001", "2000002"])
    out = {k: tmp_path / f"{k}.hsp" for k in
           ("npoint", "ntaken", "nexp", "nflagged")}
    monkeypatch.setattr(sys, "argv", [
        "merge_exposure_maps.py", "--products-dir", str(products),
        "--tile-list", "-", "--index-db", "-",
        *[arg for k, path in out.items() for arg in (f"--{k}", str(path))]])
    merge.main()

    nexp = hsp.HealSparseMap.read(str(out["nexp"]))
    nflagged = hsp.HealSparseMap.read(str(out["nflagged"]))
    assert nexp.get_values_pix(a[:500]).tolist() == [1] * 500
    assert nexp.get_values_pix(a[500:]).tolist() == [2] * 500
    assert nexp.get_values_pix(a[500:] + 500).tolist() == [1] * 500
    assert nflagged.get_values_pix(a[500:510]).tolist() == [7] * 10
    assert nflagged.get_values_pix(a[500:510] + 500).tolist() == [7] * 10
    assert nflagged.valid_pixels.size == 20
    for name in ("npoint", "ntaken"):
        count = hsp.HealSparseMap.read(str(out[name]))
        assert np.array_equal(count.get_values_pix(a), nexp.get_values_pix(a))


def test_pointing_fills_gaps_excludes_ears_and_matches_mef(maps, tmp_path,
                                                          monkeypatch):
    """A main-array grid has gaps, an unvalidated CCD, and an external ear."""
    import healsparse as hsp
    import hpgeom as hpg
    from astropy.io import fits

    # Arrange centres in the rotated WCS's own pixel plane; 80-pixel gaps.
    base = _wcs((359.99, 31.7))
    offsets = [(0, 0), (NX + 80, 0), (0, NY + 80),
               (NX + 80, NY + 80), (2 * (NX + 80), 0)]
    crvals = [base.all_pix2world(base.wcs.crpix[0] + x,
                               base.wcs.crpix[1] + y, 1) for x, y in offsets]
    ccds = [0, 1, 9, 10, 36]
    fragment, _ = _run(maps, tmp_path, monkeypatch,
                       ccds=ccds, crvals=crvals, valid=[0, 9, 10, 36])
    ra, dec = base.all_pix2world(base.wcs.crpix[0] + (NX + 80) / 2,
                                base.wcs.crpix[1], 1)
    gap = hpg.angle_to_pixel(maps.NSIDE, ra, dec)
    assert fragment[gap] & maps.COVER_MASK == 1
    ear = _sky_pixels(maps, crvals[-1], NX // 2, NY // 2)
    assert fragment[ear] & maps.COVER_MASK == 6
    unvalidated = _sky_pixels(maps, crvals[1], NX // 2, NY // 2)
    assert fragment[unvalidated] & maps.COVER_MASK == 3
    assert fragment[unvalidated] >> maps.FLAG_SHIFT == 0

    split = tmp_path / "store" / maps.SPLIT_DIR
    images, flags = [fits.PrimaryHDU()], [fits.PrimaryHDU()]
    for ccd in range(37):
        if ccd in ccds:
            header = fits.getheader(split / f"image-{EXP}-{ccd}.fits")
            images.append(fits.ImageHDU(data=np.zeros(_split(_flags()).shape),
                                        header=header))
            # Delivered floating flags follow split_exp's int16 conversion.
            flags.append(fits.ImageHDU(
                data=_split(_flags()).astype(np.float32) + 0.25))
        else:
            # Raw MEFs index every CCD; put absent grid CCDs at an existing
            # unvalidated footprint so coverage and flags remain unchanged.
            header = fits.getheader(split / f"image-{EXP}-1.fits")
            images.append(fits.ImageHDU(data=np.zeros(_split(_flags()).shape),
                                        header=header))
            flags.append(fits.ImageHDU(
                data=_split(_flags()).astype(np.float32) + 0.25))
    image_path, flag_path = tmp_path / "image.fits", tmp_path / "flag.fits"
    fits.HDUList(images).writeto(image_path)
    fits.HDUList(flags).writeto(flag_path)
    raw_path = tmp_path / "raw.hsp"
    monkeypatch.setattr(sys, "argv", [
        "exp_maps.py", "--image-mef", str(image_path),
        "--flag-mef", str(flag_path),
        "--exp", EXP, "--persist-manifest", str(tmp_path / "exp_persist.json"),
        "--fragment", str(raw_path), "--manifest", str(tmp_path / "raw.json")])
    maps.main()
    raw = hsp.HealSparseMap.read(str(raw_path))
    assert np.array_equal(raw.valid_pixels, fragment.valid_pixels)
    assert np.array_equal(raw[raw.valid_pixels],
                          fragment[fragment.valid_pixels])

    # Pointing and taken coverage do not depend on any PSF succeeding.
    manifest = tmp_path / "exp_persist.json"
    manifest.write_text(json.dumps({"files": []}))
    maps.main()
    empty_psf = hsp.HealSparseMap.read(str(raw_path))
    values = empty_psf[raw.valid_pixels]
    assert np.array_equal(values & 3, raw[raw.valid_pixels] & 3)
    assert (values & 4 == 0).all()
    assert (values >> maps.FLAG_SHIFT == 0).all()


def test_fragment_clips_flag_counts_and_preserves_mtime(maps, tmp_path,
                                                       monkeypatch):
    """Flag counts saturate at 255; identical writes preserve mtimes."""
    crval = (150.3, 31.7)
    pixel = _sky_pixels(maps, crval, NX // 2, NY // 2)
    monkeypatch.setattr(maps, "flagged_counts", lambda *_: (
        np.atleast_1d(pixel), np.array([400])))
    fragment, _ = _run(maps, tmp_path, monkeypatch,
                       ccds=[0], crvals=[crval], valid=[0])
    assert fragment[pixel] >> maps.FLAG_SHIFT == 255
    assert fragment[pixel] & maps.COVER_MASK == 7
    paths = [tmp_path / "maps" / f"maps-{EXP}.hsp", tmp_path / "exp_maps.json"]
    mtimes = [p.stat().st_mtime_ns for p in paths]
    maps.main()
    assert [p.stat().st_mtime_ns for p in paths] == mtimes


def test_merge_refuses_plain_uint8(maps, tmp_path, monkeypatch):
    """Legacy fragments require an explicit exposure-specific raw backfill."""
    import healsparse as hsp

    merge = _load("merge_exposure_maps")
    path = merge.fragment_path(tmp_path, EXP)
    path.parent.mkdir(parents=True)
    old = hsp.HealSparseMap.make_empty(
        maps.NSIDE_COVERAGE, maps.NSIDE, np.uint8)
    old[np.array([10**9])] = np.array([1], dtype=np.uint8)
    old.write(str(path))
    monkeypatch.setattr(merge.build_index, "campaign_exposures",
                        lambda *_: [EXP])
    monkeypatch.setattr(sys, "argv", [
        "merge_exposure_maps.py", "--products-dir", str(tmp_path),
        "--tile-list", "-", "--index-db", "-",
        *[arg for name in ("npoint", "ntaken", "nexp", "nflagged")
          for arg in (f"--{name}", str(tmp_path / f"{name}.hsp"))]])
    with pytest.raises(SystemExit, match=f"Exposure {EXP}.*") as error:
        merge.main()
    assert "--image-mef" in str(error.value)
    assert "--persist-manifest" in str(error.value)


def test_merge_refuses_missing_processed_fragments(maps, tmp_path, monkeypatch):
    """One good fragment cannot hide missing processed exposures."""
    import healsparse as hsp

    merge = _load("merge_exposure_maps")
    missing = ["2243882", "2243883"]
    path = merge.fragment_path(tmp_path, EXP)
    path.parent.mkdir(parents=True)
    fragment = hsp.HealSparseMap.make_empty(
        maps.NSIDE_COVERAGE, maps.NSIDE, maps.FRAGMENT_DTYPE)
    fragment[np.array([10**9])] = np.array([7], dtype=maps.FRAGMENT_DTYPE)
    fragment.write(str(path))
    for exp in missing:
        manifest = (merge.fragment_path(tmp_path, exp).parent.parent
                    / "manifests" / "exp_persist.json")
        manifest.parent.mkdir(parents=True)
        manifest.write_text('{}')
    monkeypatch.setattr(merge.build_index, "campaign_exposures",
                        lambda *_: [EXP, *missing])
    outputs = [tmp_path / f"{name}.hsp"
               for name in ("npoint", "ntaken", "nexp", "nflagged")]
    for output in outputs:
        output.write_bytes(b"previous complete map")
    monkeypatch.setattr(sys, "argv", [
        "merge_exposure_maps.py", "--products-dir", str(tmp_path),
        "--tile-list", "-", "--index-db", "-",
        *[arg for out in outputs for arg in (f"--{out.stem}", str(out))]])
    with pytest.raises(SystemExit, match="missing fragments") as error:
        merge.main()
    for exp in missing:
        assert f"--exp {exp} " in str(error.value)
        assert f"{exp}p.flag.fits.fz" in str(error.value)
    assert "--persist-manifest" in str(error.value)
    assert all(out.read_bytes() == b"previous complete map" for out in outputs)


def test_straighten_drops_nearly_collinear_hull_vertices():
    """A hull vertex turning by ~1e-6 rad (CCD corners along one side) is
    dropped, so ``query_polygon`` never sees a degenerate corner; real
    corners stay."""
    pytest.importorskip("healsparse")
    exp_maps = _load("exp_maps")
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    nearly = np.array([[0, 0], [0.5, 1e-7], [1, 0], [1, 1], [0.5, 1 + 1e-6],
                       [0, 1]], dtype=float)
    np.testing.assert_array_equal(exp_maps.straighten(nearly), square)
    np.testing.assert_array_equal(exp_maps.straighten(square), square)
