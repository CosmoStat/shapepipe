"""UNIT TESTS FOR MODULE PACKAGE: MASK_EXT.

Drives ``MaskExt`` against a synthetic healsparse map and a synthetic TAN WCS
to lock in the rasterization contract of ``mask_ext_runner``: the config-driven
bit->flag mapping (with bitwise-OR of multiple bits), chunk-seam independence,
RA wrap across 0/360, off-footprint (sentinel) handling, external-flag summing,
``int16`` output dtype, and WCS round-trip in the written FITS.

The synthetic map is a high-resolution healsparse map covering only a tiny
patch on the sky; the synthetic WCS points an image at that patch so a subset
of pixels land on masked healpix cells and the rest fall off the footprint.
"""

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from astropy.wcs import WCS

healsparse = pytest.importorskip("healsparse")
hpgeom = pytest.importorskip("hpgeom")

from shapepipe.modules.mask_ext_package.mask_ext import MaskExt


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


# Sky patch the synthetic image and mask share.
CRVAL1 = 150.0
CRVAL2 = 2.3
NSIDE_SPARSE = 131072  # ~1.6 arcsec, matching the real UNIONS masks
NSIDE_COVERAGE = 32
PIXSCALE_DEG = 0.187 / 3600.0  # UNIONS ~0.187"/px


def _make_wcs(naxis1, naxis2, crval1=CRVAL1, crval2=CRVAL2):
    """Synthetic TAN WCS header centred on the shared patch."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [crval1, crval2]
    w.wcs.crpix = [naxis1 / 2 + 0.5, naxis2 / 2 + 0.5]
    w.wcs.cd = [[-PIXSCALE_DEG, 0.0], [0.0, PIXSCALE_DEG]]
    return w


def _write_image(path, wcs_obj, naxis1, naxis2):
    """Write a zero-valued image carrying the given WCS (defines pixel grid)."""
    data = np.zeros((naxis2, naxis1), dtype=np.float32)
    hdu = fits.PrimaryHDU(data)
    # NAXIS* are set from the data shape; merge only the WCS keywords.
    hdu.header.update(wcs_obj.to_header())
    hdu.writeto(path, overwrite=True)


def _make_map(masked_ra, masked_dec, bits, sentinel=0):
    """Healsparse int32 map with the given (ra, dec) cells set to ``bits``."""
    hmap = healsparse.HealSparseMap.make_empty(
        NSIDE_COVERAGE, NSIDE_SPARSE, dtype=np.int32, sentinel=sentinel
    )
    pix = hpgeom.angle_to_pixel(NSIDE_SPARSE, masked_ra, masked_dec)
    bits = np.asarray(bits, dtype=np.int32)
    # Several fine image pixels can share one healpix cell; keep the first bit
    # value per unique cell (replace requires unique pixels).
    pix, first = np.unique(pix, return_index=True)
    hmap[pix] = bits[first]
    return hmap


def _run(tmp_path, image_wcs, naxis1, naxis2, hmap, bit_flag_map,
         off_map_flag=0, chunk_size=None, path_external_flag=None):
    """Instantiate MaskExt on written synthetic inputs and rasterize/write."""
    image_path = str(tmp_path / "image.fits")
    mask_path = str(tmp_path / "mask.hsp")
    _write_image(image_path, image_wcs, naxis1, naxis2)
    hmap.write(mask_path, clobber=True)

    inst = MaskExt(
        image_path,
        mask_path,
        bit_flag_map,
        image_num="-000-000",
        output_dir=str(tmp_path),
        w_log=_NullLogger(),
        off_map_flag=off_map_flag,
        path_external_flag=path_external_flag,
        image_prefix="pipeline",
        chunk_size=chunk_size,
    )
    return inst


def test_parse_bit_flag_map():
    """String config parses into an int->int dict; malformed entries raise."""
    assert MaskExt.parse_bit_flag_map("64:1, 2048:2") == {64: 1, 2048: 2}
    assert MaskExt.parse_bit_flag_map(" 64 : 1 ") == {64: 1}
    with pytest.raises(ValueError):
        MaskExt.parse_bit_flag_map("64")
    with pytest.raises(ValueError):
        MaskExt.parse_bit_flag_map("")


def test_bit_flag_mapping_and_off_map(tmp_path):
    """Masked cells map their bit; unmasked (off-footprint) get off_map_flag.

    The image centre pixel lands on a cell carrying bit 64 -> flag 1; the map
    covers only that centre cell, so every other pixel is off-footprint.
    """
    naxis1 = naxis2 = 16
    w = _make_wcs(naxis1, naxis2)
    # RA/Dec at the central pixel (0-based grid, origin=0).
    cx, cy = naxis1 // 2, naxis2 // 2
    ra_c, dec_c = w.all_pix2world(cx, cy, 0)
    hmap = _make_map([float(ra_c)], [float(dec_c)], [64])

    inst = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1}, off_map_flag=8)
    flags = inst.rasterize()

    assert flags.dtype == np.int16
    assert flags.shape == (naxis2, naxis1)
    # The centre cell is flagged 1; the rest of the image is off the footprint.
    n_flagged_1 = np.sum(flags == 1)
    assert n_flagged_1 >= 1
    assert np.all(flags[flags != 1] == 8)
    assert np.sum(flags == 8) == flags.size - n_flagged_1


def test_multi_bit_or(tmp_path):
    """A cell carrying two bits gets the bitwise-OR of the mapped flags."""
    naxis1 = naxis2 = 8
    w = _make_wcs(naxis1, naxis2)
    cx, cy = naxis1 // 2, naxis2 // 2
    ra_c, dec_c = w.all_pix2world(cx, cy, 0)
    # bit values 64 and 2048 set together on the centre cell.
    hmap = _make_map([float(ra_c)], [float(dec_c)], [64 | 2048])

    inst = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1, 2048: 2})
    flags = inst.rasterize()
    # 1 | 2 == 3 on the masked cell.
    assert 3 in np.unique(flags)
    assert np.all(np.isin(np.unique(flags), [0, 3]))


def test_chunk_seam_independence(tmp_path):
    """Result is independent of chunk size (no seam artifacts)."""
    naxis1 = naxis2 = 40
    w = _make_wcs(naxis1, naxis2)
    # Mask a band of cells spanning several image rows.
    ys = np.arange(0, naxis2)
    xs = np.full_like(ys, naxis1 // 2)
    ra, dec = w.all_pix2world(xs, ys, 0)
    hmap = _make_map(ra.astype(float), dec.astype(float), [64] * len(ys))

    inst_full = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1})
    full = inst_full.rasterize()

    for cs in (1, 3, 7, 40):
        inst = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1}, chunk_size=cs)
        npt.assert_array_equal(inst.rasterize(), full)
    # Sanity: some pixels actually got flagged.
    assert np.any(full == 1)


def test_ra_wrap(tmp_path):
    """RA-wrap near 0/360: an image centred on RA~0 rasterizes correctly."""
    naxis1 = naxis2 = 16
    w = _make_wcs(naxis1, naxis2, crval1=0.0, crval2=2.3)
    cx, cy = naxis1 // 2, naxis2 // 2
    ra_c, dec_c = w.all_pix2world(cx, cy, 0)
    # Straddling pixels produce RA both just below 360 and just above 0; the
    # map cell is queried by its wrapped-into-[0,360) coordinate.
    hmap = _make_map([float(np.mod(ra_c, 360.0))], [float(dec_c)], [64])

    inst = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1})
    flags = inst.rasterize()
    # No crash, correct dtype, and the centre is flagged despite the wrap.
    assert flags.dtype == np.int16
    assert np.any(flags == 1)


def test_external_flag_summing(tmp_path):
    """External instrument flag image is summed into the rasterized mask."""
    naxis1 = naxis2 = 8
    w = _make_wcs(naxis1, naxis2)
    cx, cy = naxis1 // 2, naxis2 // 2
    ra_c, dec_c = w.all_pix2world(cx, cy, 0)
    hmap = _make_map([float(ra_c)], [float(dec_c)], [64])

    # External flag: a constant field of 16 (e.g. a saturation bit).
    ext_path = str(tmp_path / "ext_flag.fits")
    ext_data = np.full((naxis2, naxis1), 16, dtype=np.int16)
    fits.PrimaryHDU(ext_data).writeto(ext_path, overwrite=True)

    inst = _run(
        tmp_path, w, naxis1, naxis2, hmap, {64: 1},
        path_external_flag=ext_path,
    )
    out_path = inst.make_mask()

    with fits.open(out_path) as hdul:
        written = hdul[0].data

    # FITS stores big-endian; still a 2-byte signed int (int16).
    assert written.dtype.kind == "i" and written.dtype.itemsize == 2
    # Every pixel gets +16 from the external flag; the centre cell also +1.
    assert np.all(written >= 16)
    assert 17 in np.unique(written)


def test_write_wcs_roundtrip(tmp_path):
    """Written FITS carries the image WCS; pix2world round-trips."""
    naxis1 = naxis2 = 12
    w = _make_wcs(naxis1, naxis2)
    cx, cy = naxis1 // 2, naxis2 // 2
    ra_c, dec_c = w.all_pix2world(cx, cy, 0)
    hmap = _make_map([float(ra_c)], [float(dec_c)], [64])

    inst = _run(tmp_path, w, naxis1, naxis2, hmap, {64: 1})
    out_path = inst.make_mask()

    with fits.open(out_path) as hdul:
        assert hdul[0].data.dtype.kind == "i"
        assert hdul[0].data.dtype.itemsize == 2
        w_out = WCS(hdul[0].header)

    ra_in, dec_in = w.all_pix2world(cx, cy, 0)
    ra_out, dec_out = w_out.all_pix2world(cx, cy, 0)
    npt.assert_allclose([ra_in, dec_in], [ra_out, dec_out], rtol=0, atol=1e-9)
