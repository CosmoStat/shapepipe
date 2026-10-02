"""Cut one real DR6 epoch of a blended object into the arrays that
prepare_ngmix_weights takes.

Inputs:
  * one MegaCam exposure with its weight and flag maps (vos:cfis/pitcairn,
    vos:cfis/weights, vos:cfis/flags), e.g. 2086786p.fits.fz;
  * the DR6 tile catalogue with VIGNET (-1e30 neighbour and off-tile
    markers) and SEG_VIGNET columns, as read_ext_sexcat writes it from
    CFIS.<tile>.r.cat and CFIS.<tile>.r.seg.fits.fz.

The epoch stamps are cut as vignetmaker does: 51x51 around the rounded
pixel of the object's (XWIN_WORLD, YWIN_WORLD) in the CCD's own WCS. The
tile VIGNET and seg stamp are MegaCam-flipped onto the epoch and split into
neighbour and off-tile pixels (split_tile_markers). The background is a
sigma-clipped median over the clean pixels of a 257x257 box around the
object, and the background RMS is constant over the stamp; the pipeline
takes both from SExtractor's maps instead. Only the survey flag map marks
defects here: the pipeline's mask step adds star halos and borders.

Run from the shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src python scripts/validation/masked_pixels/dr6_blend.py \
      scan EXP WEIGHT FLAG SEXCAT                  # rank candidate epochs
  PYTHONPATH=src python scripts/validation/masked_pixels/dr6_blend.py \
      cut EXP WEIGHT FLAG SEXCAT NUMBER OUT.npz    # write one epoch
"""
import sys

import fitsio
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

from shapepipe.modules.ngmix_package.ngmix import (
    OFF_TILE_FLAG,
    Ngmix,
    defect_mask,
    split_tile_markers,
)
from shapepipe.modules.vignetmaker_package.vignetmaker import get_stamps

RAD = 25
N = 2 * RAD + 1
_rr = np.hypot(*(np.indices((N, N)) - RAD))
BOX = 128


def read_tile(sexcat):
    with fitsio.FITS(sexcat) as f:
        hdu = [h for h in f if h.get_extname() == "LDAC_OBJECTS"][0]
        return hdu.read(columns=["NUMBER", "XWIN_WORLD", "YWIN_WORLD",
                                 "MAG_AUTO", "VIGNET", "SEG_VIGNET"])


def blend_rows(tile):
    """Rows of objects with one sizeable neighbour footprint in their stamp,
    no off-tile pixels and a target footprint of at least 40 pixels."""
    rows = []
    for i, (num, vign, seg) in enumerate(
            zip(tile["NUMBER"], tile["VIGNET"], tile["SEG_VIGNET"])):
        nb, off = split_tile_markers(vign, (N, N))
        if off.any() or (seg == num).sum() < 40:
            continue
        labels = [lab for lab in set(np.unique(seg[nb])) - {0, num}
                  if (seg == lab).sum() >= 40]
        if len(labels) == 1:
            rows.append(i)
    return np.array(rows)


def ccd_hits(hdus, ra, dec, _wcs={}):
    """(extension, [row, col]) of every CCD of the exposure holding (ra, dec)
    at least RAD pixels from its border."""
    hits = []
    for ext in range(1, len(hdus)):
        h = hdus[ext].header
        if (id(hdus), ext) not in _wcs:
            _wcs[id(hdus), ext] = WCS(h)
        x, y = _wcs[id(hdus), ext].all_world2pix(ra, dec, 1)
        if RAD <= x - 1 < h["NAXIS1"] - RAD and RAD <= y - 1 < h["NAXIS2"] - RAD:
            hits.append((ext, np.array([y - 1, x - 1])))
    return hits


def cut(img_hdus, wgt_hdus, flg_hdus, ext, pos, tile, i):
    """One epoch of tile row ``i`` from CCD extension ``ext``."""
    img, wgt, flg = (h[ext].data for h in (img_hdus, wgt_hdus, flg_hdus))
    ccd = ext - 1  # split_exp names CCD files by HDU index - 1
    stamps = {}
    for key, a, dtype in (("gal", img, float), ("weight", wgt, float),
                          ("flag", flg, np.int32)):
        s, int_pos, offset = get_stamps(a, pos[None, :], RAD)
        stamps[key] = s[0].astype(dtype)
    r0, c0 = int_pos[0]
    box = (slice(max(r0 - BOX, 0), r0 + BOX + 1),
           slice(max(c0 - BOX, 0), c0 + BOX + 1))
    good = (wgt[box] > 0) & (flg[box] == 0)
    _, bkg, rms = sigma_clipped_stats(img[box][good], sigma=3.0)
    stamps["gal"] = stamps["gal"] - bkg
    tile_vign = Ngmix.MegaCamFlip(np.copy(tile["VIGNET"][i]), ccd)
    seg = Ngmix.MegaCamFlip(np.copy(tile["SEG_VIGNET"][i]), ccd)
    neighbour, off_tile = split_tile_markers(tile_vign, (N, N))
    stamps["flag"][off_tile] = OFF_TILE_FLAG
    stamps.update(
        bkg_rms=np.full((N, N), rms), seg=seg.astype(np.int32),
        neighbour=neighbour, object_number=int(tile["NUMBER"][i]),
        ccd=ccd, offset=offset[0], int_pos=int_pos[0],
    )
    return stamps


def describe(ep):
    defect = defect_mask(ep["weight"], ep["flag"], ep["bkg_rms"])
    d_def = _rr[defect].min() if defect.any() else np.inf
    d_nb = _rr[ep["neighbour"]].min() if ep["neighbour"].any() else np.inf
    return defect.sum(), d_def, d_nb


def main():
    mode, exp, wgt, flg, sexcat = sys.argv[1:6]
    tile = read_tile(sexcat)
    img_hdus, wgt_hdus, flg_hdus = (fits.open(p) for p in (exp, wgt, flg))
    if mode == "scan":
        for i in blend_rows(tile):
            for ext, pos in ccd_hits(img_hdus, tile["XWIN_WORLD"][i],
                                     tile["YWIN_WORLD"][i]):
                ep = cut(img_hdus, wgt_hdus, flg_hdus, ext, pos, tile, i)
                n_def, d_def, d_nb = describe(ep)
                print(f"NUMBER {tile['NUMBER'][i]} mag "
                      f"{tile['MAG_AUTO'][i]:.2f} ccd {ep['ccd']} "
                      f"defects {n_def} at >= {d_def:.1f} px, "
                      f"neighbour at >= {d_nb:.1f} px", flush=True)
        return
    number, out = int(sys.argv[6]), sys.argv[7]
    i = int(np.flatnonzero(tile["NUMBER"] == number)[0])
    ext, pos = ccd_hits(img_hdus, tile["XWIN_WORLD"][i],
                        tile["YWIN_WORLD"][i])[0]
    ep = cut(img_hdus, wgt_hdus, flg_hdus, ext, pos, tile, i)
    expname = exp.split("/")[-1].split(".")[0]
    ep["label"] = f"tile 202.301, object {number}"
    ep["source"] = f"exposure {expname}, CCD {ep['ccd']}"
    ep.update(ra=tile["XWIN_WORLD"][i], dec=tile["YWIN_WORLD"][i],
              mag=tile["MAG_AUTO"][i])
    np.savez(out, **ep)
    print(ep["label"], describe(ep))


if __name__ == "__main__":
    main()
