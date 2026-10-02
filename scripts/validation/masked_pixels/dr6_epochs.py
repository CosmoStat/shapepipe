"""Cut real DR6 epoch stamps into the arrays prepare_ngmix_weights takes,
and find epochs with each kind of defect.

Inputs:
  * one MegaCam exposure with its weight and flag maps (vos:cfis/pitcairn,
    vos:cfis/weights, vos:cfis/flags), e.g. 2086786p.fits.fz;
  * the DR6 tile catalogue with VIGNET (-1e30 neighbour and off-tile
    markers) and SEG_VIGNET columns, as read_ext_sexcat writes it from
    CFIS.<tile>.r.cat and CFIS.<tile>.r.seg.fits.fz
    (sexcat-<xxx>-<yyy>.fits).

The epoch stamps are cut as vignetmaker does: 51x51 around the rounded
pixel of the object's (XWIN_WORLD, YWIN_WORLD) in the CCD's own WCS, CCDs
numbered as split_exp names them (HDU index - 1). The tile VIGNET and seg
stamp are MegaCam-flipped onto the epoch and split into neighbour and
off-tile pixels (split_tile_markers). The background is a sigma-clipped
median over the clean pixels of a 257x257 box around the object, and the
background RMS is constant over the stamp; the pipeline takes both from
SExtractor's maps instead. Only the survey flag and weight maps mark
defects here: the pipeline's mask step adds star halos and borders.

Run from the shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src python scripts/validation/masked_pixels/dr6_epochs.py \
      scan EXP WEIGHT FLAG SEXCAT > scan.txt    # one line per epoch
  PYTHONPATH=src python scripts/validation/masked_pixels/dr6_epochs.py \
      cut EXP WEIGHT FLAG SEXCAT NUMBER OUT.npz
"""
import re
import sys

import fitsio
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from scipy.ndimage import find_objects, label

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
                                 "MAG_AUTO", "FLUX_RADIUS", "VIGNET",
                                 "SEG_VIGNET"])


def tile_name(sexcat):
    return ".".join(re.search(r"(\d{3})-(\d{3})", sexcat).groups())


def ccd_positions(hdus, ra, dec):
    """For each CCD extension, the tile rows whose stamps fit on it and their
    0-indexed [row, col] positions."""
    out = []
    for ext in range(1, len(hdus)):
        h = hdus[ext].header
        x, y = WCS(h).all_world2pix(ra, dec, 1)
        inside = ((x - 1 >= RAD) & (x - 1 < h["NAXIS1"] - RAD)
                  & (y - 1 >= RAD) & (y - 1 < h["NAXIS2"] - RAD))
        rows = np.flatnonzero(inside)
        out.append((ext, rows, np.column_stack([y[rows] - 1, x[rows] - 1])))
    return out


def tile_overlay(tile, i, ccd):
    tile_vign = Ngmix.MegaCamFlip(np.copy(tile["VIGNET"][i]), ccd)
    seg = Ngmix.MegaCamFlip(np.copy(tile["SEG_VIGNET"][i]), ccd)
    neighbour, off_tile = split_tile_markers(tile_vign, (N, N))
    return seg.astype(np.int32), neighbour, off_tile


def components(defect, flag, off_tile):
    """Defect components outside the off-tile band, nearest first: (kind,
    distance, height, width, size, flag values)."""
    labels, _ = label(defect & ~off_tile, structure=np.ones((3, 3)))
    found = []
    for k, sl in enumerate(find_objects(labels), start=1):
        part = labels == k
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        size = int(part.sum())
        if h == N:
            kind = {1: "column"}.get(w, f"columns{w}" if w <= 4 else "wide")
        elif size <= 6:
            kind = "pixel"
        elif w <= 4 and h >= 8:
            kind = "trail"
        else:
            kind = "blob"
        values = sorted(set(np.unique(flag[part]).tolist()) - {0})
        found.append((kind, float(_rr[part].min()), h, w, size, values))
    return sorted(found, key=lambda c: c[1])


def cut(hdus, ext, pos, tile, i):
    """One epoch of tile row ``i`` from CCD extension ``ext``."""
    img, wgt, flg = (h[ext].data for h in hdus)
    ccd = ext - 1
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
    seg, neighbour, off_tile = tile_overlay(tile, i, ccd)
    stamps["flag"][off_tile] = OFF_TILE_FLAG
    stamps.update(
        bkg_rms=np.full((N, N), rms), seg=seg, neighbour=neighbour,
        object_number=int(tile["NUMBER"][i]), ccd=ccd,
        offset=offset[0], int_pos=int_pos[0],
    )
    return stamps


def scan(hdus, tile):
    """Print one line per epoch with any defect pixel."""
    for ext, rows, pos in ccd_positions(hdus[0], tile["XWIN_WORLD"],
                                        tile["YWIN_WORLD"]):
        if not len(rows):
            continue
        wgt = get_stamps(hdus[1][ext].data, pos, RAD)[0]
        flg = get_stamps(hdus[2][ext].data, pos, RAD)[0].astype(np.int32)
        for i, w, f in zip(rows, wgt, flg):
            seg, neighbour, off_tile = tile_overlay(tile, i, ext - 1)
            f[off_tile] = OFF_TILE_FLAG
            defect = defect_mask(w, f)
            if not defect.any():
                continue
            d_nb = _rr[neighbour].min() if neighbour.any() else np.inf
            parts = ";".join(
                f"{k}@{d:.1f}:{h}x{wd}:{n}:{'/'.join(map(str, v))}"
                for k, d, h, wd, n, v in components(defect, f, off_tile)
            )
            print(f"{tile['NUMBER'][i]} {tile['MAG_AUTO'][i]:.2f} "
                  f"{tile['FLUX_RADIUS'][i]:.2f} {ext - 1} "
                  f"{int(defect.sum())} {int(off_tile.sum())} {d_nb:.1f} "
                  f"{parts or '-'}", flush=True)


def main():
    mode, exp, wgt, flg, sexcat = sys.argv[1:6]
    tile = read_tile(sexcat)
    hdus = [fits.open(p) for p in (exp, wgt, flg)]
    if mode == "scan":
        print("NUMBER MAG_AUTO FLUX_RADIUS CCD N_DEFECT N_OFF_TILE "
              "D_NEIGHBOUR COMPONENTS(kind@dist:hxw:size:flags)")
        scan(hdus, tile)
        return
    number, out = int(sys.argv[6]), sys.argv[7]
    i = int(np.flatnonzero(tile["NUMBER"] == number)[0])
    ext, rows, pos = next(
        (e, r, p) for e, r, p in ccd_positions(
            hdus[0], tile["XWIN_WORLD"][[i]], tile["YWIN_WORLD"][[i]])
        if len(r)
    )
    ep = cut(hdus, ext, pos[0], tile, i)
    expname = exp.split("/")[-1].split(".")[0]
    ep["label"] = f"tile {tile_name(sexcat)}, object {number}"
    ep["source"] = f"exposure {expname}, CCD {ep['ccd']}"
    ep.update(ra=tile["XWIN_WORLD"][i], dec=tile["YWIN_WORLD"][i],
              mag=tile["MAG_AUTO"][i])
    np.savez(out, **ep)
    print(ep["label"], ep["source"])


if __name__ == "__main__":
    main()
