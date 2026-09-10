#!/usr/bin/env python3
"""Rasterize ONE exposure's instrument flags into a boolean healsparse fragment.

Run as the shell of the in-DAG ``exp_defect_map`` rule, never by hand.

WHY THIS EXISTS (CosmoStat/shapepipe#878). The survey footprint is built as a
positive coverage map minus the healsparse mask bits, and every masking input
the campaign has is already sky-fixed and queried per object — except one. The
instrument flag image delivered with each exposure (bad columns, saturated
pixels, bleed trails) never leaves the PIXEL domain: ``exp_split`` splits it per
CCD, SExtractor reads it as ``IMAFLAGS_ISO``, and that is the end of it. The
coverage map is built from the CCD corner WCS in the headers, so it cannot
subtract those pixels and the footprint silently includes them. The lost area is
percent-level, but it is exactly the thin, small-scale geometry an accurate
window function needs. Only ShapePipe ever opens these files, so the map has to
come from here.

WHAT IT WRITES. ``<dest>/defect-<exp>.hsp``: a boolean ``HealSparseMap``,
``nside_sparse`` 131072 over ``nside_coverage`` 128, ``bit_packed``, ``True``
where any flag bit is set. Those are not free choices — they are the
convention every other map in this campaign's ladder already follows (the
UNIONS ugriz bit maps under ``inputs.masks``, and ``config_tile_Mc.ini``'s
``MASK_EXT_PATHS``), so the fragments and the map they union into drop into that
ladder without a resolution change. ``True = masked`` likewise.

ANY BIT, NOT A BIT TABLE. The flag image is a bitmask, but the CFIS flags are
"this pixel is not to be trusted" in several flavours (the campaign's own
exposures carry values 1, 2, 3, 8 and 11), and nothing downstream distinguishes
them: SExtractor's ``IMAFLAGS_ISO`` cut is nonzero-vs-zero. So the fragment is
``!= 0`` and the map is boolean. A per-bit ladder would be a different product
answering a question nobody has asked yet.

WHERE THE WCS COMES FROM, AND WHY NOT FROM THE FLAG FILE. The flag mosaic's
per-CCD HDUs carry NO WCS at all — checked on a real exposure: ``CTYPE`` empty,
``CRVAL`` 0, the identity transform. Only the image mosaic is astrometric, so
``split_exp`` saves headers for the image suffix alone. The fragment therefore
takes each CCD's WCS from its ``image-<num>-<ccd>.fits`` split, which sits
beside the flag split and carries the full SCAMP header (``RA---TAN`` with PV
distortion). The image PIXELS are never read: only ``fits.getheader``.

NOT ``headers-<num>.npy``, which is the other thing ``exp_split`` writes and
would be one small file instead of forty header reads. It is a pickled object
array of ``astropy.wcs.WCS`` INSTANCES, so reading it unpickles astropy objects
across whatever container rebuild happens next; the FITS headers are
self-describing text and cost milliseconds. ``merge_headers`` may live with the
pickle because it is the file's author's own consumer; a second consumer should
not inherit the coupling.

HOW A PIXEL IS RASTERIZED, AND WHY IT IS SAMPLED RATHER THAN INTEGRATED. A
healpix pixel at nside 131072 is 1.61 arcsec across; a MegaCam pixel is 0.187
arcsec. So one healpix pixel covers ~74 CCD pixels, and the question is never
"which healpix pixels does this CCD pixel cover" but "which healpix pixels does
the flagged REGION touch". The answer is taken by sampling each flagged CCD
pixel on an ``oversample`` x ``oversample`` grid spanning its full extent
(``linspace(-0.5, 0.5)``, so the four corners are always sampled), converting to
sky through that CCD's WCS and binning with ``ang2pix``.

THE RASTERIZATION IS CONSERVATIVE, deliberately: a healpix pixel is masked if
any part of the flagged region falls in it. The alternative — mask when the
healpix pixel's CENTRE is flagged — erases a one-pixel bad column entirely,
which is precisely the geometry #878 exists to keep. The cost is that a thin
defect is widened to the healpix resolution; at 1.61 arcsec that is the price
of the ladder's nside.

WHAT ``oversample`` BUYS, MEASURED (exposure 2079612p, CCD 0, 430886 flagged
pixels; the reference is a 12x12 interior grid, 14229 healpix pixels):

    oversample   samples/pixel   healpix pixels   missed vs reference
    2 (corners)        4             14081             169  (1.2%)
    3                  9             14198              54  (0.4%)
    5                 25             14254               5  (0.04%)

and the cost is linear in the sample count. 3 is the default in ``config.yaml``
for that reason, and it rides on the rule's ``params`` so raising it re-rasterizes
without touching anything upstream. 2 is the geometric floor: the four corners
of a CCD pixel bound its footprint, and a healpix pixel 74 times its area cannot
sit inside it — the residual 1.2% is boundary and rounding, not a class of
missed defect.

BYTE-STABLE, tmp-then-``cmp``-then-``mv``, the pattern ``persist_exp`` uses:
healsparse's FITS output carries no timestamp (checked), so re-rasterizing an
unchanged store produces an identical file and leaves its mtime alone. mtime is
a rerun trigger and ``clean_exposure`` waits on this rule's manifest, so an
unconditional rewrite would make every reclamation look out of date once per
invocation.

THE MANIFEST IS THE ONLY DECLARED OUTPUT and it lives on the PERSISTENT root
beside the fragment (``<products_dir>/exp/<shard>/<exp>/manifests/``), not in
the exposure's scratch ``manifests/`` which ``clean_exposure`` deletes wholesale
— same placement, and same reason, as ``exp_persist``. It records the per-CCD
healpix counts, so a reader can see which CCD contributed what without opening
the map.

ORDERED BEFORE RECLAMATION. ``clean_exposure`` takes this manifest as an input,
exactly as it takes ``exp_persist``'s: the flag splits live on /scratch and go
with the store, so the fragment must be on /project before anything is deleted.
"""

import argparse
import filecmp
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import healpy as hp
import healsparse as hsp

# The split stage's run dir (RUN_NAME in config_exp_Sp.ini) and its module.
# Hardcoded for the same reason persist_exp.py hardcodes its own: this rule
# rasterizes the SPLIT stage's flag images and nothing else, and a knob here
# would be a knob for "rasterize some other stage".
RUN_NAME = "run_sp_exp_Sp"
MODULE = "split_exp_runner"

# `flag-2079612-13.fits` -> ccd 13. The number string is the exposure's
# ($SP_UNIT_NUM, `-2079612`), so the CCD is what follows the last dash.
_CCD = re.compile(r"-(\d+)\.fits$")


def split_dir(exp_dir: Path) -> Path:
    return exp_dir / "output" / RUN_NAME / MODULE / "output"


def ccd_files(exp_dir: Path) -> list:
    """``(ccd, flag path, image path)`` for every CCD this exposure split, in
    CCD order.

    Driven by the FLAG files, which are what we rasterize; the image split is
    looked up beside each one for its header alone. A flag split with no image
    beside it is a hard error rather than a skip — ``exp_split`` writes the
    three suffixes in one pass, so a missing image means the store is not what
    we think it is, and a fragment quietly short of a CCD would be a hole in the
    footprint nothing downstream could notice.
    """
    root = split_dir(exp_dir)
    out = []
    for flag in sorted(root.glob("flag-*.fits")):
        match = _CCD.search(flag.name)
        if not match:
            continue
        image = flag.with_name(flag.name.replace("flag-", "image-", 1))
        if not image.exists():
            sys.exit(f"defect_map_exp: {flag.name} has no {image.name} beside "
                     f"it in {root}; the WCS lives on the image split (see the "
                     f"module docstring)")
        out.append((int(match.group(1)), flag, image))
    return sorted(out)


def offsets(oversample: int) -> tuple:
    """Sample offsets within one CCD pixel, in pixel units.

    ``linspace`` with both endpoints, so the CORNERS are always sampled: they
    are what bounds the pixel's footprint, and the interior samples only fill
    boundary and rounding gaps (the docstring's table measures how many).
    """
    if oversample < 2:
        sys.exit(f"defect_map_exp: oversample={oversample} would sample the "
                 f"pixel centre alone and lose the pixel's extent; 2 is the "
                 f"geometric floor (its four corners)")
    step = np.linspace(-0.5, 0.5, oversample)
    grid_x, grid_y = np.meshgrid(step, step)
    return grid_x.ravel(), grid_y.ravel()


def rasterize_ccd(flag_path: Path, image_path: Path, nside: int,
                  off_x, off_y) -> np.ndarray:
    """The healpix pixel ids (NEST, ``nside``) this CCD's flags touch.

    One CCD at a time and one array at a time: the whole point of the loop in
    ``main`` is that the job's footprint is one CCD's samples, not the
    exposure's.
    """
    with warnings.catch_warnings():
        # SCAMP headers carry a deprecated RADECSYS and a redundant SIP block
        # beside the PV distortion astropy actually uses; both are FITSFixedWarning
        # noise on every one of 40 CCDs and neither changes the transform.
        warnings.simplefilter("ignore")
        wcs = WCS(fits.getheader(image_path))
    data = fits.getdata(flag_path)
    rows, cols = np.nonzero(data)
    if rows.size == 0:
        return np.empty(0, dtype=np.int64)
    # 1-based FITS pixel coordinates, sampled across each flagged pixel's extent.
    x = (cols[:, None] + 1.0 + off_x[None, :]).ravel()
    y = (rows[:, None] + 1.0 + off_y[None, :]).ravel()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ra, dec = wcs.all_pix2world(x, y, 1)
    return np.unique(hp.ang2pix(nside, ra, dec, lonlat=True, nest=True))


def write_stable(tmp: Path, dest: Path) -> None:
    """Move ``tmp`` onto ``dest``, or drop it when the bytes already match."""
    if dest.exists() and filecmp.cmp(tmp, dest, shallow=False):
        tmp.unlink()
    else:
        tmp.replace(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp-dir", required=True, type=Path,
                        help="the exposure's scratch store")
    parser.add_argument("--exp", required=True)
    parser.add_argument("--dest", required=True, type=Path,
                        help="<products_dir>/exp/<shard>/<exp>/defect; the "
                             "fragment is <dest>/defect-<exp>.hsp")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--nside", type=int, default=131072,
                        help="nside_sparse; the mask ladder's resolution")
    parser.add_argument("--nside-coverage", type=int, default=128)
    parser.add_argument("--oversample", type=int, default=3,
                        help="samples per CCD pixel per axis (see the module "
                             "docstring's measured table)")
    args = parser.parse_args()

    ccds = ccd_files(args.exp_dir)
    if not ccds:
        sys.exit(f"defect_map_exp: {args.exp}: no flag split under "
                 f"{split_dir(args.exp_dir)}")

    off_x, off_y = offsets(args.oversample)
    fragment = hsp.HealSparseMap.make_empty(
        args.nside_coverage, args.nside, np.bool_, bit_packed=True)
    per_ccd = {}
    for ccd, flag_path, image_path in ccds:
        pixels = rasterize_ccd(flag_path, image_path, args.nside, off_x, off_y)
        per_ccd[str(ccd)] = int(pixels.size)
        if pixels.size:
            fragment[pixels] = True

    args.dest.mkdir(parents=True, exist_ok=True)
    frag_path = args.dest / f"defect-{args.exp}.hsp"
    tmp = frag_path.with_name(frag_path.name + ".tmp")
    try:
        fragment.write(str(tmp), clobber=True)
        write_stable(tmp, frag_path)
    finally:
        tmp.unlink(missing_ok=True)

    body = {
        "stage": "exp_defect_map", "level": "exp", "unit": args.exp,
        "status": "complete",
        "map": str(frag_path),
        "nside": args.nside,
        "nside_coverage": args.nside_coverage,
        "oversample": args.oversample,
        "n_ccds": len(ccds),
        # Per CCD, so a reader can see WHICH CCD contributed what without
        # opening the map — a CCD at zero is a real thing (a clean chip) and a
        # whole exposure at zero is not.
        "pixels_per_ccd": per_ccd,
        "n_pixels": int(fragment.n_valid),
        "n_coverage_pixels": int(fragment.coverage_mask.sum()),
        "bytes": frag_path.stat().st_size,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.manifest.with_name(args.manifest.name + ".tmp")
    try:
        tmp.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
        write_stable(tmp, args.manifest)
    finally:
        tmp.unlink(missing_ok=True)

    print(f"[defect_map_exp] {args.exp}: {len(ccds)} CCD(s), "
          f"{body['n_pixels']} healpix pixel(s) over "
          f"{body['n_coverage_pixels']} coverage pixel(s), "
          f"{body['bytes'] / 1e6:.1f} MB -> {frag_path}")


if __name__ == "__main__":
    main()
