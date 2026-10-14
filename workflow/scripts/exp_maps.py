"""Rasterize one exposure into a packed uint16 HealSparse fragment.

Bits 0–2 store coverage: 1 = pointing hull, 2 = CCD DATASEC,
4 = validated-PSF CCD DATASEC. Bits 3–10 store nonzero flagged CCD pixel
centre counts inside validated CCDs' coverage, clipped at 255.
A scalar map uses HealSparse's compressed image storage; record-array maps
use uncompressed tables and cost about 13 times as much disk space.
The maps use nside 131072 over coverage 128, aligned with the sky masks.

Headers and flag arrays come from either the split scratch store or delivered
image/flag MEFs (CCD + 1 is the HDU index). The exp_persist manifest names the
validated CCDs. Both the fragment and exp_maps manifest are written, including
in raw-MEF backfills; durable_edge reads the fragment for reclaimed stores.
"""

import argparse
import filecmp
import json
import warnings
from contextlib import ExitStack
from fnmatch import fnmatch
from pathlib import Path

import healsparse as hsp
import hpgeom as hpg
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import ConvexHull

import persist_exp

NSIDE = 131072
NSIDE_COVERAGE = 128
FRAGMENT_DTYPE = np.dtype(np.uint16)
COVER_MASK = 7
FLAG_SHIFT = 3
SPLIT_DIR = "output/run_sp_exp_Sp/split_exp_runner/output"
PSF_PATTERN = persist_exp.resolve(persist_exp.ALWAYS)


def valid_ccds(manifest: Path, exp: str) -> list:
    """CCD numbers with a PSF model, from the exp_persist manifest's members."""
    prefix = PSF_PATTERN.split("*")[0] + f"{exp}-"
    names = [f["name"] for f in json.loads(manifest.read_text())["files"]]
    return sorted(int(n[len(prefix):].removesuffix(".fits")) for n in names
                  if fnmatch(n, PSF_PATTERN) and n.startswith(prefix))


def ccd_wcs(header):
    """Return WCS and DATASEC bounds (1-based inclusive) from a CCD header."""
    if "DATASEC" in header:
        x, y = header["DATASEC"].strip("[]").split(",")
        bounds = [int(v) for v in x.split(":") + y.split(":")]
    else:
        bounds = [1, header["NAXIS1"], 1, header["NAXIS2"]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return WCS(header), bounds


def outer_corners(wcs, bounds):
    """RA/Dec of DATASEC outer corners, half a pixel beyond its centres."""
    x0, x1, y0, y1 = bounds
    return wcs.all_pix2world([x0 - 0.5, x1 + 0.5, x1 + 0.5, x0 - 0.5],
                             [y0 - 0.5, y0 - 0.5, y1 + 0.5, y1 + 0.5], 1)


def coverage_pixels(wcs, bounds) -> np.ndarray:
    """Sky pixels (NEST) whose centres lie inside the imaging area."""
    ra, dec = outer_corners(wcs, bounds)
    return hpg.query_polygon(NSIDE, ra, dec, nest=True)


def pointing_pixels(corners):
    """Rasterize the main-array hull in a tangent plane about its centre.

    @sc [decision:masking.pointing_square]
    CCDs 0–35 supply DATASEC outer corners; ears 36–39 are excluded and
    gaps are filled. Pixel centres inside the convex hull count once.
    """
    if not corners:
        return np.empty(0, dtype=np.int64)
    ra, dec = np.concatenate(corners, axis=1)
    r, d = np.deg2rad([ra, dec])
    xyz = np.array([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])
    centre = xyz.mean(axis=1)
    centre /= np.linalg.norm(centre)
    tangent = WCS(naxis=2)
    tangent.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    tangent.wcs.crval = [np.rad2deg(np.arctan2(centre[1], centre[0])) % 360,
                         np.rad2deg(np.arcsin(centre[2]))]
    xy = np.column_stack(tangent.all_world2pix(ra, dec, 0))
    hull = xy[ConvexHull(xy).vertices]
    ra, dec = tangent.all_pix2world(hull[:, 0], hull[:, 1], 0)
    return hpg.query_polygon(NSIDE, ra, dec, nest=True)


def flagged_counts(wcs, flags: np.ndarray, bounds):
    """Sky pixels holding flagged DATASEC pixel centres, and their counts."""
    x0, x1, y0, y1 = bounds
    rows, cols = np.nonzero(flags[y0 - 1:y1, x0 - 1:x1])
    ra, dec = wcs.all_pix2world(cols + x0, rows + y0, 1)
    return np.unique(hpg.angle_to_pixel(NSIDE, ra, dec, nest=True),
                     return_counts=True)


def rasterize(ccds, validated):
    """Build a fragment from an iterable of (CCD number, header, flag array).

    @sc [decision:masking.nexp_map_valid_psf_ccds]
    All CCD DATASEC footprints count as taken; only validated CCDs count as
    successful. A footprint counts an exposure once even at overlapping edges.
    @sc [decision:masking.defect_map_from_flags]
    Flags count only inside validated CCD DATASEC coverage, not overscan.
    """
    fragment = hsp.HealSparseMap.make_empty(
        NSIDE_COVERAGE, NSIDE, FRAGMENT_DTYPE)
    records, corners = {}, []
    present = set()
    for ccd, header, flags in ccds:
        present.add(ccd)
        wcs, bounds = ccd_wcs(header)
        if ccd < 36:
            corners.append(outer_corners(wcs, bounds))
        covered = coverage_pixels(wcs, bounds)
        values = fragment.get_values_pix(covered)
        values |= 2 | (4 if ccd in validated else 0)
        fragment[covered] = values
        record = {"covered": int(covered.size), "validated": ccd in validated}
        if ccd in validated:
            pixels, counts = flagged_counts(wcs, flags, bounds)
            keep = np.isin(pixels, covered, assume_unique=True)
            values = fragment.get_values_pix(pixels[keep])
            flagged = np.minimum((values >> FLAG_SHIFT) + counts[keep], 255)
            values = (values & COVER_MASK) | (flagged.astype(np.uint16) << FLAG_SHIFT)
            fragment[pixels[keep]] = values
            record["flagged_pixels"] = int(counts.sum())
        records[ccd] = record
    if set(validated) - present:
        raise ValueError(f"Validated CCDs absent from source: {set(validated) - present}")
    if not present:
        raise ValueError("No CCD images found")
    pixels = pointing_pixels(corners)
    values = fragment.get_values_pix(pixels)
    values |= 1
    fragment[pixels] = values
    return fragment, records


def source_ccds(args):
    """Yield CCD headers and flags from split files or paired delivered MEFs."""
    if args.exp_dir is not None:
        split = args.exp_dir / SPLIT_DIR
        images = sorted(split.glob(f"image-{args.exp}-*.fits"))
        for image in images:
            ccd = int(image.stem.rsplit("-", 1)[1])
            yield ccd, fits.getheader(image), fits.getdata(
                split / f"flag-{args.exp}-{ccd}.fits")
    else:
        with ExitStack() as stack:
            images = stack.enter_context(fits.open(args.image_mef))
            flags = stack.enter_context(fits.open(args.flag_mef))
            if len(images) != len(flags):
                raise ValueError("Image and flag MEFs have different HDU counts")
            for ccd in range(len(images) - 1):
                yield ccd, images[ccd + 1].header, flags[ccd + 1].data
                # Release decompressed flag arrays between CCDs.
                del flags[ccd + 1].data


def write_stable(tmp: Path, dest: Path) -> None:
    """Replace dest unless bytes match, preserving unchanged product mtimes."""
    if dest.exists() and filecmp.cmp(tmp, dest, shallow=False):
        tmp.unlink()
    else:
        tmp.replace(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--exp-dir", type=Path)
    source.add_argument("--image-mef", type=Path)
    parser.add_argument("--flag-mef", type=Path)
    parser.add_argument("--exp", required=True)
    parser.add_argument("--persist-manifest", required=True, type=Path)
    parser.add_argument("--fragment", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    if (args.image_mef is not None) != (args.flag_mef is not None):
        parser.error("--image-mef and --flag-mef must be supplied together")
    fragment, ccds = rasterize(source_ccds(args),
                               valid_ccds(args.persist_manifest, args.exp))
    args.fragment.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.fragment.with_name(f"{args.fragment.stem}.tmp.hsp")
    fragment.write(str(tmp), clobber=True)
    write_stable(tmp, args.fragment)
    body = {"stage": "exp_maps", "level": "exp", "unit": args.exp,
            "status": "complete", "fragment": str(args.fragment),
            "format": "packed-cover-nflagged-v1", "nside": NSIDE, "ccds": ccds}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.manifest.with_name(args.manifest.name + ".tmp")
    tmp.write_text(json.dumps(body, indent=1, sort_keys=True) + "\n")
    write_stable(tmp, args.manifest)
    print(f"[exp_maps] {args.exp}: {len(ccds)} CCD(s), "
          f"{sum(c['validated'] for c in ccds.values())} validated, "
          f"{fragment.valid_pixels.size} sky pixel(s)")


if __name__ == "__main__":
    main()
