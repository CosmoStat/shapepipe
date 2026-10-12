"""Rasterize ONE exposure's footprint and instrument defects into a HealSparse
fragment: the shell of the ``exp_maps`` rule.

The fragment is a uint8 map at the UNIONS mask ladder's resolution (nside
131072 over coverage 128, 1.6" pixels, ~74 MegaCam pixels each), so it aligns
pixel for pixel with the sky masks. A sky pixel holds

    0              not covered by a CCD of this exposure with a PSF model
    1 + n          covered, and ``n`` flagged CCD pixels have their centre in it

``merge_exposure_maps.py`` sums the fragments of a campaign into the
exposure-count map (``nexp``) and the flagged-pixel map (``nflagged``). The
CCDs of one exposure do not overlap, so the first sum counts exposures.

Flagged pixels are counted rather than OR-ed into a boolean mask because most
of them are thin: bad columns and cosmic rays a single CCD pixel wide. Any-touch
rasterization widens a one-pixel column to a 1.6" strip; on real exposures it
masks ~7x more sky than the flagged pixels cover. The count keeps the area
(``nflagged / 74`` is the fraction of one exposure lost there), and
``nflagged > 0`` is still the any-touch mask.

Inputs, from the exposure's scratch store and its exp_persist manifest:

* CCDs with a PSF model: the ``validation_psf-<exp>-<ccd>.fits`` members of
  ``exp_persist.json``. psfex_interp writes that file only when the fit
  succeeds, so these are the CCDs that can contribute a shape.
* WCS and imaging area: the header of ``image-<exp>-<ccd>.fits``. The split is
  2112 x 4644 pixels with a flagged overscan border (flag value 3); only
  ``DATASEC`` (2048 x 4612) sees the sky. The flag split carries no WCS.
* defects: ``flag-<exp>-<ccd>.fits``, any nonzero bit inside ``DATASEC``.

A sky pixel is covered when its centre lies inside the polygon through the
imaging area's outer corners.
"""

import argparse
import filecmp
import json
import warnings
from fnmatch import fnmatch
from pathlib import Path

import healsparse as hsp
import hpgeom as hpg
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import persist_exp

NSIDE = 131072
NSIDE_COVERAGE = 128
SPLIT_DIR = "output/run_sp_exp_Sp/split_exp_runner/output"
PSF_PATTERN = persist_exp.resolve(persist_exp.ALWAYS)


def valid_ccds(manifest: Path, exp: str) -> list:
    """CCD numbers with a PSF model, from the exp_persist manifest's members."""
    prefix = PSF_PATTERN.split("*")[0] + f"{exp}-"
    names = [f["name"] for f in json.loads(manifest.read_text())["files"]]
    return sorted(int(n[len(prefix):].removesuffix(".fits")) for n in names
                  if fnmatch(n, PSF_PATTERN) and n.startswith(prefix))


def ccd_wcs(image: Path):
    """The CCD's WCS and imaging area ``(x0, x1, y0, y1)``, 1-based inclusive,
    from ``DATASEC`` (or the whole array without one)."""
    header = fits.getheader(image)
    if "DATASEC" in header:
        x, y = header["DATASEC"].strip("[]").split(",")
        bounds = [int(v) for v in x.split(":") + y.split(":")]
    else:
        bounds = [1, header["NAXIS1"], 1, header["NAXIS2"]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return WCS(header), bounds


def coverage_pixels(wcs, bounds) -> np.ndarray:
    """Sky pixels (NEST) whose centres lie inside the imaging area."""
    x0, x1, y0, y1 = bounds
    ra, dec = wcs.all_pix2world([x0 - 0.5, x1 + 0.5, x1 + 0.5, x0 - 0.5],
                                [y0 - 0.5, y0 - 0.5, y1 + 0.5, y1 + 0.5], 1)
    return hpg.query_polygon(NSIDE, ra, dec, nest=True)


def flagged_counts(wcs, flags: np.ndarray, bounds):
    """Sky pixels (NEST) holding flagged imaging-area pixel centres, and how
    many each holds."""
    x0, x1, y0, y1 = bounds
    rows, cols = np.nonzero(flags[y0 - 1:y1, x0 - 1:x1])
    ra, dec = wcs.all_pix2world(cols + x0, rows + y0, 1)
    return np.unique(hpg.angle_to_pixel(NSIDE, ra, dec, nest=True),
                     return_counts=True)


def write_stable(tmp: Path, dest: Path) -> None:
    """Replace ``dest`` with ``tmp`` unless the bytes match, keeping the mtime
    of an unchanged product (mtime is a snakemake rerun trigger)."""
    if dest.exists() and filecmp.cmp(tmp, dest, shallow=False):
        tmp.unlink()
    else:
        tmp.replace(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp-dir", required=True, type=Path)
    parser.add_argument("--exp", required=True)
    parser.add_argument("--persist-manifest", required=True, type=Path)
    parser.add_argument("--fragment", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    split = args.exp_dir / SPLIT_DIR
    fragment = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint8)
    ccds = {}
    for ccd in valid_ccds(args.persist_manifest, args.exp):
        wcs, bounds = ccd_wcs(split / f"image-{args.exp}-{ccd}.fits")
        covered = coverage_pixels(wcs, bounds)
        pixels, counts = flagged_counts(
            wcs, fits.getdata(split / f"flag-{args.exp}-{ccd}.fits"), bounds)
        keep = np.isin(pixels, covered, assume_unique=True)
        fragment[covered] = np.ones(covered.size, np.uint8)
        fragment[pixels[keep]] = (1 + np.minimum(counts[keep], 254)).astype(np.uint8)
        ccds[ccd] = {"covered": int(covered.size),
                     "flagged_pixels": int(counts.sum())}

    args.fragment.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.fragment.with_name(f"{args.fragment.stem}.tmp.hsp")
    fragment.write(str(tmp), clobber=True)
    write_stable(tmp, args.fragment)

    body = {"stage": "exp_maps", "level": "exp", "unit": args.exp,
            "status": "complete", "fragment": str(args.fragment),
            "nside": NSIDE, "ccds": ccds}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.manifest.with_name(args.manifest.name + ".tmp")
    tmp.write_text(json.dumps(body, indent=1, sort_keys=True) + "\n")
    write_stable(tmp, args.manifest)
    covered = sum(c["covered"] for c in ccds.values())
    flagged = sum(c["flagged_pixels"] for c in ccds.values())
    print(f"[exp_maps] {args.exp}: {len(ccds)} CCD(s) with a PSF model, "
          f"{covered} sky pixel(s), {flagged} flagged CCD pixel(s)")


if __name__ == "__main__":
    main()
