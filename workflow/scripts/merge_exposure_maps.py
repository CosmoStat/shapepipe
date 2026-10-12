"""Sum a campaign's per-exposure fragments into its two HealSparse maps: the
shell of the ``exposure_maps`` rule.

* ``nexp_<run>.hsp`` (uint8): per sky pixel, the number of exposures whose CCD
  with a valid PSF model covers it, the count behind sp_validation's
  ``npoint`` cut.
* ``nflagged_<run>.hsp`` (uint16): per sky pixel, the flagged CCD pixels of
  those exposures whose centres fall in it. A sky pixel holds ~74 MegaCam
  pixels, so ``nflagged / 74`` is the number of exposures' worth of area lost
  there, and ``nflagged > 0`` the any-touch defect mask.

Both are 0 where nothing is counted. Membership comes from
``build_index.campaign_exposures``. Missing fragments are omitted and reported;
this script does not determine why they are absent. The maps are rebuilt from
every available campaign fragment on each run. See ``exp_maps.py`` for fragment
encoding, coverage selection and count saturation.
"""

import argparse
from pathlib import Path

import healsparse as hsp
import numpy as np

import build_index
from exp_maps import NSIDE, NSIDE_COVERAGE


def fragment_path(products_dir: Path, exp: str) -> Path:
    """Where ``exp_maps`` writes an exposure's fragment (Snakefile:
    ``prod_exp_maps``)."""
    return products_dir / "exp" / exp[:2] / exp / "maps" / f"maps-{exp}.hsp"


def add(acc, pixels, values):
    """Add ``values`` to ``acc`` at ``pixels``."""
    acc[pixels] = acc.get_values_pix(pixels) + values.astype(acc.dtype)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products-dir", required=True, type=Path)
    parser.add_argument("--tile-list", required=True, type=Path)
    parser.add_argument("--index-db", required=True, type=Path)
    parser.add_argument("--nexp", required=True, type=Path)
    parser.add_argument("--nflagged", required=True, type=Path)
    args = parser.parse_args()

    exposures = build_index.campaign_exposures(args.tile_list, args.index_db)
    nexp = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint8)
    nflagged = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint16)
    missing = []
    for i, exp in enumerate(exposures):
        path = fragment_path(args.products_dir, exp)
        if not path.exists():
            missing.append(exp)
            continue
        fragment = hsp.HealSparseMap.read(str(path))
        pixels = fragment.valid_pixels
        values = fragment.get_values_pix(pixels)
        add(nexp, pixels, np.ones_like(values))
        flagged = values > 1
        add(nflagged, pixels[flagged], values[flagged] - 1)
        if i % 500 == 0:
            print(f"[exposure_maps] {i + 1}/{len(exposures)} {exp}", flush=True)

    used = len(exposures) - len(missing)
    if not used:
        raise SystemExit("merge_exposure_maps: no exposure of the campaign has "
                         "a fragment")
    for acc, out in ((nexp, args.nexp), (nflagged, args.nflagged)):
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_name(f"{out.stem}.tmp.hsp")
        acc.write(str(tmp), clobber=True)
        tmp.replace(out)
    print(f"[exposure_maps] {used} exposure(s) -> {args.nexp.name}, "
          f"{args.nflagged.name}")
    if missing:
        print(f"[exposure_maps] WARNING: {len(missing)} exposure(s) in the "
              f"campaign have no fragment (reclaimed before exp_maps ran) and "
              f"are not counted, e.g. {missing[:5]}")


if __name__ == "__main__":
    main()
