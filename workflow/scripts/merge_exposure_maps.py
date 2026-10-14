"""Sum per-exposure fragments into four campaign HealSparse maps.

* ``npoint_<run>.hsp`` (uint8): filled main-array pointing hulls.
* ``ntaken_<run>.hsp`` (uint8): CCD DATASEC coverage, including ears.

* ``nexp_<run>.hsp`` (uint8): per sky pixel, the number of exposures whose CCD
  with a valid PSF model covers it.
* ``nflagged_<run>.hsp`` (uint16): per sky pixel, the flagged CCD pixels of
  those exposures whose centres fall in it. A sky pixel holds ~74 MegaCam
  pixels, so ``nflagged / 74`` is the number of exposures' worth of area lost
  there, and ``nflagged > 0`` the any-touch defect mask.

All are 0 where nothing is counted. The campaign's exposures come from the
tile list and the index (``build_index``), as for the star-catalogue merge.
Missing fragments are omitted and reported; incompatible formats fail with a
backfill command. The maps are rebuilt from every fragment on each run.
"""

import argparse
from pathlib import Path

import build_index
import healsparse as hsp
import numpy as np
from exp_maps import FLAG_SHIFT, FRAGMENT_DTYPE, NSIDE, NSIDE_COVERAGE


def fragment_path(products_dir: Path, exp: str) -> Path:
    """Return the fragment path matching Snakefile's prod_exp_maps."""
    return products_dir / "exp" / exp[:2] / exp / "maps" / f"maps-{exp}.hsp"


def add(acc, pixels, values):
    """Add ``values`` to ``acc`` at ``pixels``."""
    acc[pixels] = acc.get_values_pix(pixels) + values.astype(acc.dtype)


def main() -> None:
    """Sum compatible fragments into the four campaign maps."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products-dir", required=True, type=Path)
    parser.add_argument("--tile-list", required=True, type=Path)
    parser.add_argument("--index-db", required=True, type=Path)
    parser.add_argument("--npoint", required=True, type=Path)
    parser.add_argument("--ntaken", required=True, type=Path)
    parser.add_argument("--nexp", required=True, type=Path)
    parser.add_argument("--nflagged", required=True, type=Path)
    args = parser.parse_args()

    exposures = build_index.campaign_exposures(args.tile_list, args.index_db)
    npoint = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint8)
    ntaken = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint8)
    nexp = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint8)
    nflagged = hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, np.uint16)
    missing = []
    for i, exp in enumerate(exposures):
        path = fragment_path(args.products_dir, exp)
        if not path.exists():
            missing.append(exp)
            continue
        fragment = hsp.HealSparseMap.read(str(path))
        if fragment.dtype != FRAGMENT_DTYPE:
            persist = path.parent.parent / "manifests" / "exp_persist.json"
            manifest = path.parent.parent / "manifests" / "exp_maps.json"
            raw = Path("/project/def-mjhudson/unions-wl/exposures")
            raise SystemExit(
                f"Exposure {exp}: incompatible fragment {path}. "
                "Backfill with:\n"
                f"python workflow/scripts/exp_maps.py --exp {exp} "
                f"--image-mef {raw / (exp + 'p.fits.fz')} "
                f"--flag-mef {raw / (exp + 'p.flag.fits.fz')} "
                f"--persist-manifest {persist} --fragment {path} "
                f"--manifest {manifest}")
        pixels = fragment.valid_pixels
        values = fragment.get_values_pix(pixels)
        for acc, bit in ((npoint, 1), (ntaken, 2), (nexp, 4)):
            selected = pixels[(values & bit) != 0]
            add(acc, selected, np.ones(selected.size, dtype=np.uint8))
        counts = values >> FLAG_SHIFT
        flagged = counts > 0
        add(nflagged, pixels[flagged], counts[flagged])
        if i % 500 == 0:
            print(f"[exposure_maps] {i + 1}/{len(exposures)} {exp}",
                  flush=True)

    used = len(exposures) - len(missing)
    if not used:
        raise SystemExit("merge_exposure_maps: no exposure of the campaign "
                         "has a fragment")
    for acc, out in ((npoint, args.npoint), (ntaken, args.ntaken),
                     (nexp, args.nexp), (nflagged, args.nflagged)):
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_name(f"{out.stem}.tmp.hsp")
        acc.write(str(tmp), clobber=True)
        tmp.replace(out)
    print(f"[exposure_maps] {used} exposure(s) -> {args.npoint.name}, "
          f"{args.ntaken.name}, {args.nexp.name}, "
          f"{args.nflagged.name}")
    if missing:
        print(f"[exposure_maps] WARNING: {len(missing)} exposure(s) in the "
              f"campaign have no fragment (reclaimed before exp_maps ran) and "
              f"are not counted, e.g. {missing[:5]}")


if __name__ == "__main__":
    main()
