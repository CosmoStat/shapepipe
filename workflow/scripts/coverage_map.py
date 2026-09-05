#!/usr/bin/env python3
"""Build the campaign's HealSparse coverage (nexp) map from the exposure footprints.

Run as the shell of the in-DAG ``coverage_map`` rule, never by hand.

WHAT THE MAP IS. Per sky pixel, the number of exposures with a VALID PSF MODEL
covering it. The CCDs of one exposure do not overlap, so stamping value 1 per CCD
polygon and accumulating counts exposures. sp_validation applies it as a
structural mask (``notebooks/demo_apply_hsp_masks.py``).

CAMPAIGN-CUMULATIVE, AND THAT IS THE POINT. This script GLOBS every
``<products_dir>/exp/*/*/manifests/exp_footprint.json`` — not just the ones the
rule declared as inputs, and INCLUDING exposures whose scratch stores have been
reclaimed. The declared inputs are the in-scope, non-tombstoned footprints, which
is what buys ordering and rerun semantics without dragging out-of-scope tiles
into the DAG; the records themselves live on the persistent root and stay valid
sky forever. So appending tiles GROWS the map rather than replacing it, which is
what a survey coverage mask should do. The rule's comment says the same thing
where a reader of the DAG will meet it.

THE STAMPING IS NOT HERE. It is ``shapepipe.utilities.coverage_map_builder.
build_map``, shared with the ``build_coverage_map`` CLI, so the RA-seam guard,
the pole guard and the polygon accumulation have exactly one implementation. This
script is the JSON-to-arrays half, and the records are raw sky: unwrapping across
RA=0 happens once, inside build_map, at the moment a polygon is stamped.

NSIDE IS NOT DEFAULTED HERE. Both values are required arguments, carried from
`coverage:` in config.yaml, because the production pair (128 / 131072) is NOT the
CoverageMapBuilder default (32 / 2048) and the difference is invisible in the
output: nside=131072 is ~0.1"/pixel, chosen to match the UNIONS bit-mask
resolution so coverage and mask align pixel-wise. A map built at the class
default would look entirely reasonable and would not align, and the consumer
would not notice.
"""

import argparse
import filecmp
import json
import sys
from pathlib import Path

import numpy as np

from shapepipe.utilities.coverage_map_builder import build_map

# Where exp_footprint writes, relative to the products root. The shard and the
# exposure id are both globbed: this map is the whole campaign's, so it is
# deliberately not built from a list of units.
FOOTPRINT_GLOB = "exp/*/*/manifests/exp_footprint.json"


def read_footprints(products_dir):
    """Every exposure footprint on the persistent root, as flat arrays.

    Returns ``(ccd_ids, ra, dec, units)``: length-M id and ``(M, 4)`` corner
    arrays over all CCDs of all exposures, plus the sorted exposure ids that
    contributed. Records are read in sorted path order so the polygon order —
    and hence the map — does not depend on readdir order.
    """
    paths = sorted(Path(products_dir).glob(FOOTPRINT_GLOB))
    if not paths:
        sys.exit(f"coverage_map: no {FOOTPRINT_GLOB} under {products_dir}; "
                 f"there is nothing to build a map from")

    ccd_ids, ra, dec, units = [], [], [], []
    for path in paths:
        body = json.loads(path.read_text())
        units.append(body["unit"])
        for ccd in body["ccds"]:
            ccd_ids.append(ccd["id"])
            ra.append(ccd["ra"])
            dec.append(ccd["dec"])
    return (np.array(ccd_ids, dtype=str),
            np.array(ra, dtype=float).reshape(-1, 4),
            np.array(dec, dtype=float).reshape(-1, 4),
            sorted(units))


def write_stable(path, body):
    """Write JSON tmp-then-``cmp``-then-``mv``; an unchanged body keeps its mtime."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
        if path.exists() and filecmp.cmp(tmp, path, shallow=False):
            tmp.unlink()
        else:
            tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--products-dir", required=True, type=Path,
                   help="the persistent root; every exposure footprint under "
                        "it goes into the map")
    p.add_argument("--out", required=True, type=Path,
                   help="the HealSparse map, <products_dir>/coverage/coverage.hsp")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--nside-coverage", required=True, type=int)
    p.add_argument("--nside", required=True, type=int)
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    ccd_ids, ra, dec, units = read_footprints(args.products_dir)
    print(f"[coverage_map] {len(units)} exposure(s), {len(ccd_ids)} CCD "
          f"footprint(s) from {args.products_dir}")

    hsp_map = build_map(ccd_ids, ra, dec, args.nside_coverage, args.nside,
                        verbose=args.verbose)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    hsp_map.write(str(args.out), clobber=True)

    # The manifest is written AFTER the map, and it is the rule's record of
    # which exposures the map contains — the question a mask's consumer asks
    # months later, and one nothing else on disk can answer once the campaign
    # has grown past it.
    write_stable(args.manifest, {
        "stage": "coverage_map", "level": "campaign", "status": "complete",
        "map": str(args.out),
        "nside_coverage": args.nside_coverage,
        "nside": args.nside,
        "n_exposures": len(units),
        "n_ccds": len(ccd_ids),
        "exposures": units,
    })
    print(f"[coverage_map] -> {args.out}")


if __name__ == "__main__":
    main()
