#!/usr/bin/env python3
"""The campaign's GSC 2.3 star catalogue, as a HEALPix-chunked sky store.

Masking needs, for every exposure, the bright stars over its focal plane. The
sky does not change between exposures, so the network cost of that is a property
of the campaign's SKY AREA, not of its exposure count: exposures overlap each
other ~7-10 deep, and a tile's exposures all look at the same square degree.

So the store is chunked by sky, not by exposure. One GSC 2.3 cone query per
HEALPix pixel of NSIDE=32, written run-independently under the ``star_cats``
config root and never fetched twice. A campaign that grows past the fetched
footprint queries only the chunks its new tiles add; one that grows within it
queries nothing.

Two numbers set the scale. A full-UNIONS footprint is ~1.5k chunks against ~25k
exposures, so the QUERY COUNT drops ~16x. The queried AREA drops ~4x: the old
design covered the footprint ~8 times over (that is just the exposure overlap
depth), the new one ~2 times, the 2x being the price of bounding a HEALPix
quadrilateral by the cone Vizier speaks (see ``pixel_cone``) — 5.6-9 deg^2 for a
3.36 deg^2 pixel, ~40-60k rows and ~3-4 MB per chunk.

Two subcommands, one module, deliberately: ``fetch`` and ``cut`` must agree
EXACTLY on which pixel holds which star, and a shared NSIDE constant in one file
is the only version of that agreement which cannot drift.

    fetch --tile-list ... --store ... --manifest ...
        The campaign side. Turns the tile list into the set of pixels its
        exposures can possibly need, fetches the missing ones, writes a manifest.

    cut --images ... --store ... --out ...
        The per-exposure side, purely local: read the focal-plane footprint from
        the exposure's image headers, load the chunks covering it, deduplicate,
        and cut to the focal-plane disc. Reproduces byte-for-byte the same sky
        selection the old one-query-per-exposure cone did.

Geometry, and why the fetch pad is what it is. Chunk-need is computed from the
TILE list rather than from exposure pointings, because tile IDs are the one thing
known before any download: a pointing center means reading a FITS header of an
image get_images has not fetched yet, and the DAG needs the chunk set at parse
time. Tiles sit on a fixed 0.5 deg grid (``cfis.get_tile_coord_from_nixy``), so
each tile is a disc of half-diagonal 0.354 deg; find_exposures gives a tile every
exposure whose footprint covers it, and the MegaCam focal plane is a disc of
radius 0.73 deg (measured on the cached catalogues). An exposure center is
therefore at most 0.354 + 0.73 deg from the tile center, and its stars 0.73 deg
beyond that: 1.81 deg, padded to ``PAD_DEG`` = 2.0. The pad is a perimeter cost
— negligible for a contiguous campaign, and paid once.

The pad is a bound, not a promise: ``cut`` verifies that every chunk covering the
exposure it was handed is on disk, and fails loudly if one is not. A missing
chunk means the geometry above is wrong, and that must not degrade quietly into
an under-masked exposure.
"""

import argparse
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS

# GSC 2.3. The same catalogue the mask module's own CDS path uses
# (mask.py: _CDS_cat_ID), so the store is a drop-in for it.
CAT_ID = "I/305/out"

# NSIDE=32 -> 3.36 deg^2 per pixel, 12288 pixels over the sky. Chosen so one
# chunk is a couple of MegaCam focal planes: small enough that a Vizier query
# stays within a small multiple of the per-exposure queries this replaces, large
# enough that a full-UNIONS footprint is ~1.5k chunks rather than ~25k. The
# cone-vs-quadrilateral overhead is scale-free, so NSIDE trades query count
# against query size and nothing else. NESTED, so a chunk id is a hierarchical
# sky address and a future NSIDE change is a subdivision.
NSIDE = 32
NEST = True

# Angular padding on the disc used to select chunks (see the module docstring).
PAD_DEG = 2.0

# The MegaCam focal-plane disc, and the margin added to a pixel's own bounding
# cone. Both in degrees.
MARGIN_DEG = 0.02

# GSC 2.3's object id: the deduplication key where chunk cones overlap.
ID_COL = "GSC2.3"


# --- the chunk store --------------------------------------------------------


def store_dir(store: Path) -> Path:
    """Chunks live under the catalogue and resolution that produced them, so a
    later NSIDE or catalogue change is a new directory beside the old one rather
    than a silent reinterpretation of files already on disk."""
    return Path(store) / CAT_ID.replace("/", "_") / f"nside{NSIDE}"


def chunk_path(store: Path, ipix: int) -> Path:
    return store_dir(store) / f"star_chunk-{ipix:06d}.fits"


def chunks_for_disc(ra_deg: float, dec_deg: float, radius_deg: float) -> list[int]:
    """Every pixel that touches the disc, as sorted ids.

    ``inclusive=True`` makes this a conservative superset — the guarantee ``cut``
    relies on is that no star inside the disc lives in a pixel this omits.
    """
    vec = hp.ang2vec(ra_deg, dec_deg, lonlat=True)
    return sorted(int(i) for i in hp.query_disc(
        NSIDE, vec, np.radians(radius_deg), inclusive=True, fact=4, nest=NEST))


def pixel_cone(ipix: int) -> tuple[float, float, float]:
    """(ra, dec, radius_arcmin) of a cone that CONTAINS pixel ``ipix``.

    Vizier speaks cones, HEALPix speaks quadrilaterals, so the query is the
    pixel's bounding cone: its center, and the largest center-to-boundary
    distance plus a margin. The cone spills over the pixel edges, which costs a
    little duplication between neighbours and buys the containment ``cut``
    depends on. The duplicates are removed on read, by ``ID_COL``.
    """
    ra_c, dec_c = hp.pix2ang(NSIDE, ipix, nest=NEST, lonlat=True)
    ra_b, dec_b = hp.vec2ang(hp.boundaries(NSIDE, ipix, step=8, nest=NEST).T,
                             lonlat=True)
    center = SkyCoord(ra_c * u.deg, dec_c * u.deg)
    radius = center.separation(SkyCoord(ra_b * u.deg, dec_b * u.deg)).deg.max()
    return float(ra_c), float(dec_c), float((radius + MARGIN_DEG) * 60.0)


def write_atomic(table: Table, path: Path) -> None:
    """Publish ``table`` at ``path`` all-or-nothing.

    The store's only cache test is ``Path.exists``, so a write killed part-way
    (timeout, OOM, node failure) would otherwise leave a truncated FITS that
    every later run trusts forever. The temp keeps the ``.fits`` suffix because
    astropy picks its writer from the extension, and is dot-prefixed and
    PID-tagged so it stays out of the ``star_chunk-*`` globs and two concurrent
    writers cannot collide.
    """
    tmp = path.parent / f".tmp-{os.getpid()}-{path.name}"
    try:
        table.write(tmp, overwrite=True)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def read_chunks(store: Path, ipixels: list[int]) -> Table:
    """Load and deduplicate the given chunks.

    A missing chunk is fatal (see the module docstring): it means the fetch
    footprint did not cover this exposure, and an under-masked exposure is worse
    than a failed job.
    """
    missing = [i for i in ipixels if not chunk_path(store, i).exists()]
    if missing:
        raise SystemExit(
            f"star chunk(s) {missing} not in {store_dir(store)}. The campaign's "
            f"star_catalogue fetch did not cover this exposure — re-run it "
            f"(and check that its tile list contains this exposure's tiles).")

    table = vstack([Table.read(chunk_path(store, i)) for i in ipixels],
                   metadata_conflicts="silent")
    _, keep = np.unique(np.asarray(table[ID_COL]), return_index=True)
    return table[np.sort(keep)]


# --- exposure footprint -----------------------------------------------------


def _wcs(header) -> WCS:
    """Build the WCS by hand, from the linear terms only.

    Same construction as ``scripts/python/create_star_cat.py``: it sidesteps
    distortion-convention incompatibilities between headers and astropy, and a
    focal-plane footprint needs nothing finer.
    """
    w = WCS(naxis=2)
    w.wcs.ctype = [header["CTYPE1"], header["CTYPE2"]]
    try:
        w.wcs.cunit = [header["CUNIT1"], header["CUNIT2"]]
    except KeyError:
        w.wcs.cunit = ["deg", "deg"]
    w.wcs.crpix = [header["CRPIX1"], header["CRPIX2"]]
    w.wcs.crval = [header["CRVAL1"], header["CRVAL2"]]
    w.wcs.cd = [[header["CD1_1"], header["CD1_2"]],
                [header["CD2_1"], header["CD2_2"]]]
    return w


def focal_plane_disc(image: Path, n_ccd: int = 40) -> tuple[float, float, float]:
    """(ra, dec, radius_deg) of the disc covering all CCDs of one exposure."""
    # ONE open for all 40 CCDs. `fits.getheader(image, ext)` opens the file,
    # walks the HDU list to `ext` and closes again, so the loop cost 40 opens and
    # O(n^2) header seeks over a compressed multi-extension exposure.
    centers, radii = [], []
    with fits.open(image) as hdul:
        for ext in range(1, n_ccd + 1):
            h = hdul[ext].header
            w = _wcs(h)
            (ra_c, dec_c), (ra_0, dec_0) = w.all_pix2world(
                [[h["NAXIS1"] / 2.0, h["NAXIS2"] / 2.0], [0, 0]], 1)
            centers.append(SkyCoord(ra_c * u.deg, dec_c * u.deg))
            radii.append(centers[-1].separation(
                SkyCoord(ra_0 * u.deg, dec_0 * u.deg)).deg)

    ras = np.array([c.ra.deg for c in centers])
    decs = np.array([c.dec.deg for c in centers])
    center = SkyCoord(ras.mean() * u.deg, decs.mean() * u.deg)
    seps = center.separation(SkyCoord(ras * u.deg, decs * u.deg)).deg
    return (float(center.ra.deg), float(center.dec.deg),
            float(np.max(seps + np.array(radii))))


def exposure_image(images_dir: Path) -> Path:
    """The one multi-extension exposure image in a get_images output dir.

    That dir is a symlink farm holding ``image-<exp>.fitsfz`` plus its weight and
    flag; only the image carries the 40 CCD WCSs.
    """
    found = sorted(p for p in Path(images_dir).iterdir() if "image" in p.name)
    if not found:
        raise SystemExit(f"no image file in {images_dir}")
    return found[0]


# --- fetch ------------------------------------------------------------------


def campaign_chunks(tile_ids: list[str]) -> list[int]:
    """Every chunk the campaign's exposures can need, from the tile list alone."""
    from shapepipe.utilities.cfis import get_tile_coord_from_nixy

    needed: set[int] = set()
    for tile_id in tile_ids:
        nix, niy = tile_id.split(".")
        ra, dec = get_tile_coord_from_nixy(nix, niy)
        needed.update(chunks_for_disc(ra.degree, dec.degree, PAD_DEG))
    return sorted(needed)


def fetch(args: argparse.Namespace) -> None:
    from shapepipe.utilities.vizier import query_vizier

    tile_ids = [ln.strip() for ln in Path(args.tile_list).read_text().splitlines()
                if ln.strip()]
    needed = campaign_chunks(tile_ids)
    out_dir = store_dir(args.store)
    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [i for i in needed if not chunk_path(args.store, i).exists()]
    print(f"star chunks: {len(needed)} needed for {len(tile_ids)} tiles, "
          f"{len(todo)} to fetch -> {out_dir}", file=sys.stderr)

    def one(ipix: int) -> int:
        ra, dec, radius_arcmin = pixel_cone(ipix)
        table = query_vizier(ra, dec, radius_arcmin, CAT_ID)
        write_atomic(table, chunk_path(args.store, ipix))
        print(f"chunk {ipix}: {len(table)} rows "
              f"(ra={ra:.4f} dec={dec:.4f} r={radius_arcmin:.1f}')",
              file=sys.stderr)
        return len(table)

    # A handful of concurrent queries, never one per exposure: the same modest
    # concurrency the per-exposure rule reached through --local-cores, now an
    # explicit number instead of an accident of the head node's CPU count.
    if todo:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(one, todo))

    write_manifest(args, tile_ids, needed, len(todo))


def write_manifest(args, tile_ids, needed, n_fetched) -> None:
    """The rule's declared output.

    Not a ``completeness.py`` verdict: this rule runs no ``shapepipe_run`` and
    has no per-runner count floors, so there is nothing to compose and no
    separate ``log:`` — under ``set -euo pipefail`` the job either completes or
    aborts at the failing query, and snakemake's captured stderr is the evidence.
    The manifest keeps the workflow's "one rule, one manifest" currency: written
    last, and only when the content changed, so an unchanged campaign leaves the
    mtime where it was rather than churning the `mtime` rerun-trigger.
    """
    body = json.dumps({
        "stage": "star_catalogue",
        "level": "campaign",
        "status": "complete",
        "catalogue": CAT_ID,
        "nside": NSIDE,
        "nest": NEST,
        "pad_deg": PAD_DEG,
        "store": str(store_dir(args.store)),
        "n_tiles": len(tile_ids),
        "n_chunks": len(needed),
        "n_fetched": n_fetched,
        "chunks": needed,
    }, indent=2, sort_keys=True)
    path = Path(args.manifest)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != body:
        path.write_text(body)


# --- cut --------------------------------------------------------------------


def cut(args: argparse.Namespace) -> None:
    image = exposure_image(args.images)
    ra, dec, radius = focal_plane_disc(image)
    ipixels = chunks_for_disc(ra, dec, radius)
    table = read_chunks(args.store, ipixels)

    center = SkyCoord(ra * u.deg, dec * u.deg)
    stars = SkyCoord(np.asarray(table["RAJ2000"]) * u.deg,
                     np.asarray(table["DEJ2000"]) * u.deg)
    inside = table[center.separation(stars).deg <= radius]

    print(f"{image.name}: ra={ra:.4f} dec={dec:.4f} r={radius:.4f} deg, "
          f"{len(ipixels)} chunks -> {len(inside)} stars", file=sys.stderr)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_atomic(inside, out)


# --- CLI --------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="fetch the campaign footprint's chunks")
    f.add_argument("--tile-list", required=True, type=Path)
    f.add_argument("--store", required=True, type=Path)
    f.add_argument("--manifest", required=True, type=Path)
    f.add_argument("--workers", type=int, default=4)
    f.set_defaults(func=fetch)

    c = sub.add_parser("cut", help="cut one exposure's catalogue from the store")
    c.add_argument("--images", required=True, type=Path,
                   help="a get_images output dir holding image-<exp>.fitsfz")
    c.add_argument("--store", required=True, type=Path)
    c.add_argument("--out", required=True, type=Path)
    c.set_defaults(func=cut)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
