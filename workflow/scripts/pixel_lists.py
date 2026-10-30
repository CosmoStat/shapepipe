"""Build order-5 tile and exposure candidate lists from enclosing caps.

Exposure caps are supplied as CSV columns exp_id,ra,dec,radius_deg, in degrees.
The run index supplies the tile-to-exposure association; every associated
exposure must have exactly one cap. Caps must enclose the complete pointing,
including ear CCDs and stars. The caller is responsible for deriving them from
the exposure's WCS footprints, not from the associated tile centres.
"""

import argparse
import csv
import math
import sqlite3
from pathlib import Path

import healpy as hp
import numpy as np

from tile_universe import read_tiles, validate_tiles

# @sc [decision:survey_publication.pixel_completeness]
TILE_CAP_DEG = 0.4


def tile_centre(tile):
    """Return the unwrapped CFIS tile centre as longitude and latitude."""
    x, y = map(int, tile.split('.'))
    dec = y / 2 - 90
    return x / (2 * math.cos(math.radians(dec))), dec


def cap_pixels(ra, dec, radius):
    """Conservatively find every order-5 pixel touched by a cap."""
    if not all(math.isfinite(v) for v in (ra, dec, radius)) or not -90 <= dec <= 90 or not 0 < radius <= 180:
        raise ValueError('invalid spherical cap')
    return hp.query_disc(32, hp.ang2vec(ra % 360, dec, lonlat=True),
                         np.deg2rad(radius), inclusive=True, fact=16, nest=True)


def tile_pairs(tiles):
    """List conservative pixel/tile pairs.

    @sc [decision:survey_publication.pixel_completeness,label:invariant] tile-cap-encloses-rows
    Every tile row must lie within the tile's 0.4 degree cap, so the
    conservative pixel query includes every pixel to which the tile contributes.
    """
    for tile in validate_tiles(tiles):
        ra, dec = tile_centre(tile)
        for pixel in cap_pixels(ra, dec, TILE_CAP_DEG):
            yield int(pixel), tile


def exposure_pairs(tiles, db, caps):
    """List pixel/exposure pairs for all exposures associated with tiles."""
    tile_set = set(tiles)
    with sqlite3.connect(f'file:{Path(db).resolve()}?mode=ro', uri=True) as conn:
        exposures = {exp for tile, exp in conn.execute('SELECT tile_id, exp_id FROM tile_exposures') if tile in tile_set}
    values = {}
    with open(caps, newline='') as f:
        for row in csv.DictReader(f):
            exp = row['exp_id']
            if exp in values:
                raise ValueError(f'duplicate exposure cap: {exp}')
            values[exp] = tuple(float(row[k]) for k in ('ra', 'dec', 'radius_deg'))
    missing = exposures - values.keys()
    if missing:
        raise ValueError(f'missing exposure caps: {sorted(missing)[:20]}')
    for exp in sorted(exposures):
        for pixel in cap_pixels(*values[exp]):
            yield int(pixel), exp


def write_pairs(path, column, pairs):
    """Write sorted, unique pairs with the shared publication schema."""
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(('pixel', column))
        writer.writerows(sorted(set(pairs)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tiles', required=True, type=Path)
    parser.add_argument('--index', required=True, type=Path)
    parser.add_argument('--exposure-caps', required=True, type=Path)
    parser.add_argument('--pixel-tiles', required=True, type=Path)
    parser.add_argument('--pixel-exposures', required=True, type=Path)
    args = parser.parse_args()
    tiles = read_tiles(args.tiles)
    # Compute both before writing either, so missing caps do not leave a partial pair.
    tp = list(tile_pairs(tiles))
    ep = list(exposure_pairs(tiles, args.index, args.exposure_caps))
    write_pairs(args.pixel_tiles, 'tile_id', tp)
    write_pairs(args.pixel_exposures, 'exp_id', ep)


if __name__ == '__main__':
    main()
