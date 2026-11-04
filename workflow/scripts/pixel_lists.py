"""Build order-5 tile and exposure candidate lists from the run index."""

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


def tile_exposure_index(tiles, db):
    """Require prepared, nonempty and consistent incidence for every tile."""
    requested = set(tiles)
    with sqlite3.connect(f'file:{Path(db).resolve()}?mode=ro', uri=True) as conn:
        ready = dict(conn.execute('SELECT tile_id, n_exp FROM tiles'))
        edges = {t: set() for t in requested}
        for tile, exp in conn.execute('SELECT tile_id, exp_id FROM tile_exposures'):
            if tile in edges:
                edges[tile].add(exp)
    missing = sorted(t for t in requested if t not in ready or not edges[t]
                     or ready[t] != len(edges[t]))
    if missing:
        raise ValueError(f'missing or incomplete tile exposure associations: {missing}')
    return edges


def exposure_pairs(tile_pixel_pairs, db):
    """Union each listed tile's exposures into its pixels.

    @sc [decision:survey_publication.pixel_completeness,label:invariant] exposure-list-covers-stars
    Every star row must land in a pixel listing its exposure: a star belongs
    to a tile using that exposure, and that tile's cap covers the row.
    """
    pairs = list(tile_pixel_pairs)
    edges = tile_exposure_index((tile for _, tile in pairs), db)
    return {(pixel, exp) for pixel, tile in pairs for exp in edges[tile]}


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
    parser.add_argument('--pixel-tiles', required=True, type=Path)
    parser.add_argument('--pixel-exposures', required=True, type=Path)
    args = parser.parse_args()
    tiles = read_tiles(args.tiles)
    # Compute both before writing either, so incomplete incidence cannot publish lists.
    tp = list(tile_pairs(tiles))
    ep = list(exposure_pairs(tp, args.index))
    write_pairs(args.pixel_tiles, 'tile_id', tp)
    write_pairs(args.pixel_exposures, 'exp_id', ep)


if __name__ == '__main__':
    main()
