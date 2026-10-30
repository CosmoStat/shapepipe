"""Create, validate and extend order-2 cell claims.

A claim is a git commit to the shared claims file. Push it before using the
cell; a rejected push is the lock. Pull and recompute the proposal after a
rejection. This script does not run git.
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import healpy as hp
import numpy as np
import yaml
from scipy.spatial import cKDTree

from pixel_lists import tile_centre, tile_exposure_index
from tile_universe import read_tiles, validate_tiles


class UniqueLoader(yaml.SafeLoader):
    """Safe YAML loader that refuses duplicate mapping keys."""


def _mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f'duplicate YAML key: {key}')
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_claims(path):
    """Validate order, cells and owners before planning.

    @sc [decision:survey_publication.seam_realization,label:invariant] one-owner-per-pixel
    An order-2 cell has at most one cluster owner; each order-5 pixel inherits
    its parent's owner. Duplicate YAML keys must not silently overwrite claims.
    """
    data = yaml.load(Path(path).read_text(), Loader=UniqueLoader)
    if not isinstance(data, dict) or data.get('order') != 2 or not isinstance(data.get('cells'), dict):
        raise ValueError('claims must contain order: 2 and a cells mapping')
    cells = data['cells']
    if any(type(k) is not int or not 0 <= k < 192 or v not in ('nibi', 'fir') for k, v in cells.items()):
        raise ValueError('invalid cell or cluster owner')
    return cells


def initial(tiles, fraction=0.6):
    """Reproduce the block-first 40:60 split in two-way/analyze.py."""
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    counts, block = Counter(), Counter()
    for tile in validate_tiles(tiles):
        ra, dec = tile_centre(tile)
        cell = int(hp.ang2pix(4, ra % 360, dec, nest=True, lonlat=True))
        counts[cell] += 1
        if (ra % 360 >= 290 or ra % 360 <= 70) and dec < 42:
            block[cell] += 1
    if not counts or any(block[c] != counts[c] for c in block):
        raise ValueError('tile centres do not form the expected equatorial block')
    neighbours = {c: set(map(int, hp.get_all_neighbours(4, c, nest=True))) & counts.keys() for c in counts}
    north = {max(counts, key=counts.get)}
    todo = list(north)
    while todo:
        new = neighbours[todo.pop()] - north
        north.update(new)
        todo.extend(new)
    south = set(counts) - north
    if not set(block) <= south:
        raise ValueError('equatorial block is not disconnected from northern component')
    vec = lambda c: np.asarray(hp.pix2vec(4, c, nest=True))
    seeds = [max(sorted(north), key=lambda c: np.dot(hp.ang2vec(r, 55, lonlat=True), vec(c))) for r in (270, 100)]
    own = {c: 0 for c in south}
    own.update({seeds[0]: 0, seeds[1]: 1})
    totals = [sum(counts[c] for c in south) + counts[seeds[0]], counts[seeds[1]]]
    free = north - set(seeds)
    target = fraction * sum(counts.values())
    while free and sum(totals) < target:
        i = min(range(2), key=lambda j: (totals[j] / (.4, .6)[j], j))
        region = [c for c in sorted(north) if own.get(c) == i]
        front = free.intersection(set().union(*(neighbours[c] for c in region)))
        if front:
            c = max(front, key=lambda c: (np.dot(vec(seeds[i]), vec(c)), counts[c], -c))
        else:
            c = max(free, key=lambda c: (max(np.dot(vec(c), vec(p)) for p in region), counts[c], -c))
        own[c] = i
        totals[i] += counts[c]
        free.remove(c)
    if fraction < 1:
        return {k: ('nibi', 'fir')[v] for k, v in own.items()}
    used = np.array(sorted(own))
    missing = np.array([p for p in range(192) if p not in own])
    if len(missing):
        vu = np.array(hp.pix2vec(4, used, nest=True)).T
        vm = np.array(hp.pix2vec(4, missing, nest=True)).T
        _, indices = cKDTree(vu).query(vm)
        for p, j in zip(missing, indices):
            own[int(p)] = own[int(used[j])]
    return {k: ('nibi', 'fir')[v] for k, v in own.items()}


def propose(cells, cluster):
    """Propose an unclaimed adjacent cell nearest the cluster's seed."""
    if cluster not in ('nibi', 'fir'):
        raise ValueError('unknown cluster')
    region = {c for c, owner in cells.items() if owner == cluster}
    if not region:
        raise ValueError('cluster has no claimed cells')
    front = {int(n) for c in region for n in hp.get_all_neighbours(4, c, nest=True) if n >= 0 and n not in cells}
    if not front:
        raise ValueError('no unclaimed adjacent cell')
    seed = hp.ang2vec(270 if cluster == 'nibi' else 100, 55, lonlat=True)
    return max(front, key=lambda c: (np.dot(seed, hp.pix2vec(4, c, nest=True)), -c))


def closure(cells, cluster, pixel_tiles, pixel_exposures, index):
    """Find owned-pixel tiles and exposures, including their input exposures."""
    def select(path, field):
        with open(path, newline='') as f:
            rows = list(csv.DictReader(f))
        result = set()
        for row in rows:
            pixel = int(row['pixel'])
            if not 0 <= pixel < 12288:
                raise ValueError(f'invalid pixel {pixel}')
            owner = cells.get(pixel >> 6)
            if owner == cluster:
                result.add(row[field])
        return result
    tiles = set(validate_tiles(select(pixel_tiles, 'tile_id')))
    exps = select(pixel_exposures, 'exp_id')
    edges = tile_exposure_index(tiles, index)
    exps.update(e for values in edges.values() for e in values)
    return sorted(tiles), sorted(exps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('initial'); p.add_argument('--tiles', required=True); p.add_argument('--output', required=True); p.add_argument('--fraction', type=float, default=0.6)
    p = commands.add_parser('propose'); p.add_argument('--claims', required=True); p.add_argument('--cluster', choices=('nibi', 'fir'), required=True)
    p = commands.add_parser('closure'); p.add_argument('--claims', required=True); p.add_argument('--cluster', choices=('nibi', 'fir'), required=True); p.add_argument('--pixel-tiles', required=True); p.add_argument('--pixel-exposures', required=True); p.add_argument('--index', required=True); p.add_argument('--tiles-out', required=True); p.add_argument('--exposures-out', required=True)
    args = parser.parse_args()
    if args.command == 'initial':
        Path(args.output).write_text(yaml.safe_dump({'order': 2, 'cells': initial(read_tiles(args.tiles), args.fraction)}, sort_keys=True))
    elif args.command == 'propose':
        print(propose(load_claims(args.claims), args.cluster))
    else:
        tiles, exps = closure(load_claims(args.claims), args.cluster, args.pixel_tiles, args.pixel_exposures, args.index)
        Path(args.tiles_out).write_text(''.join(f'{x}\n' for x in tiles))
        Path(args.exposures_out).write_text(''.join(f'{x}\n' for x in exps))


if __name__ == '__main__':
    main()
