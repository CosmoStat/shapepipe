"""Generate and validate the rule-defined run-2 tile universe."""

import argparse
import math
from pathlib import Path

# @sc [decision:survey_publication.tile_universe]
CAP_Y_MIN = 345


def included(tile):
    """Return whether a CFIS tile belongs to the run-2 catalogue.

    @sc [decision:survey_publication.tile_universe,label:invariant] tile-universe-by-rule
    The catalogue excludes polar-cap and dateline-redundant tiles by their
    unwrapped tile ID before planning, rather than by cluster-specific skips.
    """
    x, y = map(int, tile.split('.'))
    return y < CAP_Y_MIN and x < 720 * math.cos(math.radians(y / 2 - 90))


def read_tiles(path):
    """Read tile IDs, ignoring blank and comment lines."""
    return list(dict.fromkeys(line.strip() for line in Path(path).read_text().splitlines()
                              if line.strip() and not line.lstrip().startswith('#')))


def validate_tiles(tiles):
    """Refuse a plan containing any excluded tile."""
    excluded = [tile for tile in tiles if not included(tile)]
    if excluded:
        raise ValueError(f"tile universe excludes {len(excluded)} tile(s): "
                         + ', '.join(excluded[:20]))
    return tiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='DR6 tile list')
    parser.add_argument('output', type=Path, help='sorted run-2 tile list')
    args = parser.parse_args()
    tiles = read_tiles(args.input)
    selected = sorted(tile for tile in tiles if included(tile))
    args.output.write_text(''.join(f'{tile}\n' for tile in selected))
    print(f'Excluded {len(tiles) - len(selected)} of {len(tiles)} tiles; wrote {args.output}')


if __name__ == '__main__':
    main()
