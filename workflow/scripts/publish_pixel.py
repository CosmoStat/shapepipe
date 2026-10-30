#!/usr/bin/env python3
"""Build and publish one owned pixel's catalogues and maps to VOSpace."""

import argparse
from pathlib import Path
import tempfile

import yaml

from build_hats_pixel import build_pixel, write_catalogue_metadata
from build_map_pixel import build_maps
from index_vos import reindex
from publish_vos import exists, upload_once


class UniqueLoader(yaml.SafeLoader):
    """Reject duplicate ownership claims rather than silently choosing one."""


def _mapping(loader, node):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node)
        if key in result:
            raise ValueError(f'Duplicate claim: {key}')
        result[key] = loader.construct_object(value_node)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def require_owner(pixel, cluster, claims):
    """Refuse invalid, unclaimed or differently owned pixels."""
    data = yaml.load(Path(claims).read_text(), Loader=UniqueLoader)
    if data.get('order') != 2 or cluster not in ('nibi', 'fir'):
        raise ValueError('Invalid claims order or cluster')
    if not 0 <= pixel < 12 * 4 ** 5 or data.get('cells', {}).get(pixel >> 6) != cluster:
        raise ValueError(f'{cluster} does not own pixel {pixel}')


def publish(client, args):
    """Build, upload and index the catalogues, then publish the maps."""
    require_owner(args.pixel, args.cluster, args.claims)
    output = Path(args.output)
    root = args.root.rstrip('/')
    for kind, name, source, schema in [
        ('tiles', 'galaxies', args.tile_list, args.galaxy_schema_from),
        ('stars', 'stars', args.exposure_list, args.star_schema_from),
    ]:
        catalogue = output / name
        sample = build_pixel(args.pixel, kind, source, args.products_dir, args.run_dir,
                             catalogue, args.sha, args.mask, schema_from=schema)
        remote = root + '/' + name
        upload_once(client, sample, remote + '/' + sample.relative_to(catalogue).as_posix())
        with tempfile.TemporaryDirectory() as tmp:
            write_catalogue_metadata(tmp, sample, kind, name)
            for relative in ('properties', 'dataset/_common_metadata'):
                uri = remote + '/' + relative
                if not exists(client, uri):
                    upload_once(client, Path(tmp) / relative, uri)
        reindex(client, remote)
    maps = build_maps(args.pixel, args.exposure_list, args.products_dir, args.run_dir,
                      output / 'maps', args.sha, args.mask)
    for name, path in maps.items():
        upload_once(client, path, root + '/maps/' + name + '/' + path.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pixel', type=int, required=True)
    parser.add_argument('--cluster', choices=['nibi', 'fir'], required=True)
    for name in ('claims', 'tile-list', 'exposure-list', 'products-dir', 'run-dir', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--root', required=True, help='VOSpace publication root')
    parser.add_argument('--sha', required=True)
    parser.add_argument('--mask', type=Path, action='append', required=True)
    parser.add_argument('--galaxy-schema-from', type=Path)
    parser.add_argument('--star-schema-from', type=Path)
    args = parser.parse_args()
    from vos import Client
    publish(Client(), args)


if __name__ == '__main__':
    main()
