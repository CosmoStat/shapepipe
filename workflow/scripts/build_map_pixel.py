#!/usr/bin/env python3
"""Sum exposure fragments into one order-5 pixel per coverage map."""

import argparse
from pathlib import Path

import healsparse as hsp
import numpy as np

from build_hats_pixel import HATS_ORDER, metadata, write_once
from exp_maps import FLAG_SHIFT, FRAGMENT_DTYPE, NSIDE, NSIDE_COVERAGE
from merge_exposure_maps import add, fragment_path
from publication_sources import check_exposure_pixels, finished_sources


# @sc [decision:survey_publication.row_partitioning,decision:survey_publication.pixel_completeness]
def build_maps(pixel, source_list, products_dir, run_dir, output, sha, mask):
    """Sum only the listed finished exposures and restrict to this pixel."""
    if not 0 <= pixel < 12 * 4 ** HATS_ORDER:
        raise ValueError(f'Invalid order-5 pixel {pixel}')
    units = finished_sources(pixel, 'maps', source_list, products_dir, run_dir)
    meta = metadata(sha, mask)
    maps = {name: hsp.HealSparseMap.make_empty(NSIDE_COVERAGE, NSIDE, dtype,
                                                metadata=meta)
            for name, dtype in [('npoint', np.uint8), ('ntaken', np.uint8),
                                ('nexp', np.uint8), ('nflagged', np.uint16)]}
    for unit in units:
        path = fragment_path(Path(products_dir), unit)
        fragment = hsp.HealSparseMap.read(str(path))
        if fragment.dtype != FRAGMENT_DTYPE or fragment.nside_sparse != NSIDE:
            raise ValueError(f'Exposure {unit}: incompatible fragment {path}')
        pixels = fragment.valid_pixels
        parents = pixels >> (2 * (17 - HATS_ORDER))
        check_exposure_pixels(source_list, unit, np.unique(parents))
        pixels = pixels[parents == pixel]
        values = fragment.get_values_pix(pixels)
        for name, bit in [('npoint', 1), ('ntaken', 2), ('nexp', 4)]:
            selected = pixels[(values & bit) != 0]
            add(maps[name], selected, np.ones(len(selected), dtype=np.uint8))
        counts = values >> FLAG_SHIFT
        selected = counts > 0
        add(maps['nflagged'], pixels[selected], counts[selected])
    output = Path(output)
    paths = {}
    for name, acc in maps.items():
        path = output / name / f'Npix={pixel}.hsp'
        write_once(path, lambda tmp, m=acc: m.write(str(tmp), clobber=True))
        paths[name] = path
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pixel', type=int, required=True)
    for arg in ('source-list', 'products-dir', 'run-dir', 'output'):
        parser.add_argument('--' + arg, type=Path, required=True)
    parser.add_argument('--mask', type=Path, action='append', required=True,
                        help='Configured mask file; repeat for every mask')
    parser.add_argument('--sha', required=True)
    args = parser.parse_args()
    build_maps(args.pixel, args.source_list, args.products_dir, args.run_dir,
               args.output, args.sha, args.mask)


if __name__ == '__main__':
    main()
