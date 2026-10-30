#!/usr/bin/env python3
"""Build one fixed-order HATS pixel from a run's durable products.

The source list is the planning output, not a directory scan. Pass every
configured external mask file; each is read only to compute its SHA-256.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from hats.io.file_io import get_parquet_write_table_kwargs
from hats.pixel_math import compute_spatial_index
from shapepipe.utilities.final_cat import read_final_cat

import merge_star_cat
from publication_sources import check_exposure_pixels, finished_sources

# @sc [decision:survey_publication.row_partitioning]
HATS_ORDER = 5


def metadata(sha, mask):
    """Return provenance, hashing each external mask without copying it.

    @sc [label:provenance] code-sha-in-file
    Every output carries the ShapePipe code SHA and the paths and SHA-256s of
    the external masks. tests/unit/test_publication_pixel.py reads them back.
    """
    if not sha:
        raise ValueError("ShapePipe code SHA is required")
    masks = [mask] if isinstance(mask, (str, Path)) else list(mask)
    if not masks:
        raise ValueError('At least one configured mask is required')
    records = []
    for path in masks:
        digest = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        records.append({'path': str(path), 'sha256': digest.hexdigest()})
    return {'SPSHA': sha, 'MASKS': json.dumps(records, separators=(',', ':'))}


def write_once(path, writer):
    """Publish complete files without overwriting an existing destination."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.publication-', dir=path.parent, suffix=path.suffix)
    os.close(fd)
    try:
        writer(Path(tmp))
        os.link(tmp, path)
    finally:
        Path(tmp).unlink(missing_ok=True)


def pixel_path(root, pixel):
    """Return the HATS order-5 path for a pixel."""
    return Path(root) / 'dataset' / 'Norder=5' / f'Dir={pixel // 10000 * 10000}' / f'Npix={pixel}.parquet'


def _array_table(data):
    columns = {}
    for name in data.dtype.names:
        values = data[name]
        values = values.astype(values.dtype.newbyteorder('='), copy=False)
        if values.ndim == 1:
            columns[name] = pa.array(values, from_pandas=False)
        else:
            element_type = pa.from_numpy_dtype(values.dtype)
            for _ in values.shape[1:]:
                element_type = pa.list_(element_type)
            columns[name] = pa.array(values.tolist(), type=element_type, from_pandas=False)
    return pa.table(columns)


def _star_rows(products, exp):
    manifest = products / 'exp' / exp[:2] / exp / 'manifests' / 'exp_persist.json'
    record = json.loads(manifest.read_text())
    data = merge_star_cat.read_exposure(exp, Path(record['tar']))
    table = _array_table(data)
    # The planning branch adds these to read_exposure. Until rebased, derive
    # the same key from the per-CCD contiguous rows it returns.
    if 'EXP_ID' not in table.column_names:
        table = table.append_column('EXP_ID', pa.array([exp] * len(data)))
    if 'SOURCE_ORDINAL' not in table.column_names:
        ccd = data['CCD_NB']
        ordinal = np.empty(len(data), dtype=np.int64)
        counts = {}
        for i, value in enumerate(ccd):
            key = int(value)
            ordinal[i] = counts.get(key, 0)
            counts[key] = ordinal[i] + 1
        table = table.append_column('SOURCE_ORDINAL', pa.array(ordinal))
    return table


def build_pixel(pixel, kind, source_list, products_dir, run_dir, output, sha, mask,
                columns=None, schema_from=None):
    """Build a pixel after all sources have succeeded or been abandoned.

    @sc [decision:survey_publication.row_partitioning,decision:survey_publication.pixel_completeness,label:selection] pixel-rows-exact
    Each listed successful source contributes every row whose own position is
    in this pixel exactly once. tests/unit/test_publication_pixel.py compares
    the reassembled pixels with the input rows, including overlaps.
    """
    if pixel < 0 or pixel >= 12 * 4 ** HATS_ORDER:
        raise ValueError(f'Invalid order-5 pixel {pixel}')
    if kind not in ('tiles', 'stars'):
        raise ValueError(f'Unknown catalogue kind {kind}')
    products = Path(products_dir)
    units = finished_sources(pixel, kind, source_list, products, run_dir)
    ra, dec = ('XWIN_WORLD', 'YWIN_WORLD') if kind == 'tiles' else ('RA', 'DEC')
    pieces = []
    for unit in units:
        if kind == 'tiles':
            path = products / 'tiles' / unit[:2] / unit / f'final_cat-{unit}.hdf5'
            table = _array_table(read_final_cat(path, columns))
        else:
            table = _star_rows(products, unit)
        spatial = compute_spatial_index(table[ra].to_numpy(), table[dec].to_numpy())
        if kind == 'stars':
            check_exposure_pixels(source_list, unit, np.unique(spatial >> 48))
        table = table.append_column('_healpix_29', pa.array(spatial, type=pa.int64()))
        pieces.append(table.filter(pa.array((spatial >> 48) == pixel)))
    if pieces:
        result = pa.concat_tables(pieces).sort_by([('_healpix_29', 'ascending')])
    elif schema_from:
        schema = pq.read_schema(schema_from)
        result = pa.Table.from_batches([], schema=schema)
    else:
        raise ValueError('No successful sources; --schema-from is required for an empty pixel')
    result = result.replace_schema_metadata({k.encode(): v.encode() for k, v in metadata(sha, mask).items()})
    path = pixel_path(output, pixel)
    write_once(path, lambda tmp: pq.write_table(result, tmp, **get_parquet_write_table_kwargs(
        result.schema, write_table_kwargs={'row_group_size': 32768})))
    return path


def write_catalogue_metadata(root, sample, kind, name):
    """Write HATS properties and common schema once, without row counts."""
    if kind not in ('tiles', 'stars'):
        raise ValueError(kind)
    ra, dec = ('XWIN_WORLD', 'YWIN_WORLD') if kind == 'tiles' else ('RA', 'DEC')
    schema = pq.read_schema(sample).remove_metadata()
    props = (f'#HATS catalog\nobs_collection={name}\ndataproduct_type=object\n'
             f'hats_col_ra={ra}\nhats_col_dec={dec}\nhats_col_healpix=_healpix_29\n'
             'hats_col_healpix_order=29\nhats_npix_suffix=.parquet\nhats_order=5\n')
    root = Path(root)
    write_once(root / 'properties', lambda p: p.write_text(props))
    write_once(root / 'dataset' / '_common_metadata', lambda p: pq.write_metadata(schema, p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('build')
    build.add_argument('--pixel', type=int, required=True)
    build.add_argument('--kind', choices=['tiles', 'stars'], required=True)
    for arg in ('source-list', 'products-dir', 'run-dir', 'output'):
        build.add_argument('--' + arg, type=Path, required=True)
    build.add_argument('--mask', type=Path, action='append', required=True,
                       help='Configured mask file; repeat for every mask')
    build.add_argument('--sha', required=True)
    build.add_argument('--param-file', type=Path)
    build.add_argument('--schema-from', type=Path)
    info = sub.add_parser('metadata')
    info.add_argument('--root', type=Path, required=True)
    info.add_argument('--sample', type=Path, required=True)
    info.add_argument('--kind', choices=['tiles', 'stars'], required=True)
    info.add_argument('--name', required=True)
    args = parser.parse_args()
    if args.command == 'metadata':
        write_catalogue_metadata(args.root, args.sample, args.kind, args.name)
    else:
        columns = None
        if args.param_file:
            columns = [line.split('#', 1)[0].strip() for line in args.param_file.read_text().splitlines()]
            columns = [line for line in columns if line]
        build_pixel(args.pixel, args.kind, args.source_list, args.products_dir, args.run_dir,
                    args.output, args.sha, args.mask, columns, args.schema_from)


if __name__ == '__main__':
    main()
