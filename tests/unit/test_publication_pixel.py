"""Synthetic publication patch, including boundaries and source gating."""

import hashlib
import sys
from pathlib import Path

import h5py
import healsparse as hsp
import numpy as np
import pyarrow.parquet as pq
import pytest
from hats.pixel_math import compute_spatial_index

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
import build_hats_pixel as builder
import build_map_pixel as maps


@pytest.fixture
def patch(tmp_path):
    products = tmp_path / 'products'
    mask = tmp_path / 'external.hsp'
    mask.write_bytes(b'external mask')
    # Points are spread over several order-5 pixels; each tile has overlapping
    # sky with its neighbours and distinct values, which must not be deduped.
    coords = [(10., 10.), (11., 10.), (12., 10.), (13., 10.)]
    pixels = [int(x >> 48) for x in compute_spatial_index(
        np.array([x for x, y in coords]), np.array([y for x, y in coords]))]
    tiles = ['100.200', '101.200', '102.200', '103.200']
    expected = []
    for n, tile in enumerate(tiles):
        d = products / 'tiles' / tile[:2] / tile
        d.mkdir(parents=True)
        (d / 'tile_make_cat.json').write_text('{}')
        with h5py.File(d / f'final_cat-{tile}.hdf5', 'w') as f:
            f['XWIN_WORLD'] = [coords[n][0], coords[(n + 1) % 4][0]]
            f['YWIN_WORLD'] = [10., 10.]
            f['VALUE'] = [n * 2, n * 2 + 1]
        expected.extend([(n * 2, coords[n][0], 10.),
                         (n * 2 + 1, coords[(n + 1) % 4][0], 10.)])
    listing = tmp_path / 'pixel_tiles.csv'
    listing.write_text('pixel,tile_id\n' + ''.join(
        f'{p},{tile}\n' for p in sorted(set(pixels)) for tile in tiles))
    return tmp_path, products, mask, listing, sorted(set(pixels)), expected


def test_pixel_rows_exact(patch):
    root, products, mask, listing, pixels, expected = patch
    seen = []
    for pixel in pixels:
        path = builder.build_pixel(pixel, 'tiles', listing, products, root,
                                   root / 'hats', 'abc123', mask)
        tab = pq.ParquetFile(path).read()
        hp = tab['_healpix_29'].to_numpy()
        assert np.all(hp >> 48 == pixel)
        assert np.all(hp[:-1] <= hp[1:])
        seen.extend(zip(tab['VALUE'].to_pylist(), tab['XWIN_WORLD'].to_pylist(),
                        tab['YWIN_WORLD'].to_pylist()))
    assert sorted(seen) == sorted(expected)
    builder.write_catalogue_metadata(root / 'hats', path, 'tiles', 'test')
    assert 'hats_nrows' not in (root / 'hats/properties').read_text()
    assert pq.read_schema(root / 'hats/dataset/_common_metadata').names == tab.schema.names
    with pytest.raises(FileExistsError):
        builder.write_catalogue_metadata(root / 'hats', path, 'tiles', 'test')


def test_missing_source_then_abandoned(patch):
    root, products, mask, listing, pixels, _ = patch
    missing = '100.200'
    (products / 'tiles/10' / missing / 'tile_make_cat.json').unlink()
    with pytest.raises(RuntimeError, match=missing):
        builder.build_pixel(pixels[0], 'tiles', listing, products, root, root / 'hats', 'abc', mask)
    (root / 'abandoned.txt').write_text(missing + ' # no retry\n')
    p = builder.build_pixel(pixels[0], 'tiles', listing, products, root, root / 'hats', 'abc', mask)
    assert not set(pq.ParquetFile(p).read()['VALUE'].to_pylist()) & {0, 1}


def test_code_sha_in_file(patch):
    root, products, mask, listing, pixels, _ = patch
    p = builder.build_pixel(pixels[0], 'tiles', listing, products, root, root / 'hats', 'abc123', mask)
    meta = pq.read_metadata(p).metadata
    assert meta[b'SPSHA'] == b'abc123'
    assert meta[b'MASKPATH'] == str(mask).encode()
    assert meta[b'MASKSHA'] == hashlib.sha256(mask.read_bytes()).hexdigest().encode()


def test_star_partition_and_unlisted_pixel(patch, monkeypatch):
    root, products, mask, _, pixels, _ = patch
    coords = np.array([10., 11., 12., 13.])
    sp = compute_spatial_index(coords, np.full(4, 10.)) >> 48
    import pyarrow as pa
    monkeypatch.setattr(builder, '_star_rows', lambda products, exp: pa.table(
        {'RA': coords, 'DEC': np.full(4, 10.), 'EXP_ID': [exp] * 4,
         'CCD_NB': [1] * 4, 'SOURCE_ORDINAL': np.arange(4)}))
    manifest = products / 'exp/12/123456/manifests/exp_persist.json'
    manifest.parent.mkdir(parents=True)
    manifest.write_text('{}')
    listing = root / 'exposures.csv'
    listing.write_text('pixel,exp_id\n' + ''.join(f'{p},123456\n' for p in pixels))
    collected = []
    for p in pixels:
        path = builder.build_pixel(p, 'stars', listing, products, root, root / 'stars', 'abc', mask)
        collected.extend(pq.ParquetFile(path).read()['SOURCE_ORDINAL'].to_pylist())
    assert sorted(collected) == list(range(4))
    listing.write_text('pixel,exp_id\n' + ''.join(f'{p},123456\n' for p in pixels[:-1]))
    with pytest.raises(ValueError, match=f'123456.*{pixels[-1]}'):
        builder.build_pixel(pixels[0], 'stars', listing, products, root, root / 'other', 'abc', mask)


def test_stars_read_persisted_tar(tmp_path):
    import io
    import json
    import tarfile
    from astropy.io import fits
    import merge_star_cat

    products = tmp_path / 'products'
    base = products / 'exp/12/123456'
    (base / 'psf').mkdir(parents=True)
    (base / 'manifests').mkdir()
    tar_path = base / 'psf/123456.tar'
    with tarfile.open(tar_path, 'w') as tar:
        for ccd in (2, 5):
            columns = [fits.Column(name=name, format='D', array=[10., 11.])
                       for name in merge_star_cat.COLUMNS]
            data = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(),
                                 fits.BinTableHDU.from_columns(columns)])
            buf = io.BytesIO()
            data.writeto(buf)
            raw = buf.getvalue()
            info = tarfile.TarInfo(f'validation_psf-123456-{ccd}.fits')
            info.size = len(raw)
            tar.addfile(info, io.BytesIO(raw))
    (base / 'manifests/exp_persist.json').write_text(json.dumps({'tar': str(tar_path)}))
    rows = builder._star_rows(products, '123456')
    assert rows['EXP_ID'].to_pylist() == ['123456'] * 4
    assert rows['CCD_NB'].to_pylist() == [2, 2, 5, 5]
    assert rows['SOURCE_ORDINAL'].to_pylist() == [0, 1, 0, 1]


def test_maps_match_sum_and_provenance(tmp_path, monkeypatch):
    products = tmp_path / 'products'
    mask = tmp_path / 'mask.hsp'
    mask.write_bytes(b'mask')
    pixel1, pixel2 = 100, 101
    shift = 24
    children = np.array([(pixel1 << shift) + 4, (pixel2 << shift) + 9], dtype=np.int64)
    listing = tmp_path / 'exposures.csv'
    listing.write_text('pixel,exp_id\n100,123\n101,123\n100,124\n101,124\n')
    for exp, values in [('123', [7 | (2 << 3), 3 | (1 << 3)]),
                        ('124', [5 | (3 << 3), 7 | (4 << 3)])]:
        d = products / 'exp' / exp[:2] / exp
        (d / 'manifests').mkdir(parents=True)
        (d / 'maps').mkdir()
        (d / 'manifests/exp_maps.json').write_text('{}')
        frag = hsp.HealSparseMap.make_empty(128, 131072, np.uint16)
        frag[children] = np.array(values, dtype=np.uint16)
        frag.write(str(d / 'maps' / f'maps-{exp}.hsp'))
    import merge_exposure_maps as merge
    monkeypatch.setattr(merge.build_index, 'campaign_exposures', lambda *args: ['123', '124'])
    merged = {name: tmp_path / f'{name}.hsp' for name in ('npoint', 'ntaken', 'nexp', 'nflagged')}
    monkeypatch.setattr(sys, 'argv', ['merge_exposure_maps.py', '--products-dir', str(products),
                                     '--tile-list', 'unused', '--index-db', 'unused',
                                     *sum(([f'--{name}', str(path)] for name, path in merged.items()), [])])
    merge.main()
    for i, p in enumerate([pixel1, pixel2]):
        result = maps.build_maps(p, listing, products, tmp_path, tmp_path / 'out', 'abc123', mask)
        for name, path in result.items():
            m = hsp.HealSparseMap.read(str(path))
            whole = hsp.HealSparseMap.read(str(merged[name]))
            assert m.get_values_pix(children[i]) == whole.get_values_pix(children[i])
            assert set(m.valid_pixels) <= {int(children[i])}
            assert m.metadata['SPSHA'] == 'abc123'
            assert m.metadata['MASKPATH'] == str(mask)
            assert m.metadata['MASKSHA'] == hashlib.sha256(mask.read_bytes()).hexdigest()
    listing.write_text('pixel,exp_id\n100,123\n100,124\n101,124\n')
    with pytest.raises(ValueError, match='123.*101'):
        maps.build_maps(100, listing, products, tmp_path, tmp_path / 'other', 'abc', mask)
