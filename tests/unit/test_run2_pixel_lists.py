"""Check conservative candidate lists and real tile-row enclosure."""

import csv
import sqlite3
import sys
from pathlib import Path

import healpy as hp
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
import pixel_lists as p


@pytest.mark.decision('survey_publication.pixel_completeness')
def test_real_tile_rows_within_cap():
    path = Path(__file__).resolve().parents[1] / 'fixtures/run2/real_rows.csv'
    with path.open() as f:
        rows = list(csv.DictReader(f))
    separations = []
    for row in rows:
        centre = hp.ang2vec(*p.tile_centre(row['tile_id']), lonlat=True)
        point = hp.ang2vec(float(row['ra']), float(row['dec']), lonlat=True)
        separations.append(np.rad2deg(hp.rotator.angdist(centre, point))[0])
    assert len({r['tile_id'] for r in rows}) == 4
    assert max(separations) < p.TILE_CAP_DEG


def make_index(db, rows):
    with sqlite3.connect(db) as conn:
        conn.execute('CREATE TABLE tiles (tile_id TEXT, ra_dir TEXT, n_exp INTEGER)')
        conn.execute('CREATE TABLE tile_exposures (tile_id TEXT, exp_id TEXT)')
        conn.executemany('INSERT INTO tiles VALUES (?, ?, ?)',
                         [(tile, tile[:3], 1) for tile, _ in rows])
        conn.executemany('INSERT INTO tile_exposures VALUES (?, ?)', rows)


@pytest.mark.decision('survey_publication.pixel_completeness')
def test_exposure_list_covers_star_rows(tmp_path):
    path = Path(__file__).resolve().parents[1] / 'fixtures/run2/real_rows.csv'
    with path.open() as f:
        rows = list(csv.DictReader(f))
    tiles = sorted({r['tile_id'] for r in rows})
    db = tmp_path / 'index.sqlite'
    make_index(db, [(t, 'e' + t) for t in tiles])
    pairs = p.exposure_pairs(list(p.tile_pairs(tiles)), db)
    for row in rows:
        pixel = int(hp.ang2pix(32, float(row['ra']), float(row['dec']), nest=True, lonlat=True))
        assert (pixel, 'e' + row['tile_id']) in pairs


def test_missing_or_stale_associations_rejected(tmp_path):
    db = tmp_path / 'index.sqlite'
    make_index(db, [('001.200', 'e1')])
    with sqlite3.connect(db) as conn:
        conn.execute("INSERT INTO tile_exposures VALUES ('002.200', 'stale')")
    with pytest.raises(ValueError, match='002.200'):
        p.exposure_pairs(list(p.tile_pairs(['001.200', '002.200'])), db)


def test_cli_missing_tile_writes_neither_list(tmp_path, monkeypatch):
    db = tmp_path / 'index.sqlite'
    make_index(db, [('001.200', 'e1')])
    tiles = tmp_path / 'tiles.txt'
    tiles.write_text('001.200\n002.200\n')
    tp, ep = tmp_path / 'tp.csv', tmp_path / 'ep.csv'
    monkeypatch.setattr(sys, 'argv', ['pixel_lists', '--tiles', str(tiles),
                        '--index', str(db), '--pixel-tiles', str(tp),
                        '--pixel-exposures', str(ep)])
    with pytest.raises(ValueError, match='002.200'):
        p.main()
    assert not tp.exists() and not ep.exists()
