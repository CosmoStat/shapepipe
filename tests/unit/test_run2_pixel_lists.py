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


def test_pixel_lists_use_tile_exposure_index(tmp_path):
    db = tmp_path / 'index.sqlite'
    with sqlite3.connect(db) as conn:
        conn.execute('CREATE TABLE tile_exposures (tile_id TEXT, exp_id TEXT)')
        conn.executemany('INSERT INTO tile_exposures VALUES (?, ?)',
                         [('001.200', 'e1'), ('002.200', 'e2')])
    caps = tmp_path / 'caps.csv'
    caps.write_text('exp_id,ra,dec,radius_deg\ne1,1,10,0.6\ne2,2,10,0.6\n')
    assert {e for _, e in p.exposure_pairs(['001.200'], db, caps)} == {'e1'}
    assert {t for _, t in p.tile_pairs(['001.200'])} == {'001.200'}
    caps.write_text('exp_id,ra,dec,radius_deg\ne2,2,10,0.6\n')
    with pytest.raises(ValueError, match='missing exposure caps'):
        list(p.exposure_pairs(['001.200'], db, caps))
