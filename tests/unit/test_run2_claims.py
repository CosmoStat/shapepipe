"""Ownership claims and closure use unique order-2 owners."""

import sqlite3
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
import claims


@pytest.mark.decision('survey_publication.seam_realization')
def test_duplicate_cell_rejected(tmp_path):
    path = tmp_path / 'claims.yaml'
    path.write_text('order: 2\ncells:\n  3: nibi\n')
    assert claims.load_claims(path) == {3: 'nibi'}
    path.write_text('order: 2\ncells:\n  3: nibi\n  3: fir\n')
    with pytest.raises(ValueError, match='duplicate'):
        claims.load_claims(path)


def test_closure_includes_owned_pixels_and_tile_exposures(tmp_path):
    tiles = tmp_path / 'pixel_tiles.csv'
    tiles.write_text('pixel,tile_id\n64,001.200\n128,002.200\n')
    exps = tmp_path / 'pixel_exposures.csv'
    exps.write_text('pixel,exp_id\n64,e2\n128,e3\n')
    db = tmp_path / 'index.sqlite'
    with sqlite3.connect(db) as con:
        con.execute('CREATE TABLE tile_exposures (tile_id TEXT, exp_id TEXT)')
        con.execute("INSERT INTO tile_exposures VALUES ('001.200', 'e1')")
    assert claims.closure({1: 'nibi', 2: 'fir'}, 'nibi', tiles, exps, db) == (['001.200'], ['e1', 'e2'])


def test_proposal_is_unclaimed_neighbour():
    cell = claims.propose({39: 'nibi', 27: 'fir'}, 'nibi')
    assert cell not in {39, 27}
    assert cell in claims.hp.get_all_neighbours(4, 39, nest=True)
