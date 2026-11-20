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
        con.execute('CREATE TABLE tiles (tile_id TEXT, ra_dir TEXT, n_exp INTEGER)')
        con.execute("INSERT INTO tiles VALUES ('001.200', '001', 1)")
        con.execute('CREATE TABLE tile_exposures (tile_id TEXT, exp_id TEXT)')
        con.execute("INSERT INTO tile_exposures VALUES ('001.200', 'e1')")
    assert claims.closure({1: 'nibi', 2: 'fir'}, 'nibi', tiles, exps, db) == (['001.200'], ['e1', 'e2'])


def test_proposal_is_unclaimed_neighbour():
    cell = claims.propose({39: 'nibi', 27: 'fir'}, 'nibi')
    assert cell not in {39, 27}
    assert cell in claims.hp.get_all_neighbours(4, 39, nest=True)


def test_closure_rejects_excluded_tile(tmp_path):
    tiles = tmp_path / 'tiles.csv'
    tiles.write_text('pixel,tile_id\n64,579.253\n')
    exps = tmp_path / 'exps.csv'
    exps.write_text('pixel,exp_id\n')
    with pytest.raises(ValueError):
        claims.closure({1: 'nibi'}, 'nibi', tiles, exps, tmp_path / 'missing.sqlite')


def test_closure_rejects_unprepared_tile(tmp_path):
    tiles = tmp_path / 'tiles.csv'
    tiles.write_text('pixel,tile_id\n64,001.200\n')
    exps = tmp_path / 'exps.csv'
    exps.write_text('pixel,exp_id\n')
    db = tmp_path / 'index.sqlite'
    with sqlite3.connect(db) as con:
        con.execute('CREATE TABLE tiles (tile_id TEXT, ra_dir TEXT, n_exp INTEGER)')
        con.execute('CREATE TABLE tile_exposures (tile_id TEXT, exp_id TEXT)')
    with pytest.raises(ValueError, match='001.200'):
        claims.closure({1: 'nibi'}, 'nibi', tiles, exps, db)


def test_partial_initial_can_propose_for_both_clusters():
    from tile_universe import included, read_tiles
    fixture = Path(__file__).resolve().parents[1] / 'fixtures/run2/included_reference.txt'
    tiles = [t for t in read_tiles(fixture) if included(t)]
    cells = claims.initial(tiles)
    for cluster in ('nibi', 'fir'):
        proposal = claims.propose(cells, cluster)
        front = {int(n) for c, owner in cells.items() if owner == cluster
                 for n in claims.hp.get_all_neighbours(4, c, nest=True)
                 if n >= 0 and n not in cells}
        seed = claims.hp.ang2vec(270 if cluster == 'nibi' else 100, 55, lonlat=True)
        assert proposal == max(front, key=lambda c: (claims.np.dot(seed, claims.hp.pix2vec(4, c, nest=True)), -c))
    assert len(claims.initial(tiles, 1)) == 192
