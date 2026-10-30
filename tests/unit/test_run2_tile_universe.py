"""Run-2 tile selection against the reconstructed original DR6 list."""

import sys
from pathlib import Path

import healpy as hp
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
from tile_universe import included, read_tiles, validate_tiles

FIXTURES = Path(__file__).resolve().parents[1] / 'fixtures/run2'
# Sixteen centre-wrap tiles that #969's whole-footprint skip list retained.
EXTRA = set('579.253 590.250 597.248 604.246 611.244 614.243 624.240 '
            '627.239 630.238 633.237 636.236 639.235 642.234 650.231 '
            '653.230 658.228'.split())


@pytest.mark.decision('survey_publication.tile_universe')
def test_rule_reproduces_dr6_reference():
    old_included = set(read_tiles(FIXTURES / 'included_reference.txt'))
    old_excluded = set(read_tiles(FIXTURES / 'excluded_reference.txt'))
    all_tiles = old_included | old_excluded
    assert len(all_tiles) == 23114
    assert {t for t in all_tiles if not included(t)} == old_excluded | EXTRA
    assert len({t for t in all_tiles if included(t)}) == 22885


def test_batch_refuses_excluded_tile():
    validate_tiles(['001.200'])
    with pytest.raises(ValueError, match='579.253'):
        validate_tiles(['001.200', '579.253'])


def test_dateline_caps_lose_no_order5_pixels():
    tiles = set(read_tiles(FIXTURES / 'included_reference.txt')) | set(read_tiles(FIXTURES / 'excluded_reference.txt'))
    def cap(tile):
        x, y = map(int, tile.split('.'))
        dec = y / 2 - 90
        ra = x / (2 * np.cos(np.deg2rad(dec)))
        return set(hp.query_disc(32, hp.ang2vec(ra % 360, dec, lonlat=True),
                                 np.deg2rad(.4), inclusive=True, fact=16, nest=True))
    kept = set().union(*(cap(t) for t in tiles if included(t)))
    dateline = set().union(*(cap(t) for t in tiles if not included(t) and int(t.split('.')[1]) < 345))
    assert dateline <= kept, sorted(dateline - kept)
