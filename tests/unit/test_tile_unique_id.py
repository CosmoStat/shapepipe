"""Survey-wide object ID: ``tile_id * 10**6 + NUMBER``."""

import numpy as np
import pytest

from shapepipe.utilities import cfis


@pytest.mark.parametrize(
    "name", ["-301-279", "CFIS.301.279.r", "CFIS.301.279.r.weight.fits.fz"]
)
def test_tile_id_from_name(name):
    assert cfis.get_tile_id(name) == 301279


def test_tile_id_keeps_leading_zeros():
    assert cfis.get_tile_id("-004-012") == 4012


@pytest.mark.parametrize("name", ["-1-2", "-1301-279", "-301-2790"])
def test_tile_id_rejects_non_three_digit_components(name):
    with pytest.raises(cfis.CfisError):
        cfis.get_tile_id(name)


def test_unique_id_example():
    assert cfis.get_tile_unique_id(301279, 42) == 301279000042


def test_unique_id_round_trip():
    tile_id = 999999
    number = np.array([0, 1, 12345, 999999])
    unique_id = cfis.get_tile_unique_id(tile_id, number)
    assert unique_id.dtype == np.int64
    assert len(np.unique(unique_id)) == len(number)
    tile_back, number_back = cfis.split_tile_unique_id(unique_id)
    np.testing.assert_array_equal(tile_back, tile_id)
    np.testing.assert_array_equal(number_back, number)


@pytest.mark.parametrize("number", [10**6, -1, [5, 10**6]])
def test_unique_id_rejects_number_out_of_range(number):
    with pytest.raises(cfis.CfisError):
        cfis.get_tile_unique_id(301279, number)


def test_unique_id_rejects_tile_id_out_of_range():
    with pytest.raises(cfis.CfisError):
        cfis.get_tile_unique_id(10**6, 1)
