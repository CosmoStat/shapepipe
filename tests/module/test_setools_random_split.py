"""Unit tests for SETools' file-number-seeded star split."""

import numpy as np
import pytest

from shapepipe.modules.setools_package.setools import SETools


def _split_for_global_seed(selected_catalogue, file_number, global_seed):
    """Build the star split after setting NumPy's unrelated global state."""
    np.random.seed(global_seed)
    splitter = object.__new__(SETools)
    splitter._cat_size = len(selected_catalogue)
    splitter._file_number_string = file_number
    splitter._rand_split = {
        "star_split": ["RATIO=20", "MASK=selected"],
    }
    splitter.mask = {"selected": selected_catalogue["selected"]}
    splitter._make_rand_split()
    return splitter.rand_split["star_split"]


@pytest.mark.decision("star_selection_psf.psf_train_validation_split")
def test_star_split_is_seeded_disjoint_and_exhaustive():
    """The file number alone determines the 20/80 selected-star split."""
    n_selected = 37
    selected_catalogue = np.zeros(n_selected, dtype=[("selected", bool)])
    selected_catalogue["selected"] = True
    file_number = "-123456-07"

    global_state = np.random.get_state()
    try:
        first = _split_for_global_seed(selected_catalogue, file_number, 19)
        second = _split_for_global_seed(selected_catalogue, file_number, 9876)
    finally:
        np.random.set_state(global_state)

    ratio_20 = first["ratio_20"]
    ratio_80 = first["ratio_80"]
    np.testing.assert_array_equal(ratio_20, second["ratio_20"])
    np.testing.assert_array_equal(ratio_80, second["ratio_80"])
    assert len(ratio_20) == int(np.ceil(0.2 * n_selected))
    assert set(ratio_20).isdisjoint(ratio_80)
    np.testing.assert_array_equal(
        np.sort(np.concatenate((ratio_20, ratio_80))),
        np.arange(n_selected),
    )
