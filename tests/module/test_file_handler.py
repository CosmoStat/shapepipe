"""UNIT TESTS FOR FILE_HANDLER.

This module contains unit tests for the shapepipe.pipeline.file_handler
module, in particular the early validation of NUMBER_LIST entries
against the input file numbers actually found on disk (#746): a typo in
NUMBER_LIST must fail at start-up, not when a module first tries to
open the non-existent files.

:Author: Claude (on behalf of Cail Daley) <cail.daley@cea.fr>

"""

import numpy as np
import pytest

from shapepipe.pipeline.file_handler import FileHandler

RE_PATTERN = r"-\d{3}-\d{3}"
NUM_SCHEME = f"RE:{RE_PATTERN}"


def _format_process_list(tmp_path, number_list, scanned):
    """Run FileHandler._format_process_list on a bare instance."""
    handler = FileHandler.__new__(FileHandler)
    handler._number_list = number_list
    handler._verbose = False

    memory_map = str(tmp_path / "match.npy")
    np.save(memory_map, np.array(scanned))

    return handler._format_process_list(
        [("image", ".fits", str(tmp_path))],
        memory_map,
        RE_PATTERN,
        NUM_SCHEME,
        "parallel",
    )


def test_number_list_subset_of_scanned_passes(tmp_path):
    """NUMBER_LIST entries found on disk yield exactly those processes."""
    process_list = _format_process_list(
        tmp_path,
        number_list=["-277-282"],
        scanned=["-277-282", "-300-301"],
    )
    assert process_list == [
        ["-277-282", f"{tmp_path}/image-277-282.fits"],
    ]


def test_number_list_typo_fails_early(tmp_path):
    """A NUMBER_LIST entry absent from disk must raise at start-up."""
    with pytest.raises(ValueError, match="NUMBER_LIST"):
        _format_process_list(
            tmp_path,
            number_list=["-999-999"],
            scanned=["-277-282"],
        )


def test_no_number_list_uses_scanned(tmp_path):
    """Without NUMBER_LIST the scanned numbers drive the process list."""
    process_list = _format_process_list(
        tmp_path,
        number_list=None,
        scanned=["-277-282"],
    )
    assert process_list == [
        ["-277-282", f"{tmp_path}/image-277-282.fits"],
    ]
