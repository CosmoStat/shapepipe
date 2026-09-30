"""The FWHM mode that centres the PSF-star box follows the stellar locus.

``mode(FWHM_IMAGE{preselect})`` centres a +-0.2 px window, so its estimator
must land on the stellar locus and must not move when a few objects enter or
leave the preselection at its edges.
"""

import numpy as np
import pytest

from shapepipe.pipeline.str_handler import StrInterpreter

pytestmark = pytest.mark.decision("star_selection_psf.star_selection_box")

LOCUS = 2.85  # px, stellar FWHM
LOCUS_SIGMA = 0.06  # px, locus width measured on 2114045p
BOX_HALF_WIDTH = 0.2  # px, the star_selection.setools window
FWHM_MIN, FWHM_MAX = 0.3 / 0.187, 1.5 / 0.187  # the preselection bounds
SEEDS = range(20)


def mode(values):
    return StrInterpreter._mode(None, np.asarray(values, dtype=np.float32))


def preselection(seed, n_stars=70, n_galaxies=120, n_small=6):
    """A narrow stellar Gaussian, a broad galaxy tail, a few artefacts."""
    rng = np.random.default_rng(seed)
    stars = rng.normal(LOCUS, LOCUS_SIGMA, n_stars)
    galaxies = LOCUS + 0.4 + rng.gamma(2.0, 1.0, n_galaxies)
    small = rng.uniform(FWHM_MIN, LOCUS - 0.4, n_small)
    values = np.concatenate([stars, galaxies, small])
    return values[(values > FWHM_MIN) & (values < FWHM_MAX)]


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_recovers_the_stellar_locus(seed):
    assert abs(mode(preselection(seed)) - LOCUS) < 0.03


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_is_stable_to_objects_added_at_the_edges(seed):
    values = preselection(seed)
    edges = [FWHM_MIN + 0.01, FWHM_MAX - 0.01, FWHM_MAX - 0.02]
    shift = mode(np.concatenate([values, edges])) - mode(values)
    assert abs(shift) < BOX_HALF_WIDTH / 20


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_is_stable_to_objects_removed_at_the_edges(seed):
    values = np.sort(preselection(seed))
    shift = mode(values[2:-1]) - mode(values)
    assert abs(shift) < BOX_HALF_WIDTH / 20


def test_small_samples_fall_back_to_the_median_and_empty_to_minus_one():
    values = np.array([2.7, 2.8, 2.9, 5.0])
    assert mode(values) == pytest.approx(np.median(values))
    assert mode(np.array([], dtype=np.float32)) == -1


def test_duplicated_values_do_not_break_the_estimate():
    """Bootstrap-like resamples repeat values; the peak still lands."""
    values = np.repeat(preselection(0)[:40], 5)
    assert abs(mode(values) - LOCUS) < 0.05
