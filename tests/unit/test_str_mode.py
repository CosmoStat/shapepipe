"""The PSF-star box follows the stellar locus, in position and in width.

``star_selection.setools`` keeps objects whose FWHM lies within
``3 * locus_width`` of ``mode`` of the preselection. ``mode`` must land on the
stellar locus and not move when a few objects enter or leave the preselection
at its edges; ``locus_width`` must track the width of the locus and not the
galaxies beside it.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from shapepipe.pipeline.str_handler import StrInterpreter, stellar_locus

pytestmark = pytest.mark.decision("star_selection_psf.star_selection_box")

SETOOLS = (
    Path(__file__).resolve().parents[2]
    / "workflow/config/cfis/star_selection.setools"
)
LOCUS = 2.85  # px, stellar FWHM
LOCUS_SIGMA = 0.07  # px, median locus width on 2114045p
FWHM_MIN, FWHM_MAX = 0.3 / 0.187, 1.5 / 0.187  # the preselection bounds
SEEDS = range(20)


def mode(values):
    return StrInterpreter._mode(None, np.asarray(values, dtype=np.float32))


def width(values):
    return StrInterpreter._locus_width(
        None, np.asarray(values, dtype=np.float32)
    )


def preselection(
    seed, sigma=LOCUS_SIGMA, n_stars=70, n_galaxies=120, gap=0.4
):
    """A stellar Gaussian, a broad galaxy tail ``gap`` above it, artefacts."""
    rng = np.random.default_rng(seed)
    stars = rng.normal(LOCUS, sigma, n_stars)
    galaxies = LOCUS + gap + rng.gamma(2.0, 1.0, n_galaxies)
    small = rng.uniform(FWHM_MIN, LOCUS - 0.4, 6)
    values = np.concatenate([stars, galaxies, small])
    return values[(values > FWHM_MIN) & (values < FWHM_MAX)]


def star_box():
    """The FWHM_IMAGE lines of [MASK:star_selection], as SETools reads them."""
    text = SETOOLS.read_text().split("[MASK:star_selection]")[1]
    section = text.split("\n[")[0]
    return [
        line.replace(" ", "")
        for line in section.splitlines()
        if re.match(r"FWHM_IMAGE\s*[<>]", line)
    ]


def selected(values, preselect):
    catalogue = {"FWHM_IMAGE": np.asarray(values, dtype=np.float32)}
    keep = np.ones(len(values), dtype=bool)
    for line in star_box():
        keep &= StrInterpreter(
            line, catalogue, make_compare=True,
            mask_dict={"preselect": preselect},
        ).result
    return keep


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_recovers_the_stellar_locus(seed):
    assert abs(mode(preselection(seed)) - LOCUS) < 0.03


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_is_stable_to_objects_added_at_the_edges(seed):
    values = preselection(seed)
    edges = [FWHM_MIN + 0.01, FWHM_MAX - 0.01, FWHM_MAX - 0.02]
    shift = mode(np.concatenate([values, edges])) - mode(values)
    assert abs(shift) < LOCUS_SIGMA / 7


@pytest.mark.parametrize("seed", SEEDS)
def test_mode_is_stable_to_objects_removed_at_the_edges(seed):
    values = np.sort(preselection(seed))
    shift = mode(values[2:-1]) - mode(values)
    assert abs(shift) < LOCUS_SIGMA / 7


@pytest.mark.parametrize("sigma", [0.05, 0.08, 0.12, 0.2])
def test_locus_width_recovers_the_locus_sigma(sigma):
    """Narrower and broader than the 0.1 px kernel, with galaxies beside."""
    ratios = [width(preselection(seed, sigma)) / sigma for seed in SEEDS]
    assert abs(np.median(ratios) - 1) < 0.1
    assert np.all(np.abs(np.array(ratios) - 1) < 0.5)


def test_box_scales_with_the_locus():
    """A locus three times broader gets a box three times wider."""
    narrow = [width(preselection(seed, 0.05)) for seed in SEEDS]
    broad = [width(preselection(seed, 0.15)) for seed in SEEDS]
    assert np.median(broad) / np.median(narrow) == pytest.approx(3, rel=0.15)


@pytest.mark.parametrize("seed", SEEDS)
def test_galaxy_tail_does_not_inflate_the_width(seed):
    """Galaxies from just outside the box (3.5 sigma), 6 per star."""
    alone = width(preselection(seed, n_galaxies=0))
    crowded = width(preselection(seed, n_galaxies=400, gap=0.25))
    assert crowded / alone == pytest.approx(1, abs=0.2)


def test_a_galaxy_bump_beside_a_broad_locus_does_not_merge_into_it():
    """Poor-seeing corner CCDs: a galaxy bump 3 sigma above a broad locus."""
    widths = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        stars = rng.normal(LOCUS, 0.15, 60)
        bump = rng.normal(LOCUS + 0.45, 0.15, 45)
        widths.append(width(np.concatenate([stars, bump])))
    assert np.median(widths) < 0.25


@pytest.mark.parametrize("seed", SEEDS[:5])
def test_setools_box_is_three_locus_widths_about_the_mode(seed):
    """The committed config selects mode +- 3 * locus_width."""
    for sigma in (0.05, 0.15):
        values = preselection(seed, sigma)
        values = np.concatenate([values, np.linspace(1.0, 9.0, 400)])
        preselect = (values > FWHM_MIN) & (values < FWHM_MAX)
        centre, sig = stellar_locus(values[preselect].astype(np.float32))
        keep = selected(values, preselect)
        kept = values[keep].astype(np.float32)
        assert kept.min() >= centre - 3 * sig - 1e-4
        assert kept.max() <= centre + 3 * sig + 1e-4
        assert kept.min() < centre - 2.5 * sig
        assert kept.max() > centre + 2.5 * sig


def test_small_samples_fall_back_to_the_median_and_mad():
    values = np.array([2.7, 2.8, 2.9, 5.0])
    assert mode(values) == pytest.approx(np.median(values))
    assert width(values) == pytest.approx(
        1.4826 * np.median(np.abs(values - np.median(values))), rel=1e-3
    )


def test_empty_sample_returns_minus_one_and_selects_nothing():
    assert mode(np.array([], dtype=np.float32)) == -1
    assert width(np.array([], dtype=np.float32)) == -1
    values = np.array([2.8, 2.9])
    assert not selected(values, np.zeros(2, dtype=bool)).any()


def test_identical_values_get_the_width_floor():
    assert width(np.full(30, 2.85)) == pytest.approx(0.02)


def test_duplicated_values_do_not_break_the_estimate():
    """Bootstrap-like resamples repeat values; the peak still lands."""
    values = np.repeat(preselection(0)[:40], 5)
    assert abs(mode(values) - LOCUS) < 0.05
