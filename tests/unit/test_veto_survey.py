"""Counterfactual arm semantics, independent of the scientific bias model."""
import importlib.util
from pathlib import Path

import numpy as np

PATH = Path(__file__).resolve().parents[2] / 'scripts/validation/masked_pixels/survey_prevalence.py'


def load():
    spec = importlib.util.spec_from_file_location('survey_prevalence', PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_arm_boundaries_and_footprint_only_rejection():
    m = load()
    # An interpolated pixel outside the ISO footprint: A drops, D keeps.
    r = dict(d_interp=4., d_noise=99., hit0=False, hit2=False,
             interp_hit0=False, interp_hit2=False, known=True)
    k = m.epoch_rules(r)
    assert not k['A'] and k['B'] and k['D0'] and k['D2'] and k['E0']
    # A noise-filled hole outside the footprint: D keeps, hybrid drops.
    r.update(d_interp=99., d_noise=10.)
    k = m.epoch_rules(r)
    assert not k['A'] and k['D0'] and not k['E0']
    # Veto boundaries are strict, not inclusive.
    r.update(d_interp=7., d_noise=13.)
    assert m.epoch_rules(r)['A']
    r.update(hit2=True, interp_hit2=True)
    assert m.epoch_rules(r)['D0'] and not m.epoch_rules(r)['D2']
    assert not m.epoch_rules(r)['E2']


def test_unknown_footprint_falls_back_to_fixed_rule_not_clean():
    m = load()
    r = dict(d_interp=4., d_noise=99., hit0=False, hit2=False,
             interp_hit0=False, interp_hit2=False, known=False)
    k = m.epoch_rules(r)
    assert not k['D0'] and not k['D2'] and not k['E0']


def test_classifier_uses_production_fill_and_euclidean_dilation():
    m = load()
    flag = np.zeros((51, 51), dtype=np.int16)
    flag[25, 29] = 8
    flag[:, 35:39] = 2
    seg = np.zeros_like(flag)
    seg[25, 27] = 42
    r, components = m.classify_epoch(flag, np.ones_like(flag),
                                    np.ones_like(flag), None, seg, 42)
    assert r['d_interp'] == 4 and r['d_noise'] == 10
    assert not r['hit0'] and r['hit2'] and r['interp_hit2']
    assert r['iso_bits'] == 0
    assert set(c['fill'] for c in components) == {'interpolate', 'noise'}
