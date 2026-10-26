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


def test_multi_epoch_normalization_and_shared_calibration_error(monkeypatch):
    monkeypatch.syspath_prepend(str(PATH.parent))
    import survey_combine as m
    from survey_prevalence import EP_DTYPE, CP_DTYPE
    monkeypatch.setattr(m, 'ARMS', ('B',))
    objects = np.zeros(2, dtype=[('X_IMAGE', 'f8'), ('Y_IMAGE', 'f8')])
    objects['X_IMAGE'] = [0, 1]
    objects['Y_IMAGE'] = [0, 1]
    epochs = np.zeros(6, dtype=EP_DTYPE)
    epochs['oi'] = [0, 0, 0, 1, 1, 1]
    epochs['eligible'] = True
    epochs['ivar'] = 1.
    components = np.zeros(2, dtype=CP_DTYPE)
    components['ei'] = [0, 3]
    components['distance'] = 4
    components['fill'] = 'interpolate'
    components['kind'] = 'bleed'
    lib = [dict(kind='bleed', fill='interpolate', d=4, hlr=.3,
                psf_shear=[0., 0.], n_epochs=n) for n in (1, 3)]
    # Effective coefficient N*delta_N is half the single-epoch delta.
    q = np.array([[.3, 0, .03, 0], [.15, 0, .015, 0]])
    qb = np.zeros((2, 10, 4))
    qb[1, :, 0] = [.0, .3] * 5
    out = m.combine(epochs, objects, components, np.array([0, 0]), q, qb,
                    np.abs(q), {'B': np.ones(6, bool)}, np.ones(2, bool), lib)['B']
    np.testing.assert_allclose(out['mean'], [.05, 0, .005, 0])
    assert out['multi_corrected_components'] == 2
    # One shared fixture error does NOT shrink by sqrt(number of objects).
    np.testing.assert_allclose(out['fixture_sigma'][0], .05)
    assert out['objects_lost'] == out['epochs_lost'] == 0


def test_object_cut_reads_pre_veto_flags_not_fixed_veto_survivors(monkeypatch):
    monkeypatch.syspath_prepend(str(PATH.parent))
    import survey_combine as m
    from survey_prevalence import EP_DTYPE
    e = np.zeros(2, dtype=EP_DTYPE)
    e['eligible'] = e['known'] = True
    e['d_interp'] = [4, np.inf]
    e['d_noise'] = np.inf
    e['iso_bits'] = [8, 0]
    e['hit0'][0] = e['interp_hit0'][0] = True
    keep = m.build_keep(e, 1)
    assert keep['A'].tolist() == [False, True]
    assert keep['B'].all()
    assert not keep['C'].any()


def test_973_object_cut_does_not_invent_noise_only_zero_weight_flags(monkeypatch):
    monkeypatch.syspath_prepend(str(PATH.parent))
    import survey_combine as m
    from survey_prevalence import EP_DTYPE
    e = np.zeros(1, dtype=EP_DTYPE)
    e['eligible'] = e['known'] = e['hit0'] = True
    e['d_interp'], e['d_noise'] = np.inf, 10.
    keep = m.build_keep(e, 1)
    assert keep['C'].all()  # neither instrument ISO bits nor interpolation
    assert not keep['D0'].any()  # D sees all production defects


def test_size_matched_envelope_preserves_size_footprint_correlation(monkeypatch):
    monkeypatch.syspath_prepend(str(PATH.parent))
    import survey_combine as m
    from survey_prevalence import CP_DTYPE
    c = np.zeros(1, dtype=CP_DTYPE)
    c['fill'], c['distance'] = 'interpolate', 4
    lib = [dict(fill='interpolate', hlr=h, d=4) for h in (.3, .7)]
    q = np.array([[.01, .001, .0001, .00001], [1., .1, .01, .001]])
    envelope = m.matched_size_envelope(c, np.array([0]), lib, q)
    np.testing.assert_allclose(envelope, q[[0]])
