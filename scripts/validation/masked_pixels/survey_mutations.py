"""In-memory negative checks for survey transport; run in the SIF allocation."""
import inspect
from pathlib import Path
import subprocess
import sys

MUTATIONS = ('ignore-joint-fit-dilution', 'average-away-shared-fixture-error', 'cut-flags-after-fixed-veto')


def run(mutation):
    import pytest
    import survey_combine as module
    target = 'test_multi_epoch_normalization_and_shared_calibration_error'
    function = module.build_keep if mutation == 'cut-flags-after-fixed-veto' else module.combine
    source = inspect.getsource(function)
    if mutation == 'cut-flags-after-fixed-veto':
        old, new = "flag = e['eligible'] &", "flag = masks['A'] &"
        target = 'test_object_cut_reads_pre_veto_flags_not_fixed_veto_survivors'
    elif mutation == 'ignore-joint-fit-dilution':
        old, new = 'if use_multi:', 'if False:'
    else:
        old = "weights=alpha[c['ei']] / alive.sum()"
        new = "weights=alpha[c['ei']] / alive.sum() ** 1.5"
    assert source.count(old) == 1
    exec(source.replace(old, new), module.__dict__)
    class Tracker:
        failures = 0

        def pytest_runtest_logreport(self, report):
            if report.when == 'call' and report.failed:
                self.failures += 1

    tracker = Tracker()
    rc = pytest.main(['tests/unit/test_veto_survey.py::' + target,
                      '-o', 'addopts=', '-q', '-p', 'no:cacheprovider'], plugins=[tracker])
    assert rc == 1 and tracker.failures == 1, (mutation, rc, tracker.failures)
    print('KILLED', mutation)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] not in MUTATIONS:
            raise SystemExit(MUTATIONS)
        run(sys.argv[1])
    else:
        for mutation in MUTATIONS:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), mutation], check=True)
