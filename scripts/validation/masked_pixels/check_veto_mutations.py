"""Negative checks for the Napoli veto decision; run in an allocated SIF.

Each mutation must produce exactly one assertion failure.
No tracked file is modified; each mutation gets a fresh pytest process.
"""
import subprocess
import sys

MUTATIONS = ("remove-interpolation-veto", "remove-noise-veto",
             "ignore-additive-bias", "noise-as-interpolation")


def check(mutation):
    import pytest
    from shapepipe.modules.ngmix_package import ngmix as ngm

    class Injector:
        failures = 0

        def pytest_collection_modifyitems(self, items):
            # Mutate after collection so tests keep their original geometry
            # constants; only production selection changes.
            if mutation == "remove-interpolation-veto":
                ngm.EPOCH_INTERPOLATED_DEFECT_RADIUS = 0
                return
            if mutation == "remove-noise-veto":
                ngm.EPOCH_CENTRAL_DEFECT_RADIUS = 0
                return
            for module in {item.module for item in items}:
                original = module.defect_response

                def response(*args, _original=original, **kwargs):
                    if mutation == "noise-as-interpolation":
                        kwargs["options"] = {"defect_fill": "interpolate"}
                    result = _original(*args, **kwargs)
                    if mutation == "ignore-additive-bias":
                        result["c"] = [0., 0.]
                    return result

                module.defect_response = response

        def pytest_runtest_logreport(self, report):
            if report.when == "call" and report.failed:
                self.failures += 1

    if mutation == "remove-interpolation-veto":
        target = ("tests/module/test_ngmix_defect_fill.py::"
                  "test_a_flagged_column_near_the_centre_is_still_vetoed")
    elif mutation == "remove-noise-veto":
        target = ("tests/module/test_ngmix_defect_fill.py::"
                  "test_under_noise_fill_every_defect_is_vetoed_at_the_noise_radius")
    else:
        target = ("tests/science/test_defect_recovery.py::"
                  "test_inside_radius_bias_is_significant")
    plugin = Injector()
    rc = pytest.main([target, "-o", "addopts=", "-q", "-p", "no:cacheprovider"],
                     plugins=[plugin])
    assert rc == 1 and plugin.failures == 1, (mutation, rc, plugin.failures)
    print(f"KILLED {mutation}: exactly one assertion failure", flush=True)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        if sys.argv[1] not in MUTATIONS:
            raise SystemExit(f"Expected one of {MUTATIONS}")
        check(sys.argv[1])
    else:
        for mutation in MUTATIONS:
            subprocess.run([sys.executable, __file__, mutation], check=True)
