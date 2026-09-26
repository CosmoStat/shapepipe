"""FAKE INTERP RUNNER.

Module runner for ``fake_psf`` under the ``<psf>_interp_runner`` name.

The Snakemake workflow's configs read the galaxy PSF from
``${SP_PSF}_interp_runner`` for every PSF model, so with ``psf_model: fake``
(image simulations, true PSF) this runner stands where ``psfex_interp_runner``
and ``mccd_interp_runner`` stand for the real data. It writes the same
``galaxy_psf`` SqliteDict, taken from the simulation's PSF dictionary instead
of a fitted model. ``fake_psf_runner`` is the same module under its original
name, used by the legacy bash job scripts.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.fake_psf_package import fake_psf
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.0",
    file_pattern=["sexcat"],
    file_ext=".fits",
    depends=["numpy", "astropy", "sqlitedict"],
    numbering_scheme="-000-000",
)
def fake_interp_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Fake Interp Runner."""
    sexcat_path = input_file_list[0]
    psf_dict_path = config.getexpanded(module_config_sec, "PSF_DICT_PATH")
    output_path = f'{run_dirs["output"]}/galaxy_psf{file_number_string}.sqlite'

    inst = fake_psf.FakePsf(sexcat_path, psf_dict_path, output_path, w_log)
    inst.process()

    return None, None
