"""FAKE INTERP RUNNER.

Module runner for ``fake_psf`` under the ``<psf>_interp_runner`` name.

The Snakemake workflow's configs read the galaxy PSF from
``${SP_PSF}_interp_runner`` for every PSF model, so with ``psf_model: fake``
(image simulations, true PSF) this runner stands where ``psfex_interp_runner``
and ``mccd_interp_runner`` stand for the real data. It writes the same
``galaxy_psf`` SqliteDict from the simulation's true PSF: its PSFEx models
(``PSFEX_MODELS_DIR``), evaluated at each object's position on each epoch by
``psfex_interp``'s multi-epoch interpolation without its gate, or else its PSF
dictionary (``PSF_DICT_PATH``), one stamp per exposure and CCD.
``fake_psf_runner`` is the dictionary path under its original name, used by
the legacy bash job scripts.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.fake_psf_package import fake_psf
from shapepipe.modules.psfex_interp_package import psfex_interp
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    file_pattern=["sexcat", "log_exp_headers"],
    file_ext=[".fits", ".sqlite"],
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
    sexcat_path, headers_path = input_file_list[:2]
    psfex_models = config.getexpanded(module_config_sec, "PSFEX_MODELS_DIR")
    if psfex_models:
        inst = psfex_interp.PSFExInterpolator(
            None,
            sexcat_path,
            run_dirs["output"],
            file_number_string,
            w_log,
            config.getlist(module_config_sec, "POSITION_PARAMS"),
            True,
            0,
            float("inf"),
        )
        inst.process_me(
            [psfex_models],
            config.get(module_config_sec, "ME_DOT_PSF_PATTERN"),
            headers_path,
        )
        return None, None

    psf_dict_path = config.getexpanded(module_config_sec, "PSF_DICT_PATH")
    output_path = f'{run_dirs["output"]}/galaxy_psf{file_number_string}.sqlite'
    inst = fake_psf.FakePsf(sexcat_path, psf_dict_path, output_path, w_log)
    inst.process()
    return None, None
