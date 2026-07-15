"""FAKE PSF RUNNER.

Module runner for ``fake_psf``.

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
def fake_psf_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Fake PSF Runner."""
    sexcat_path = input_file_list[0]
    psf_dict_path = config.getexpanded(module_config_sec, "PSF_DICT_PATH")
    output_path = f'{run_dirs["output"]}/galaxy_psf{file_number_string}.sqlite'

    inst = fake_psf.FakePsf(sexcat_path, psf_dict_path, output_path, w_log)
    inst.process()

    return None, None
