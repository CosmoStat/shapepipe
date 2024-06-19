"""GET IMAGES RUNNER.

Module runner for ``get_images``

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.modules.get_images_package.get_images import GetImages
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    version="1.1",
    depends=["numpy"],
    run_method="serial",
    numbering_scheme="_0",
)
def get_images_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Get Images Runner."""
    # Read config file section

    # Copy/download method
    retrieve_method = config.get(module_config_sec, "RETRIEVE")
    retrieve_ok = ["vos", "symlink"]
    if retrieve_method not in retrieve_ok:
        raise ValueError(
            f"key RETRIEVE={retrieve_method} is invalid, "
            + f"must be in {retrieve_ok}"
        )

    if config.has_option(module_config_sec, "RETRIEVE_OPTIONS"):
        retrieve_options = config.getexpanded(
            module_config_sec, "RETRIEVE_OPTIONS"
        )
    else:
        retrieve_options = None

    if config.has_option(module_config_sec, "N_TRY"):
        n_try = config.getint(module_config_sec, "N_TRY")
    else:
        n_try = 3

    # Paths
    input_dir = config.getlist(module_config_sec, "INPUT_PATH")
    input_file_pattern = config.getlist(module_config_sec, "INPUT_FILE_PATTERN")
    input_file_ext = config.getlist(module_config_sec, "INPUT_FILE_EXT")
    output_file_pattern = config.getlist(
        module_config_sec, "OUTPUT_FILE_PATTERN"
    )

    input_numbering = config.get(module_config_sec, "INPUT_NUMBERING")

    if config.has_option(module_config_sec, "OUTPUT_PATH"):
        output_dir = config.getexpanded(module_config_sec, "OUTPUT_PATH")
    else:
        output_dir = run_dirs["output"]

    # Flags to check for already retrieved files
    if config.has_option(module_config_sec, "CHECK_EXISTING_DIR"):
        check_existing_dir = config.getexpanded(
            module_config_sec, "CHECK_EXISTING_DIR"
        )
        if config.has_option(module_config_sec, "N_EXPECTED"):
            n_expected = config.getint(module_config_sec, "N_EXPECTED")
        else:
            n_expected = 1
    else:
        check_existing_dir = None
        n_expected = None

    # Create get images class instance
    get_images_inst = GetImages(
        retrieve_method,
        retrieve_options,
        input_file_list,
        input_numbering,
        input_file_pattern,
        input_file_ext,
        output_file_pattern,
        w_log,
        check_existing_dir=check_existing_dir,
        n_expected=n_expected,
        n_try=n_try,
    )

    # Run processing
    get_images_inst.process(input_dir, output_dir)

    # No return objects
    return None, None
