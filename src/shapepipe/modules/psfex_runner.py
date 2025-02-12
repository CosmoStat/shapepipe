"""PSFEX RUNNER.

Module runner for ``psfex``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.psfex_package.psfex_script import PSFExCaller
from shapepipe.pipeline.execute import execute


@module_runner(
    input_module="setools_runner",
    version="1.0",
    file_pattern=["star_selection"],
    file_ext=[".fits"],
    executes="psfex",
)
def psfex_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The PSFEx Runner."""
    # Extract psfex run configurations
    if config.has_option(module_config_sec, "EXEC_PATH"):
        psfex_executable_path = config.getexpanded(
            module_config_sec, "EXEC_PATH"
        )
    else:
        psfex_executable_path = "psfex"
    output_dir = run_dirs["output"]

    outcatalog_name = f"{output_dir}/psfex_cat{file_number_string}.cat"

    psfex_config_file = config.getexpanded(module_config_sec, "DOT_PSFEX_FILE")

    input_file_path = input_file_list[0]

    # Check image options
    if config.has_option(module_config_sec, "CHECKIMAGE"):
        check_image_list = config.getlist(module_config_sec, "CHECKIMAGE")
    else:
        check_image_list = [""]

    # Create psfex caller class instance
    psfex_inst = PSFExCaller(
        psfex_executable_path,
        input_file_path,
        psfex_config_file,
        output_dir,
        outcatalog_name,
        check_image_list,
    )

    # Generate psfex command line
    command_line = psfex_inst.generate_command()
    w_log.info(f"Running command '{command_line}'")

    # Execute command line
    stderr, stdout = execute(command_line)

    # Parse psfex errors
    stdout, stderr = PSFExCaller.parse_errors(stderr, stdout)

    # Return stdout and stderr
    return stdout, stderr
