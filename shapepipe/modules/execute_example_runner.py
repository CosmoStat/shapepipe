"""EXECUTE MODULE EXAMPLE.

This module defines methods for an example command line execution module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline.execute import execute


@module_runner(
    input_module="python_example_runner",
    version="1.0",
    file_pattern="pyex_output",
    file_ext=".cat",
    executes="head",
    run_method="parallel",
)
def execute_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Execute Example Runner."""
    command_line = f"head {input_file_list[0]}"
    output_file_name = (
        f'{run_dirs["output"]}/head_output{file_number_string}.txt'
    )

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, "w")
    text_file.write(stdout)

    return stdout, stderr
