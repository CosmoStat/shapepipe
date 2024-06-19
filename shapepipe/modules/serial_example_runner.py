"""SERIAL MODULE EXAMPLE.

This module defines methods for an example serial module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.modules.module_decorator import module_runner


class Dummy(object):
    """Dummy Class."""

    def __init__(self):

        pass

    def _read_file(self, file_name):
        """Read File."""
        with open(file_name) as data_file:
            content = data_file.read().replace("\n", "")

        return content

    def read_files(self, input_file_list):
        """Read Files."""
        self.content = ""

        for file_list in input_file_list:
            for file_name in file_list:
                self.content += self._read_file(file_name)

    def write_file(self, file_name):
        """Write Files."""
        text_file = open(file_name, "w")
        text_file.write(self.content)
        text_file.close()


@module_runner(
    version="1.1",
    input_module="python_example_runner",
    file_pattern=["numbers", "letters", "pyex_output"],
    file_ext=[".txt", ".txt", ".cat"],
    depends="numpy",
    run_method="serial",
)
def serial_example_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The Serial Example Runner."""
    output_file_name = (
        f'{run_dirs["output"]}/serial_outputfile_number_string.cat'
    )

    inst = Dummy()
    inst.read_files(input_file_list)
    inst.write_file(output_file_name)

    return inst.content, None
