
"""PSFEX PACKAGE
This module contains a wrapper class to prepare the psfex command line.
:Author: Axel Guinot
"""
import os
import re

class PSFex_caller:
    """This class contains functions to generate a PSFex command line.
    Parameters:
    ----------
    psfex_executable_path: str
        Full path to the PSFEx executable
    input_file_path: str
        Full path of the input file
    psfex_config_file: str
        Full path of the PSFEx configuration file (.psfex)
    output_dir: str
        Full path of the output directory
    outcatalog_name: str
        Full path pf the output catalogue
    check_image_list: list of str
        List of check images

    Returns:
    command_line: str
        Command line with correct and complete paths for PSFEx execution.
    """


    def __init__(self,
        psfex_executable_path,
        input_file_path,
        psfex_config_file,
        output_dir,
        outcatalog_name,
        check_image_list):
        self.psfex_executable_path = psfex_executable_path
        self.input_file_path = input_file_path
        self.psfex_config_file = psfex_config_file
        self.output_dir = output_dir
        self.outcatalog_name = outcatalog_name
        self.check_image_list = check_image_list

    def generate_command(self):
        """ Generates a command line for running PSFEx """


        #prepare command line
        command_line = (f'{self.psfex_executable_path} '
                        + f'{self.input_file_path} '
                        + f'-c {self.psfex_config_file} '
                        + f'-PSF_DIR {self.output_dir} '
                        + f'-OUTCAT_NAME {self.outcatalog_name}')

        if ((len(self.check_image_list) == 1)
             & (self.check_image_list[0] == '')) :
            check_type_list = ['NONE']
            check_name_list = ['none']
        else:
            #get pattern for filenaming from file in list
            input_file_name = os.path.split(input_file_path)[-1]
            input_file_noext = os.path.splitext(input_file_name)[0]
            suffix = re.split(file_number_string,input_file_noext)[0]

            check_type_list = []
            check_name_list = []
            for check_image in self.check_image_list:
                check_type_list.append(check_image.upper())
                check_name_list.append(f'{output_dir}'
                                + f'/{suffix}'
                                + f'{check_image.lower()}'
                                + f'{file_number_string}.fits')

        #add checks to command line
        command_line += (f' -CHECKIMAGE_TYPE {",".join(check_type_list)}'
                       + f' -CHECKIMAGE_NAME {",".join(check_name_list)}')

        self.command_line = command_line
        return command_line

    def exec(self):
        """
        Executes the command line.
        """
        self.generate_command()
        stderr, stdout = execute(self.command_line)

        return stderr, stdout

    def parse_errors(stderr, stdout):
        """
        Move errors from output from psfex to errors.
        """
        check_error = re.findall('error', stdout.lower())
        check_error2 = re.findall('all done', stdout.lower())

        if check_error == []:
            stderr2 = ''
        else:
            stderr2 = stdout

        if check_error2 == []:
            stderr2 = stdout

        return stdout, stderr2
