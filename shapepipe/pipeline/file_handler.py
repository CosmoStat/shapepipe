# -*- coding: utf-8 -*-

"""FILE HANDLER

This module defines a class for handling pipeline files.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
import re
import numpy as np
from shutil import copyfile
from glob import glob
from functools import reduce, partial
from shapepipe.pipeline.run_log import RunLog
from shapepipe.pipeline.shared import split_module_run
from shapepipe.modules.module_runners import get_module_runners
from shapepipe.utilities.file_system import mkdir


def find_files(path, pattern='*', ext='*'):
    """Find Files

    This method recursively retrieves file names from a given path that
    match a given pattern and/or have a given extension.

    Parameters
    ----------
    path : str
        Full path to files
    pattern : str, optional
        File pattern, default is '*'
    ext : str, optional
        File extension, default is '*'

    Returns
    -------
    list
        List of file names

    Raises
    ------
    ValueError
        For '*' in pattern
    ValueError
        For '*' in extension
    ValueError
        For invalid extension format

    """
    dot = '.'
    star = '*'

    if pattern != star and star in pattern:
        raise ValueError('Do not include "*" in pattern.')

    if ext != star and star in ext:
        raise ValueError('Do not include "*" in extension.')

    if (not ext.startswith(dot) and dot in ext) or (ext.count(dot) > 1):
        raise ValueError(f'Invalid extension format: "{ext}".')

    if ext != star and not ext.startswith(dot):
        ext = dot + ext

    search_string = f'{path}/**/*{pattern}*{ext}'

    return glob(search_string, recursive=True)


def check_duplicate(input_list):
    """Check Duplicate.

    Check whether input list contains at least one duplicate.

    Parameters
    ----------
    input_list : list
        input list

    Returns
    -------
    ok : bool
        True (False) if does (does not) contain at least one duplicate

    """
    input_set = set()

    for elem in input_list:
        if elem in input_set:
            return True
        else:
            input_set.add(elem)

    return False


class FileHandler(object):
    """File Handler

    This class manages the files used and produced during a pipeline run.

    Parameters
    ----------
    run_name : str
        Run name
    module_list : list
        List of modules to be run
    config : CustomParser
        Configuaration parser instance

    """

    def __init__(self, run_name, modules, config):

        self._run_name = run_name

        self._module_list = modules
        if self._module_list[-1] == '':
            raise ValueError('Invalid module list, check for a trailing comma')

        self._config = config
        self.module_runners = get_module_runners(self._module_list)
        self._input_list = config.getlist('FILE', 'INPUT_DIR')
        self._output_dir = config.getexpanded('FILE', 'OUTPUT_DIR')
        self._log_name = config.get('FILE', 'LOG_NAME')
        self._correct_pattern = config.getboolean(
            'FILE',
            'CORRECT_FILE_PATTERN',
        )
        self._run_log_file = self.setpath(
            self._output_dir,
            config.get('FILE', 'RUN_LOG_NAME'),
            '.txt',
        )
        self._module_dict = {}

        if config.has_option('FILE', 'FILE_PATTERN'):
            self._file_pattern = config.getlist('FILE', 'FILE_PATTERN')
        if config.has_option('FILE', 'FILE_EXT'):
            self._file_ext = config.getlist('FILE', 'FILE_EXT')
        if config.has_option('FILE', 'NUMBERING_SCHEME'):
            self._numbering_scheme = config.get('FILE', 'NUMBERING_SCHEME')
        else:
            self._numbering_scheme = r'RE:\_\d+'
        if config.has_option('FILE', 'NUMBER_LIST'):
            if os.path.isfile(config.get('FILE', 'NUMBER_LIST')):
                self._number_list = (
                    self.read_number_list(config.get('FILE', 'NUMBER_LIST'))
                )
            else:
                self._number_list = config.getlist('FILE', 'NUMBER_LIST')
        else:
            self._number_list = None

    @property
    def run_dir(self):
        """Run Directory

        This method defines the run directory.

        """
        return self._run_dir

    @run_dir.setter
    def run_dir(self, value):

        self._run_dir = self.check_dir(value, check_exists=True)

    @property
    def _input_dir(self):
        """Input Directory

        This method defines the input directories.

        """
        return self.__input_dir

    @_input_dir.setter
    def _input_dir(self, value):

        self.__input_dir = self.check_dirs(value)

    @property
    def _output_dir(self):
        """Output Directory

        This method defines the output directory.

        """
        return self.__output_dir

    @_output_dir.setter
    def _output_dir(self, value):

        self.__output_dir = self.check_dir(value)

    @staticmethod
    def read_number_list(file_name):
        """Read Number List

        Extract number strings to be processed from a file.

        Parameters
        ----------
        file_name : str
            Number list file name

        """
        with open(file_name) as data_file:
            number_list = data_file.readlines()

        return [value.rstrip('\n') for value in number_list]

    @classmethod
    def check_dir(cls, dir_name, check_exists=False):
        """Check Directory

        Raise error if directory exists.

        Parameters
        ----------
        dir_name : str
            Directory name

        Raises
        ------
        OSError
            If directory already exists

        """
        if check_exists and os.path.isdir(dir_name):
            raise OSError('Directory {dir_name} already exists.')

        return cls.strip_slash(dir_name)

    @classmethod
    def check_dirs(cls, dir_list):
        """Check Directories

        Check directories in list

        Parameters
        ----------
        dir_list : list
            Directory list

        """
        return [cls.check_dir(dir) for dir in dir_list]

    @classmethod
    def mkdir(cls, dir_name):
        """Make Directory

        This method creates a directory at the specified path.

        Parameters
        ----------
        dir_name : str`
            Directory name with full path

        """
        cls.check_dir(dir_name, check_exists=True)
        mkdir(dir_name)

    @staticmethod
    def setpath(path, name, ext=''):
        """Set Path Name

        This method appends the file/directory name to the input path.

        Parameters
        ----------
        path : str
            Full path
        name : str
            File or directory name
        ext : str, optional
            File extension, default is ''

        Returns
        -------
        str
            Formated path

        """
        return f'{path}/{name}{ext}'

    @staticmethod
    def strip_slash(path):
        """Strip Slash

        This method removes the trailing slash from a path.

        Parameters
        ----------
        path : str
            Full path

        Returns
        -------
        str
            Updated path

        """
        return path.rstrip('/')

    @classmethod
    def strip_slash_list(cls, path_list):
        """Strip Slash List

        This method removes the trailing slash from a list of paths.

        Parameters
        ----------
        path_list : list
            List of paths

        Returns
        -------
        list
            Updated paths

        """
        return [cls.strip_slash(path) for path in path_list]

    @staticmethod
    def flatten_list(input_list):
        """Flatten List

        Flatten a list of lists.

        Parameters
        ----------
        input_list : list
            A list of lists

        Returns
        -------
        list
            Flattened list

        """
        return [item for sublist in input_list for item in sublist]

    @staticmethod
    def _get_module_run_name(dir):
        """Get Module Run Name

        Retrieve module run name, module name and search string from input
        string.

        Parameters
        ----------
        dir : str
            Input directory string

        Returns
        -------
        tuple
            Module run name, module name, search string

        """
        string, module_run = dir.split(':')
        module = module_run.split('/')[0]

        return module_run.lower(), module.lower(), string

    def _check_input_dir_list(self, dir_list):
        """Check Input Directory List

        Check an input list to see if the directories exist or if the the run
        log should be serarched for an appropriate output directory.

        Parameters
        ----------
        dir_list : list
            List of directories

        Raises
        ------
        ValueError
            For invalid input directory value

        """
        input_dir = []

        for dir in dir_list:

            if os.path.isdir(dir):
                input_dir.append(dir)

            elif 'last' in dir.lower():
                module_run, module, _ = self._get_module_run_name(dir)
                last_module = self._run_log.get_last(module)
                input_dir.append(
                    self.setpath(
                        self.setpath(last_module, module_run),
                        'output',
                    )
                )

            elif 'all' in dir.lower():
                module_run, module, _ = self._get_module_run_name(dir)
                all_runs = self._run_log.get_all(module)
                input_dir.extend([
                    self.setpath(
                        self.setpath(run.split(' ')[0], module_run),
                        'output'
                    )
                    for run in all_runs
                ])

            elif ':' in dir.lower():
                module_run, _, string = self._get_module_run_name(dir)
                input_dir.append(
                    self.setpath(
                        self.setpath(
                            self._run_log.get_run(string),
                            module_run,
                        ),
                        'output',
                    )
                )

            else:
                raise ValueError(
                    f'Invalid INPUT_DIR ({dir}). Make sure the paths '
                    + 'provided are valid directories or use the '
                    + 'allowed special keys.'
                )

        return input_dir

    def _get_input_dir(self):
        """Get Input Directory

        This method sets the module input directory

        """
        self._input_dir = self._check_input_dir_list(self._input_list)

    def create_global_run_dirs(self):
        """Create Global Run Directories

        This method creates the pipeline output directories for a given run.

        """
        self.run_dir = self.setpath(self._output_dir, self._run_name)
        self._log_dir = self.setpath(self.run_dir, 'logs')
        self._tmp_dir = self.setpath(self.run_dir, 'tmp')
        self.log_name = self.setpath(self._log_dir, self._log_name)
        self._run_log = RunLog(
            self._run_log_file,
            self._module_list,
            self.run_dir,
        )

        self.mkdir(self.run_dir)
        self.mkdir(self._log_dir)
        self.mkdir(self._tmp_dir)

        self._get_input_dir()
        self._copy_config_to_log()

    def _copy_config_to_log(self):
        """Copy Config to Log

        Copy configuration file to run log directory.

        """
        config_file_name = os.path.basename(self._config.file_name)

        copyfile(
            self._config.file_name,
            f'{self._log_dir}/{config_file_name}',
        )

    def get_module_current_run(self, module):
        """Get Module Current Run

        Get the run current run count for the module.

        Parameters
        ----------
        module : str
            Module name

        Returns
        -------
        str
            Module run count

        """
        return str(self._module_dict[module]['run_count'])

    def get_module_config_sec(self, module):
        """Get Module Configuration Section

        Get the name of section name in the configuration file for the module.

        Parameters
        ----------
        module : str
            Module name

        Returns
        -------
        str
            Configuration file section name

        """
        return self._module_dict[module]['latest'].upper()

    def get_add_module_property(self, run_name, property):
        """Get Additional Module Properties

        Get a list of additional module property values.

        Parameters
        ----------
        run_name : str
            Module run name
        property : str
            Property name

        Returns
        -------
        list
            Additional module property values

        """
        if (self._config.has_option(
            run_name.upper(),
            f'ADD_{property.upper()}',
        )):

            return self._config.getlist(
                run_name.upper(),
                f'ADD_{property.upper()}',
            )

    def _set_module_property(self, module, run_name, property, get_type):
        """Set Module Property

        Set a module property from either the configuration file or the module
        runner.

        Parameters
        ----------
        module : str
            Module name
        property : str
            Property name
        get_type : str
            Type of object to get from config file

        Notes
        -----
        Module properties are set with the following hierarchy:

        1. Look for properties in module section of the config file,
        2. look for default properties in the file handler, which sources the
        ``[FILE]`` section of the config file,
        3. look for default properties in the module runner definition.

        In other words, module runner definitions will only be used if the
        properties are not set in the confg file. File handler properties will
        be used for all modules, overriding all module runner definitions.
        Module specific properties from the config file will override all other
        definitions, but only for the module in question.

        """
        # 1) Check for parameter value in module section of config file
        if self._config.has_option(run_name.upper(), property.upper()):
            if get_type == 'str':
                prop_val = self._config.get(
                    run_name.upper(),
                    property.upper(),
                )
            elif get_type == 'list':
                prop_val = self._config.getlist(
                    run_name.upper(),
                    property.upper(),
                )
            else:
                raise ValueError(f'{get_type} is not a valid get type')

        # 2) Check for default parameter values in file handler
        elif hasattr(self, f'_{property}'):
            prop_val = getattr(self, f'_{property}')

        # 3) Check for default parameter values in module runner
        elif hasattr(self.module_runners[module], property):
            prop_val = getattr(self.module_runners[module], property)

        else:
            raise ValueError(
                f'No value for {property} in {module} could be found.'
            )

        # Look for additional module properties for list objects
        if (
            isinstance(prop_val, list)
            and self.get_add_module_property(run_name, property)
        ):
            prop_val += self.get_add_module_property(run_name, property)

        self._module_dict[module][run_name][property] = prop_val

    def _set_module_properties(self, module, run_name):
        """Get Module Properties

        Get module properties defined in module runner wrapper.

        Parameters
        ----------
        module : str
            Module name

        """
        module_props = {
            'numbering_scheme': 'str',
            'run_method': 'str',
            'input_module': 'list',
            'file_pattern': 'list',
            'file_ext': 'list',
            'depends': 'list',
            'executes': 'list'
        }

        for property, get_type in module_props.items():
            self._set_module_property(module, run_name, property, get_type)

        # Make sure the number of patterns and extensions match
        if (
            (len(self._module_dict[module][run_name]['file_ext']) == 1)
            and (len(self._module_dict[module][run_name]['file_pattern']) > 1)
        ):
            self._module_dict[module][run_name]['file_ext'] = [
                self._module_dict[module][run_name]['file_ext'][0]
                for i in self._module_dict[module][run_name]['file_pattern']
            ]

        if (
            len(self._module_dict[module][run_name]['file_ext'])
            != len(self._module_dict[module][run_name]['file_pattern'])
        ):
            n_fext = len(self._module_dict[module][run_name]['file_ext'])
            n_fpat = len(self._module_dict[module][run_name]['file_pattern'])
            raise ValueError(
                f'The number of file_ext values ({n_fext}) does not match the '
                + f'number of file_pattern values ({n_fpat}).'
            )

    def _create_module_run_dirs(self, module, run_name):
        """Create Module Run Directories

        This method creates the module output directories for a given run.

        Parameters
        ----------
        module : str
            Module name

        """
        self._module_dict[module][run_name]['run_dir'] = (
            self.setpath(self._run_dir, run_name)
        )
        self._module_dict[module][run_name]['log_dir'] = (
            self.setpath(
                self._module_dict[module][run_name]['run_dir'],
                'logs',
            )
        )
        self._module_dict[module][run_name]['output_dir'] = (
            self.setpath(
                self._module_dict[module][run_name]['run_dir'],
                'output',
            )
        )

        self.mkdir(self._module_dict[module][run_name]['run_dir'])
        self.mkdir(self._module_dict[module][run_name]['log_dir'])
        self.mkdir(self._module_dict[module][run_name]['output_dir'])

        # Set current output directory to module output directory
        self.output_dir = self._module_dict[module][run_name]['output_dir']
        self.module_run_dirs = {
            'run': self.run_dir,
            'log': self._log_dir,
            'tmp': self._tmp_dir,
            'output': self.output_dir
        }

    def _set_module_input_dir(self, module, run_name):
        """Set Module Input Directory

        Set the module input directory. If the module specified is the
        first module in the pipeline or does not have any input modules then
        only the `INPUT_DIR` from `[FILE]` is used, otherwise the output
        directory from the preceding module is used.

        Additional input directories can be specified with `INPUT_DIR` from
        `[MODULE]`.

        Parameters
        ----------
        module : str
            Module name

        """
        if (
            isinstance(
                self._module_dict[module][run_name]['input_module'],
                type(None),
            )
            or len(self._module_dict) == 1
        ):
            input_dir = self._input_dir

        else:

            print(self._module_dict)

            input_dir = []
            for input_module in (
                self._module_dict[module][run_name]['input_module']
            ):
                input_module, in_mod_run = split_module_run(input_module)
                if in_mod_run == input_module:
                    in_mod_run = self._module_dict[input_module]['latest']
                if input_module in self._module_dict:
                    input_dir.append(
                        self._module_dict[input_module][in_mod_run][
                            'output_dir'
                        ]
                    )

        if self._config.has_option(run_name.upper(), 'INPUT_DIR'):
            input_dir = self._check_input_dir_list(
                self._config.getlist(run_name.upper(), 'INPUT_DIR')
            )

        if self.get_add_module_property(run_name, 'input_dir'):
            input_dir += self.get_add_module_property(run_name, 'input_dir')

        if not input_dir:
            raise RuntimeError(
                'Could not find appropriate input directory '
                + f'for module {run_name}.'
            )

        self._module_dict[module][run_name]['input_dir'] = (
            self.check_dirs(input_dir)
        )

    @staticmethod
    def _generate_re_pattern(match_pattern):
        """Generate Regular Expression Pattern

        Generate a regular expression pattern from an input string.

        Parameters
        ----------
        match_pattern : str
            Pattern string

        Returns
        -------
        _sre.SRE_Pattern
            Regular expression pattern

        Raises
        ------
        TypeError
            For invalid input type

        """
        if not isinstance(match_pattern, str):
            TypeError('Match pattern must be a string.')

        chars = [char for char in match_pattern if not char.isalnum()]
        split_pattern = '|'.join(chars).replace('.', r'\.')
        chars = [f'\\{char}' for char in chars] + ['']
        num_length = [
            f'\\d{{{len(digits)}}}'
            for digits in re.split(split_pattern, match_pattern)
        ]
        re_pattern = r''.join(
            [a for b in zip(num_length, chars) for a in b]
        ).replace('{1}', '+')

        return re.compile(re_pattern)

    @staticmethod
    def _strip_dir_from_file(file_name, dir_list):
        """Strip Directory from File Name

        Remove the directory string from the file name.

        Parameters
        ----------
        file_name : str
            File name
        dir_list : list
            Input directory list

        Returns
        -------
        str
            File name

        """
        return [
            file_name.replace(_dir + '/', '')
            for _dir in dir_list if _dir in file_name
        ][0]

    @classmethod
    def _get_re(cls, num_scheme):
        """Get Regular Expression

        Return the regular expression corresponding to the numbering scheme.

        Parameters
        ----------
        num_scheme : str
            Numbering scheme

        Raises
        ------
        ValueError
            if num_scheme is None

        Returns
        -------
        str
            Regular Expression

        """
        # Raise an error if num_scheme is None,
        # but not if num_scheme==''
        if not num_scheme and num_scheme != '':

            raise ValueError('No numbering scheme adapted')

        elif num_scheme.startswith('RE:'):

            re_pattern = num_scheme.replace('RE:', '')

        else:

            re_pattern = cls._generate_re_pattern(num_scheme)

        return re_pattern

    def _save_num_patterns(
        self,
        dir_list,
        re_pattern,
        pattern,
        ext,
        output_file
    ):
        """Save Number Patterns

        Save file number patterns to numpy binary, update file patterns and
        get correct file paths.

        Parameters
        ----------
        dir_list : list
            List of input directories
        re_pattern : str
            Regular expression pattern
        pattern : str
            File pattern
        ext : str
            File extension
        output_file : str
            Output file name

        """
        # Find all files matching the input pattern and extension from the
        # available input directories and identify the correct path
        true_file_list = None
        true_path = None

        for path in dir_list:

            file_list = find_files(path, pattern, ext)

            if file_list:
                true_file_list = file_list
                true_path = path
                del file_list
                break

        if not true_file_list:
            raise RuntimeError(
                f'No files found matching "{pattern}" and "{ext}" in the '
                + f'directories {dir_list}.'
            )

        # Correct the extension if necessary
        new_ext = '.' + ext if not ext.startswith('.') else ext

        if new_ext != ext:
            print(f'Updating file extension from "{ext}" to "{new_ext}".')
            print()

        # Select files matching the numbering scheme
        final_file_list = []
        found_match = False
        correct_pattern = self._correct_pattern
        new_pattern = pattern

        for file in true_file_list:

            striped = self._strip_dir_from_file(file, dir_list)
            search_res = re.search(re_pattern, striped)

            if search_res:
                file_name = search_res.group()
                final_file_list.append(file_name)
                found_match = True

            # Correct the pattern if necessary
            if found_match and correct_pattern:

                new_pattern = striped

                for substring in (new_ext, file_name):
                    new_pattern = new_pattern.replace(substring, '')

                if new_pattern != pattern:
                    print(
                        f'Updating pattern from "{pattern}" to '
                        + f'"{new_pattern}".'
                    )
                    print()

                correct_pattern = False

        if not found_match:
            raise RuntimeError(
                f'Could not match numbering scheme "{self._numbering_scheme}" '
                + f'to any of the input files matching "{pattern}" and '
                + f'"{ext}" in the directories {dir_list}.'
            )

        if check_duplicate(final_file_list):
            raise RuntimeError(
                'Input file list contains at least two elements that match '
                + 'file pattern and numbering scheme, leading to identical '
                + 'input files.  Make sure that the correct input '
                + 'directory is used.'
            )

        # Save file list
        np.save(output_file, np.array(final_file_list))

        del true_file_list, final_file_list

        return new_pattern, new_ext, true_path

    @staticmethod
    def _save_match_patterns(output_file, mmap_list):
        """Save Match Patterns

        Save matching number patterns to numpy binary.

        Parameters
        ----------
        output_file : str
            Output file name
        mmap_list : list
            List of memory maps

        """
        num_pattern_list = [np.load(mmap, mmap_mode='r') for mmap in mmap_list]

        np.save(
            output_file,
            reduce(
                partial(np.intersect1d, assume_unique=True),
                num_pattern_list
            )
        )

        del num_pattern_list

    @staticmethod
    def _get_file_name(path, pattern, number, ext):
        """Get File Name

        Get file name corresponding to the path, file pattern, number pattern
        and file extension.

        Parameters
        ----------
        path : str
            Path to file
        pattern : str
            File pattern
        number : str
            Number pattern
        ext : str
            File extension

        Retunrs
        -------
        str
            File name

        """

        return f'{path}/{pattern}{number}{ext}'

    @staticmethod
    def _remove_mmaps(mmap_list):
        """Remove Memory Maps

        Remove memory map files in input list.

        Parameters
        ----------
        mmap_list : list or str
            List of memory map files

        """

        if isinstance(mmap_list, str):
            mmap_list = [str]

        for mmap in set(mmap_list):
            os.remove(mmap)

    def _format_process_list(
        self,
        patterns,
        memory_map,
        re_pattern,
        num_scheme,
        run_method,
    ):
        """Format Process List

        Format the list of files to be processed.

        Parameters
        ----------
        patterns : list
            List of file patterns
        memory_map : str
            Name of memory map file
        re_pattern : str
            Regular expression for numbering scheme
        num_scheme : str
            Numbering scheme
        run_method : str
            Run method

        Returns
        -------
        list
            List of processes

        """
        pattern_list, ext_list, path_list = list(zip(*patterns))

        if isinstance(self._number_list, type(None)):
            number_list = np.load(memory_map, mmap_mode='r')
        else:
            number_list = self._number_list

        process_list = []

        for number in number_list:

            if not re.search(re_pattern, number):
                raise ValueError(
                    f'The string "{number}" does not match the '
                    + f'numbering scheme "{num_scheme}".'
                )

            if run_method == 'serial':
                process_items = []
            else:
                process_items = [number]
            process_items.extend([
                self._get_file_name(path, fp, number, ext)
                for path, fp, ext in zip(path_list, pattern_list, ext_list)
            ])
            process_list.append(process_items)

        return process_list

    def _save_process_list(
        self,
        dir_list,
        pattern_list,
        ext_list,
        num_scheme,
        run_method,
    ):
        """Save Process List

        Save list of processes to a numpy binary.

        Parameters
        ----------
        dir_list : list
            List of input directories
        pattern_list : list
            List of file patterns
        ext_list : list
            List of file extensions
        num_scheme : str
            Numbering scheme

        """
        np_mmap_list = [
            self.setpath(self._tmp_dir, f'nums_{pattern}_{ext}.npy')
            for pattern, ext in zip(pattern_list, ext_list)
        ]
        match_mmap = self.setpath(self._tmp_dir, 'matching_num_patterns.npy')
        self.process_mmap = self.setpath(self._tmp_dir, 'process_list.npy')

        re_pattern = self._get_re(num_scheme)

        temp = [
            self._save_num_patterns(
                dir_list,
                re_pattern,
                pattern,
                ext,
                np_mmap
            )
            for pattern, ext, np_mmap in
            zip(pattern_list, ext_list, np_mmap_list)
        ]

        self._save_match_patterns(match_mmap, np_mmap_list)

        process_list = self._format_process_list(
            temp,
            match_mmap,
            re_pattern,
            num_scheme,
            run_method,
        )

        np.save(self.process_mmap, np.array(process_list))
        del process_list

        self.process_list = np.load(self.process_mmap, mmap_mode='r')

        self.missed = []

        self._remove_mmaps(np_mmap_list + [match_mmap])

    def remove_process_mmap(self):
        """Remove Process MMAP

        Remove process list memory map.

        """
        self._remove_mmaps([self.process_mmap])

    def _get_module_input_files(self, module, run_name):
        """Get Module Input Files

        Retrieve the module input files names from the input directory.

        Parameters
        ----------
        module : str
            Module name

        """
        dir_list = self._module_dict[module][run_name]['input_dir']
        pattern_list = self._module_dict[module][run_name]['file_pattern']
        ext_list = self._module_dict[module][run_name]['file_ext']
        num_scheme = self._module_dict[module][run_name]['numbering_scheme']
        run_method = self._module_dict[module][run_name]['run_method']

        self._save_process_list(
            dir_list,
            pattern_list,
            ext_list,
            num_scheme,
            run_method,
        )

    def set_up_module(self, module):
        """Set Up Module

        Set up module parameters for file handler.

        Parameters
        ----------
        module : str
            Module name

        """

        multi_call = self._module_list.count(module) > 1

        if module in self._module_dict.keys():
            self._module_dict[module]['run_count'] += 1
        else:
            self._module_dict[module] = {}
            self._module_dict[module]['run_count'] = 1

        if multi_call:
            call_num = self._module_dict[module]['run_count']
            run_name = f'{module}_run_{call_num}'
        else:
            run_name = module

        self._module_dict[module]['latest'] = run_name
        self._module_dict[module][run_name] = {}

        self._set_module_properties(module, run_name)
        self._create_module_run_dirs(module, run_name)
        self._set_module_input_dir(module, run_name)
        self._get_module_input_files(module, run_name)

    def get_worker_log_name(self, module, file_number_string):
        """Get Worker Log Name

        This method generates a worker log name.

        Parameters
        ----------
        module : str
            Module name
        file_number_string : str
            File numbering in output

        Returns
        -------
        str
            Worker log file name

        """
        run_name = self._module_dict[module]['latest']
        log_dir = self._module_dict[module][run_name]['log_dir']

        return f'{log_dir}/process{file_number_string}'
